"""Re-relax interface slabs with cell freedom (paper-style, Tatsumi 2026 Eq. 4 recipe).

Mask: relax in-plane (xx, yy, xy=gamma); fix out-of-plane (zz, yz, xz) to preserve vacuum.

Restart from prior atom-only-relaxed CIFs in relax/{tag}_interface_PBE.cif
(these had cell pinned at alpha's lattice; now the cell will adjust to a
compromise between alpha and beta natural lattices).

Then compute paper-style W_ad using existing E_alpha-slab (at alpha natural) and
E_beta-slab,natural (from script 19).

Output:
  relax_cell/{tag}_interface_cellrelax_PBE.{cif,json}
  results/wad_paper_full.json (final compilation)
"""
import json, time
from pathlib import Path
import numpy as np
from ase.io import read, write
from ase.optimize import BFGS
from ase.filters import FrechetCellFilter
from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator
from pfp_api_client.pfp.estimator import Estimator, EstimatorCalcMode

EV_TO_J = 1.602176634e-19; ANG2_TO_M2 = 1e-20

WORK = Path("/home/jovyan/work_dir/cosn3_betasn_adhesion")
RELAX = WORK / "relax"             # prior atom-only relaxed
RELAX_CR = WORK / "relax_cell"     # new cell-relaxed
RELAX_CR.mkdir(exist_ok=True)
RES = WORK / "results"

calc = ASECalculator(Estimator(calc_mode=EstimatorCalcMode.PBE))

# Mask order in FrechetCellFilter: 6 components in Voigt notation
# [eps_xx, eps_yy, eps_zz, eps_yz, eps_xz, eps_xy]
# True = strain component allowed to vary, False = fixed
# We relax: xx, yy (in-plane lengths), xy (gamma between in-plane vectors)
# We fix: zz (vacuum direction length), yz, xz (no out-of-plane mixing)
MASK = [True, True, False, False, False, True]

idx = json.loads((WORK / "slab_split_index.json").read_text())

print("=" * 130)
print(f"{'pair':22s} {'term':>4s} {'N':>5s} "
      f"{'a_old':>7s} {'b_old':>7s} {'γ_old':>6s} "
      f"{'a_new':>7s} {'b_new':>7s} {'γ_new':>6s} "
      f"{'A_old':>7s} {'A_new':>7s} {'ΔA%':>6s} "
      f"{'E_old':>10s} {'E_new':>10s} {'ΔE [eV]':>9s} {'steps':>5s} {'conv':>5s} {'t[s]':>5s}")
print("=" * 130)

records = {}
for key, info in idx.items():
    cif_in = RELAX / f"{key}_interface_PBE.cif"
    atoms = read(cif_in)
    n = len(atoms)
    a0, b0, c0 = atoms.cell.lengths()
    gamma0 = atoms.cell.angles()[2]
    A_old = a0 * b0 * np.sin(np.radians(gamma0))

    atoms.calc = calc
    e0 = atoms.get_potential_energy()

    flt = FrechetCellFilter(atoms, mask=MASK)
    t0 = time.time()
    opt = BFGS(flt, logfile=str(RELAX_CR / f"{key}_cellrelax.log"))
    converged = opt.run(fmax=0.015, steps=1500)
    dt = time.time() - t0

    e1 = atoms.get_potential_energy()
    f1 = abs(atoms.get_forces()).max()
    a1, b1, c1 = atoms.cell.lengths()
    gamma1 = atoms.cell.angles()[2]
    A_new = a1 * b1 * np.sin(np.radians(gamma1))
    dA_pct = (A_new - A_old) / A_old * 100

    write(RELAX_CR / f"{key}_interface_cellrelax_PBE.cif", atoms)
    rec = {
        "tag": key, "n_atoms": n,
        "lattice_old": [float(a0), float(b0), float(c0), float(gamma0)],
        "lattice_new": [float(a1), float(b1), float(c1), float(gamma1)],
        "area_old_A2": float(A_old), "area_new_A2": float(A_new),
        "area_change_pct": float(dA_pct),
        "E_old_eV": float(e0), "E_new_eV": float(e1), "delta_E_eV": float(e1 - e0),
        "n_steps": int(opt.nsteps), "converged": bool(converged),
        "fmax_final": float(f1), "wall_seconds": float(dt),
    }
    records[key] = rec
    pair, term = key.rsplit("_term", 1)
    print(f"{pair:22s} {term:>4s} {n:>5d} "
          f"{a0:>7.3f} {b0:>7.3f} {gamma0:>6.2f} "
          f"{a1:>7.3f} {b1:>7.3f} {gamma1:>6.2f} "
          f"{A_old:>7.2f} {A_new:>7.2f} {dA_pct:>+6.2f} "
          f"{e0:>10.4f} {e1:>10.4f} {e1-e0:>+9.4f} "
          f"{opt.nsteps:>5d} {str(converged):>5s} {dt:>5.1f}")

(RES / "interface_cellrelax_PBE.json").write_text(json.dumps(records, indent=2))
print(f"\nSaved: {RES / 'interface_cellrelax_PBE.json'}")

# Now paper-style W_ad with cell-relaxed interface
print()
print("=" * 130)
print(" Paper-style W_ad with FULL cell relaxation")
print("=" * 130)

# Load alpha slab energies (these were relaxed at alpha natural lattice; reuse as is)
def load_e_alpha(key):
    return json.loads((RELAX / f"{key}_alpha_PBE.json").read_text())["energy_eV"]

# Load beta slab energies at beta NATURAL lattice from script 19
beta_nat = json.loads((RES / "wad_paper_style.json").read_text())["natural_beta_slabs"]

print(f"{'pair':22s} {'term':>4s} {'N':>5s} {'A_int_relaxed':>13s} "
      f"{'E_α':>10s} {'E_β,nat':>10s} {'E_int':>11s} "
      f"{'W_ad,paper':>11s} {'W_ad,prev':>10s}")
print("-" * 130)

table = []
for key, rec in records.items():
    pair, term = key.rsplit("_term", 1)
    e_alpha = load_e_alpha(key)
    e_beta_nat = beta_nat[pair]["E_beta_natural"]
    e_int = rec["E_new_eV"]
    A_new = rec["area_new_A2"]
    n = rec["n_atoms"]
    wad_paper = (e_alpha + e_beta_nat - e_int) * EV_TO_J / (A_new * ANG2_TO_M2)
    # Previous "paper-style" W_ad with α-fixed cell (from script 19)
    e_int_old = json.loads((RELAX / f"{key}_interface_PBE.json").read_text())["energy_eV"]
    A_old = rec["area_old_A2"]
    wad_prev = (e_alpha + e_beta_nat - e_int_old) * EV_TO_J / (A_old * ANG2_TO_M2)
    print(f"{pair:22s} {term:>4s} {n:>5d} {A_new:>13.2f} "
          f"{e_alpha:>10.4f} {e_beta_nat:>10.4f} {e_int:>11.4f} "
          f"{wad_paper:>11.3f} {wad_prev:>10.3f}")
    table.append({
        "pair": pair, "term": int(term), "method": "PBE_paper",
        "n_atoms_int": n,
        "A_int_relaxed_A2": A_new, "A_int_old_A2": A_old, "area_change_pct": rec["area_change_pct"],
        "E_alpha_eV": e_alpha, "E_beta_natural_eV": e_beta_nat,
        "E_int_relaxed_eV": e_int, "E_int_old_eV": e_int_old,
        "Wad_paper_full_Jm2": float(wad_paper),
        "Wad_paper_alpha_fixed_Jm2": float(wad_prev),
    })

print()
print("=" * 100)
print(" Pair-level max W_ad,paper-full (orientation preference)")
print("=" * 100)
by_pair = {}
for r in table:
    by_pair.setdefault(r["pair"], []).append(r)
ranked = sorted(((max(r["Wad_paper_full_Jm2"] for r in rs), p, rs) for p, rs in by_pair.items()), reverse=True)

print(f"{'rank':>4s} {'pair':22s} {'W_ad max':>10s} {'(prev)':>8s}")
for i, (w, p, rs) in enumerate(ranked):
    w_prev = max(r["Wad_paper_alpha_fixed_Jm2"] for r in rs)
    print(f"{i+1:>4d} {p:22s} {w:>10.3f} {w_prev:>+8.3f}")

(RES / "wad_paper_full.json").write_text(json.dumps({
    "interface_cellrelax_records": records,
    "table": table,
    "ranking_paper_full": [{"rank": i+1, "pair": p,
                           "Wad_paper_full_max": w} for i, (w, p, _) in enumerate(ranked)],
}, indent=2))
print(f"\nSaved: {RES / 'wad_paper_full.json'}")
