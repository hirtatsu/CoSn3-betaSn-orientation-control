"""Cell-relax all standalone slabs (alpha and beta) for paper-style consistency.

(a) Alpha slabs: 9 (3 faces × 3 terminations) — restart from atom-relaxed CIFs
    in alpha_surfaces_pmg/{tag}.cif
(b) Beta natural slabs: 6 (one per pair) — restart from atom-relaxed CIFs
    in relax_beta_natural/{pair}_beta_natural.cif

Mask: [T,T,F,F,F,T] — relax in-plane (xx, yy, xy=gamma); fix vacuum (zz, yz, xz)

Output:
  alpha_surfaces_pmg_cellrelax/{tag}.{cif,json}
  relax_beta_natural_cellrelax/{pair}_beta_natural.{cif,json}
  results/alpha_surfaces_cellrelax_PBE.json
  results/wad_paper_full_consistent.json (final W_ad,paper after slab cell relax)
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
ALPHA_OLD = WORK / "alpha_surfaces_pmg"
ALPHA_NEW = WORK / "alpha_surfaces_pmg_cellrelax"; ALPHA_NEW.mkdir(exist_ok=True)
BETA_OLD = WORK / "relax_beta_natural"
BETA_NEW = WORK / "relax_beta_natural_cellrelax"; BETA_NEW.mkdir(exist_ok=True)
RES = WORK / "results"
RELAX_CR = WORK / "relax_cell"   # cell-relaxed interfaces
RELAX = WORK / "relax"            # atom-only relaxed (for E_alpha-slab)

calc = ASECalculator(Estimator(calc_mode=EstimatorCalcMode.PBE))
MASK = [True, True, False, False, False, True]

def cellrelax_slab(cif_in, cif_out, log_path, fmax=0.015, max_steps=1500):
    atoms = read(cif_in); atoms.calc = calc
    a0, b0, c0 = atoms.cell.lengths()
    gamma0 = atoms.cell.angles()[2]
    e0 = atoms.get_potential_energy()
    flt = FrechetCellFilter(atoms, mask=MASK)
    t0 = time.time()
    opt = BFGS(flt, logfile=str(log_path))
    converged = opt.run(fmax=fmax, steps=max_steps)
    dt = time.time() - t0
    e1 = atoms.get_potential_energy()
    f1 = abs(atoms.get_forces()).max()
    a1, b1, c1 = atoms.cell.lengths()
    gamma1 = atoms.cell.angles()[2]
    write(cif_out, atoms)
    return {
        "n_atoms": len(atoms),
        "lattice_old": [float(a0), float(b0), float(c0), float(gamma0)],
        "lattice_new": [float(a1), float(b1), float(c1), float(gamma1)],
        "area_old": float(a0*b0*np.sin(np.radians(gamma0))),
        "area_new": float(a1*b1*np.sin(np.radians(gamma1))),
        "E_old": float(e0), "E_new": float(e1), "delta_E": float(e1-e0),
        "n_steps": int(opt.nsteps), "converged": bool(converged),
        "fmax": float(f1), "wall_s": float(dt),
    }

# (a) Alpha slabs
print("=" * 110)
print(" Alpha slab cell relaxation (3 faces × 3 terms)")
print("=" * 110)
print(f"{'tag':25s} {'N':>4s} {'A_old':>7s} {'A_new':>7s} {'ΔA%':>6s} "
      f"{'E_old':>11s} {'E_new':>11s} {'ΔE':>9s} {'steps':>5s} {'conv':>5s} {'γ_PBE':>7s}")
print("-" * 110)

bulk_alpha = read(RES / "bulk_alpha_CoSn3_PBE.cif")
bulk_alpha.calc = calc
mu_alpha = bulk_alpha.get_potential_energy() / len(bulk_alpha)

alpha_results = {}
for face_dir in ["600", "010", "301"]:
    for ti in [0, 1, 2]:
        tag = f"alpha_{face_dir}_t{ti}"
        cif_in = ALPHA_OLD / f"{tag}.cif"
        if not cif_in.exists():
            print(f"  [missing] {cif_in.name}"); continue
        cif_out = ALPHA_NEW / f"{tag}.cif"
        log = ALPHA_NEW / f"{tag}.log"
        r = cellrelax_slab(cif_in, cif_out, log)
        # gamma_min = (E_slab - N*mu_alpha) / (2*A) at NEW lattice
        n = r["n_atoms"]; A = r["area_new"]; E = r["E_new"]
        gamma_J = (E - n * mu_alpha) * EV_TO_J / (2 * A * ANG2_TO_M2)
        r["gamma_Jm2"] = float(gamma_J)
        alpha_results[tag] = r
        dA = (r["area_new"] - r["area_old"]) / r["area_old"] * 100
        print(f"{tag:25s} {n:>4d} {r['area_old']:>7.2f} {r['area_new']:>7.2f} "
              f"{dA:>+6.2f} {r['E_old']:>11.4f} {r['E_new']:>11.4f} {r['delta_E']:>+9.4f} "
              f"{r['n_steps']:>5d} {str(r['converged']):>5s} {gamma_J:>7.3f}")

(RES / "alpha_surfaces_cellrelax_PBE.json").write_text(json.dumps(alpha_results, indent=2))

# Summary: gamma_min per face
print()
print(" Alpha γ_min after cell relaxation (vs prior atom-only):")
import re
prior = json.loads((RES / "alpha_surfaces_pmg_PBE.json").read_text())
for face_key in ["(6, 0, 0)", "(0, 1, 0)", "(3, 0, 1)"]:
    face_dir = ''.join(re.findall(r'\d', face_key))
    gammas = [v["gamma_Jm2"] for k, v in alpha_results.items() if k.startswith(f"alpha_{face_dir}_")]
    if not gammas: continue
    g_min_new = min(gammas); g_max_new = max(gammas)
    g_min_old = prior[face_key]["gamma_min_Jm2"]
    g_max_old = prior[face_key]["gamma_max_Jm2"]
    ref = prior[face_key]["gamma_ref_Table1"]
    ref_str = f"{ref:.2f}" if ref else "—"
    print(f"   {face_key}: γ_min new={g_min_new:.3f} (old {g_min_old:.3f}, ref {ref_str}); "
          f"γ_max new={g_max_new:.3f} (old {g_max_old:.3f})")

# (b) Beta natural slabs
print()
print("=" * 110)
print(" Beta natural slab cell relaxation (1 per pair)")
print("=" * 110)
print(f"{'pair':22s} {'N':>4s} {'A_old':>7s} {'A_new':>7s} {'ΔA%':>6s} "
      f"{'E_old':>11s} {'E_new':>11s} {'ΔE':>9s} {'steps':>5s} {'conv':>5s}")
print("-" * 110)

bulk_beta = read(RES / "bulk_beta_Sn_PBE.cif")
bulk_beta.calc = calc
mu_beta = bulk_beta.get_potential_energy() / len(bulk_beta)

idx = json.loads((WORK / "slab_split_index.json").read_text())
unique_pairs = {}
for key in idx:
    pair = key.rsplit("_term", 1)[0]
    if pair not in unique_pairs:
        unique_pairs[pair] = key

beta_natural_new = {}
for pair, key in unique_pairs.items():
    cif_in = BETA_OLD / f"{pair}_beta_natural.cif"
    if not cif_in.exists():
        print(f"  [missing] {cif_in.name}"); continue
    cif_out = BETA_NEW / f"{pair}_beta_natural.cif"
    log = BETA_NEW / f"{pair}_beta_natural.log"
    r = cellrelax_slab(cif_in, cif_out, log)
    beta_natural_new[pair] = r
    dA = (r["area_new"] - r["area_old"]) / r["area_old"] * 100
    print(f"{pair:22s} {r['n_atoms']:>4d} {r['area_old']:>7.2f} {r['area_new']:>7.2f} "
          f"{dA:>+6.2f} {r['E_old']:>11.4f} {r['E_new']:>11.4f} {r['delta_E']:>+9.4f} "
          f"{r['n_steps']:>5d} {str(r['converged']):>5s}")

(RES / "beta_natural_cellrelax_PBE.json").write_text(json.dumps(beta_natural_new, indent=2))

# (c) Recompute paper-style W_ad with new slab energies
print()
print("=" * 130)
print(" Updated W_ad,paper with cell-relaxed standalone slabs")
print("=" * 130)
print(f"{'pair':22s} {'term':>4s} {'A_int':>8s} {'E_α':>10s} {'E_β,nat':>10s} {'E_int':>11s} "
      f"{'W_ad,prev':>10s} {'W_ad,new':>10s} {'Δ':>7s}")
print("-" * 130)

ifc_records = json.loads((RES / "interface_cellrelax_PBE.json").read_text())
prev_table = json.loads((RES / "wad_paper_full.json").read_text())["table"]
prev_lookup = {(r["pair"], r["term"]): r["Wad_paper_full_Jm2"] for r in prev_table}

# Need alpha slab energies at cell-relaxed lattice. The alpha-only slabs from interface
# split (in relax/{key}_alpha_PBE.cif) had cell pinned at alpha lattice. Cell-relax those.
# But those are the SAME α(600) etc. terminations as in alpha_surfaces_pmg.
# Use the interface-side α slab (same termination as appears at the interface).
# For consistency, we cell-relax the interface-side α slabs too.
print("\n  (Cell-relaxing interface-side α slabs for proper consistency...)")
ALPHA_IFC_NEW = WORK / "relax_alpha_cellrelax"; ALPHA_IFC_NEW.mkdir(exist_ok=True)
alpha_ifc_results = {}
for key in idx:
    cif_in = RELAX / f"{key}_alpha_PBE.cif"
    cif_out = ALPHA_IFC_NEW / f"{key}_alpha_cellrelax.cif"
    log = ALPHA_IFC_NEW / f"{key}_alpha_cellrelax.log"
    r = cellrelax_slab(cif_in, cif_out, log)
    alpha_ifc_results[key] = r

(RES / "alpha_ifc_cellrelax_PBE.json").write_text(json.dumps(alpha_ifc_results, indent=2))

table = []
for key, ifc_rec in ifc_records.items():
    pair, term = key.rsplit("_term", 1)
    e_alpha = alpha_ifc_results[key]["E_new"]
    e_beta_nat = beta_natural_new[pair]["E_new"]
    e_int = ifc_rec["E_new_eV"]
    A_int = ifc_rec["area_new_A2"]
    wad_new = (e_alpha + e_beta_nat - e_int) * EV_TO_J / (A_int * ANG2_TO_M2)
    wad_prev = prev_lookup.get((pair, int(term)), None)
    diff = wad_new - wad_prev if wad_prev is not None else None
    print(f"{pair:22s} {term:>4s} {A_int:>8.2f} {e_alpha:>10.4f} {e_beta_nat:>10.4f} {e_int:>11.4f} "
          f"{wad_prev if wad_prev is not None else 0:>10.3f} {wad_new:>10.3f} "
          f"{diff if diff is not None else 0:>+7.3f}")
    table.append({
        "pair": pair, "term": int(term), "method": "PBE_paper_full_consistent",
        "A_int": float(A_int), "E_alpha": float(e_alpha), "E_beta_nat": float(e_beta_nat),
        "E_int": float(e_int), "Wad_paper_full_consistent": float(wad_new),
        "Wad_paper_full_prev": wad_prev, "delta": diff,
    })

print()
print("=" * 80)
print(" Pair-level max W_ad,paper (fully cell-relaxed) ranking")
print("=" * 80)
by_pair = {}
for r in table:
    by_pair.setdefault(r["pair"], []).append(r)
ranked = sorted(((max(r["Wad_paper_full_consistent"] for r in rs), p, rs) for p, rs in by_pair.items()), reverse=True)

print(f"{'rank':>4s} {'pair':22s} {'W_ad max':>10s} {'(prev step)':>12s} {'c-axis':>8s}")
c_axis = {"alpha600_beta100": "∥", "alpha600_beta001": "⊥", "alpha600_beta110": "∥",
          "alpha312_beta100": "∥", "alpha312_beta001": "⊥", "alpha010_beta100": "∥"}
for i, (w, p, rs) in enumerate(ranked):
    w_prev = max(r["Wad_paper_full_prev"] for r in rs) if rs[0]["Wad_paper_full_prev"] else None
    print(f"{i+1:>4d} {p:22s} {w:>10.3f} {(w_prev if w_prev else 0):>12.3f} {c_axis.get(p, '?'):>8s}")

(RES / "wad_paper_full_consistent.json").write_text(json.dumps({
    "alpha_ifc_cellrelax": alpha_ifc_results,
    "beta_natural_cellrelax": beta_natural_new,
    "alpha_surfaces_cellrelax": alpha_results,
    "table": table,
    "ranking": [{"rank": i+1, "pair": p,
                 "Wad_paper_full_consistent_max": w} for i, (w, p, _) in enumerate(ranked)],
}, indent=2))
print(f"\nSaved: {RES / 'wad_paper_full_consistent.json'}")
