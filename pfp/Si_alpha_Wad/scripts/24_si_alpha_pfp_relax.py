"""PFP/PBE cell-relax for Si/α-CoSn3 W_ad (manuscript convention).

Following Tatsumi 2026 Eq.4 paper-style W_ad recipe (same as scripts 20+21):
  - Cell-relax interface with FrechetCellFilter mask=[T,T,F,F,F,T]
  - Extract Si and α subsets from interface, cell-relax each independently
    (FrechetCellFilter mask=[T,T,F,F,F,T]) → each finds own equilibrium
  - W_ad = (E_α + E_Si − E_int) × 1.602e-19 / (A_int × 1e-20)  [J/m²]

Stages:
  A. Si bulk relax (ExpCellFilter, all 6) — for μ_Si reference and lattice check
  B. Si(100) free slab cell-relax — cross-check; not used in W_ad directly
  C. Cell-relax 3 interfaces: Si/α(600) t0, t1, Si/α(312) t0
  D. Cell-relax 6 split slabs (3 stacks × 2 slabs)
  E. Compute W_ad → results/wad_Si_alpha_PBE.json

Run on Matlantis. Inputs in: ~/work_dir/cosn3_betasn_adhesion_pfp/Si_alpha_wad/stacks/
Results in:                    ~/work_dir/cosn3_betasn_adhesion_pfp/Si_alpha_wad/results/
"""
import json, time
from pathlib import Path
import numpy as np
from ase.io import read, write
from ase.optimize import BFGS
from ase.filters import FrechetCellFilter, ExpCellFilter
from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator
from pfp_api_client.pfp.estimator import Estimator, EstimatorCalcMode

EV_TO_J = 1.602176634e-19
ANG2_TO_M2 = 1e-20

WORK = Path("/home/jovyan/work_dir/cosn3_betasn_adhesion_pfp/Si_alpha_wad")
STK = WORK / "stacks"
REL = WORK / "relaxed"; REL.mkdir(exist_ok=True, parents=True)
RES = WORK / "results"; RES.mkdir(exist_ok=True)
LOG = WORK / "logs";    LOG.mkdir(exist_ok=True)

calc = ASECalculator(Estimator(calc_mode=EstimatorCalcMode.PBE))

MASK_SLAB = [True, True, False, False, False, True]   # in-plane only (xx,yy,xy)
FMAX = 0.015
MAX_STEPS = 2000

def cellrelax(cif_in, cif_out, log_path, mask=MASK_SLAB, full=False, fmax=FMAX, steps=MAX_STEPS):
    atoms = read(cif_in); atoms.calc = calc
    a0,b0,c0 = atoms.cell.lengths()
    g0 = atoms.cell.angles()[2]
    e0 = atoms.get_potential_energy()
    flt = ExpCellFilter(atoms) if full else FrechetCellFilter(atoms, mask=mask)
    t0 = time.time()
    opt = BFGS(flt, logfile=str(log_path))
    conv = opt.run(fmax=fmax, steps=steps)
    dt = time.time() - t0
    e1 = atoms.get_potential_energy()
    f1 = abs(atoms.get_forces()).max()
    a1,b1,c1 = atoms.cell.lengths()
    g1 = atoms.cell.angles()[2]
    write(cif_out, atoms)
    return {
        "n_atoms": len(atoms),
        "lattice_old": [float(a0), float(b0), float(c0), float(g0)],
        "lattice_new": [float(a1), float(b1), float(c1), float(g1)],
        "area_old_A2": float(a0*b0*np.sin(np.radians(g0))),
        "area_new_A2": float(a1*b1*np.sin(np.radians(g1))),
        "E_old_eV": float(e0), "E_new_eV": float(e1),
        "delta_E_eV": float(e1 - e0),
        "n_steps": int(opt.nsteps), "converged": bool(conv),
        "fmax_final": float(f1), "wall_seconds": float(dt),
    }

results = {"stage_A_bulk": {}, "stage_B_si_freeslab": {},
           "stage_C_interfaces": {}, "stage_D_split_slabs": {}, "stage_E_wad": {}}

# ----- Stage A -----
print("="*100); print(" Stage A: Si bulk full-cell relax (ExpCellFilter)"); print("="*100)
rec = cellrelax(STK / "Si_bulk.cif", REL / "Si_bulk.cif",
                LOG / "Si_bulk.log", full=True)
mu_Si = rec["E_new_eV"] / rec["n_atoms"]
rec["mu_Si_eV_per_atom"] = mu_Si
results["stage_A_bulk"] = rec
print(f"  Si bulk: a={rec['lattice_new'][0]:.4f} (was {rec['lattice_new'][0]:.4f}), "
      f"E={rec['E_new_eV']:.4f} eV, μ_Si={mu_Si:.4f} eV/atom, "
      f"steps={rec['n_steps']}, conv={rec['converged']}, t={rec['wall_seconds']:.1f}s")

# ----- Stage B -----
print(); print("="*100); print(" Stage B: Si(100) free slab cell-relax"); print("="*100)
rec = cellrelax(STK / "Si100_free_slab.cif", REL / "Si100_free_slab.cif",
                LOG / "Si100_free_slab.log")
results["stage_B_si_freeslab"] = rec
print(f"  Si free slab: N={rec['n_atoms']}, A={rec['area_new_A2']:.2f} Å², E={rec['E_new_eV']:.4f}, "
      f"steps={rec['n_steps']}, conv={rec['converged']}, t={rec['wall_seconds']:.1f}s")

# ----- Stage C -----
print(); print("="*100); print(" Stage C: Interface cell-relax"); print("="*100)
print(f"{'tag':30s} {'N':>5s} {'a_old':>6s} {'b_old':>6s} {'γ_old':>6s} "
      f"{'a_new':>6s} {'b_new':>6s} {'γ_new':>6s} {'A_new':>7s} {'ΔA%':>6s} "
      f"{'E_eV':>11s} {'steps':>5s} {'conv':>5s} {'t[s]':>5s}")
print("-"*120)
INTERFACE_TAGS = [
    "Si100_alpha600_term0", "Si100_alpha600_term1",
    "Si100_alpha312_term0",
]
for tag in INTERFACE_TAGS:
    rec = cellrelax(STK / f"{tag}.cif", REL / f"{tag}_int.cif",
                    LOG / f"{tag}_int.log")
    results["stage_C_interfaces"][tag] = rec
    da = (rec["area_new_A2"]-rec["area_old_A2"])/rec["area_old_A2"]*100
    print(f"{tag:30s} {rec['n_atoms']:>5d} "
          f"{rec['lattice_old'][0]:>6.2f} {rec['lattice_old'][1]:>6.2f} {rec['lattice_old'][3]:>6.1f} "
          f"{rec['lattice_new'][0]:>6.2f} {rec['lattice_new'][1]:>6.2f} {rec['lattice_new'][3]:>6.1f} "
          f"{rec['area_new_A2']:>7.2f} {da:>+6.2f} "
          f"{rec['E_new_eV']:>11.4f} {rec['n_steps']:>5d} {str(rec['converged']):>5s} {rec['wall_seconds']:>5.1f}")

# ----- Stage D -----
print(); print("="*100); print(" Stage D: Split slabs cell-relax (Si side + α side per stack)"); print("="*100)
print(f"{'tag':35s} {'N':>4s} {'A_old':>7s} {'A_new':>7s} {'ΔA%':>6s} "
      f"{'E_eV':>11s} {'steps':>5s} {'conv':>5s} {'t[s]':>5s}")
print("-"*100)
for tag in INTERFACE_TAGS:
    for side in ["Si", "alpha"]:
        sub = f"{tag}_{side}"
        rec = cellrelax(STK / f"{sub}.cif", REL / f"{sub}.cif",
                        LOG / f"{sub}.log")
        results["stage_D_split_slabs"][sub] = rec
        da = (rec["area_new_A2"]-rec["area_old_A2"])/rec["area_old_A2"]*100
        print(f"{sub:35s} {rec['n_atoms']:>4d} "
              f"{rec['area_old_A2']:>7.2f} {rec['area_new_A2']:>7.2f} {da:>+6.2f} "
              f"{rec['E_new_eV']:>11.4f} {rec['n_steps']:>5d} {str(rec['converged']):>5s} {rec['wall_seconds']:>5.1f}")

# ----- Stage E -----
print(); print("="*100); print(" Stage E: W_ad (paper-style, Tatsumi 2026 Eq.4)"); print("="*100)
print(f"  W_ad = (E_α + E_Si − E_int) × eV→J / (A_int × Å²→m²)")
print(f"{'tag':30s} {'A_int':>7s} {'E_int':>11s} {'E_Si':>11s} {'E_α':>11s} {'W_ad [J/m²]':>12s}")
print("-"*100)
wad_table = []
for tag in INTERFACE_TAGS:
    e_int = results["stage_C_interfaces"][tag]["E_new_eV"]
    A_int = results["stage_C_interfaces"][tag]["area_new_A2"]
    e_Si  = results["stage_D_split_slabs"][f"{tag}_Si"]["E_new_eV"]
    e_a   = results["stage_D_split_slabs"][f"{tag}_alpha"]["E_new_eV"]
    wad = (e_a + e_Si - e_int) * EV_TO_J / (A_int * ANG2_TO_M2)
    wad_table.append({"tag": tag, "A_int_A2": A_int, "E_int_eV": e_int,
                      "E_Si_eV": e_Si, "E_alpha_eV": e_a, "Wad_Jm2": float(wad)})
    print(f"{tag:30s} {A_int:>7.2f} {e_int:>11.4f} {e_Si:>11.4f} {e_a:>11.4f} {wad:>12.4f}")

# Pair-level max W_ad (max over terminations)
by_pair = {}
for r in wad_table:
    pair = r["tag"].rsplit("_term", 1)[0]
    by_pair.setdefault(pair, []).append(r)
ranking = []
for pair, rs in by_pair.items():
    wmax = max(r["Wad_Jm2"] for r in rs)
    wmin = min(r["Wad_Jm2"] for r in rs)
    ranking.append({"pair": pair, "Wad_max_Jm2": wmax, "Wad_min_Jm2": wmin,
                    "n_terms": len(rs)})
ranking.sort(key=lambda r: -r["Wad_max_Jm2"])

print()
print(" Pair-level summary (max/min over terminations):")
print(f"{'rank':>4s} {'pair':30s} {'W_ad max [J/m²]':>16s} {'W_ad min':>10s} {'#term':>5s}")
for i, r in enumerate(ranking):
    print(f"{i+1:>4d} {r['pair']:30s} {r['Wad_max_Jm2']:>16.4f} {r['Wad_min_Jm2']:>10.4f} {r['n_terms']:>5d}")

results["stage_E_wad"]["table"] = wad_table
results["stage_E_wad"]["ranking"] = ranking

# Save
out_json = RES / "wad_Si_alpha_PBE.json"
out_json.write_text(json.dumps(results, indent=2))
print(f"\nSaved: {out_json}")
