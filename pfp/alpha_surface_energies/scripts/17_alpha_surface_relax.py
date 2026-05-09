"""Relax pymatgen-generated alpha-CoSn3 slabs (600, 010, 301; all terminations)
plus reuse ASE results for (312) and (321) which gave sensible values.
"""
import json, time
from pathlib import Path
import numpy as np
from ase.io import read, write
from ase.optimize import BFGS
from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator
from pfp_api_client.pfp.estimator import Estimator, EstimatorCalcMode

EV_TO_J = 1.602176634e-19; ANG2_TO_M2 = 1e-20

WORK = Path("/home/jovyan/work_dir/cosn3_betasn_adhesion")
SLABS = WORK / "alpha_slabs"
SURF = WORK / "alpha_surfaces_pmg"; SURF.mkdir(exist_ok=True)
RES = WORK / "results"

calc = ASECalculator(Estimator(calc_mode=EstimatorCalcMode.PBE))

bulk = read(RES / "bulk_alpha_CoSn3_PBE.cif")
bulk.calc = calc
e_bulk_per_atom = bulk.get_potential_energy() / len(bulk)
print(f"Bulk PBE: E/atom = {e_bulk_per_atom:.6f} eV")

INDEX = json.loads((WORK / "alpha_slabs_index.json").read_text())
TARGETS = ["(6, 0, 0)", "(0, 1, 0)", "(3, 0, 1)"]  # only small cells

results = {}
print()
print("=" * 110)
print(f"{'face':>8s} {'term':>4s} {'N':>5s} {'A [A^2]':>9s} {'E_slab [eV]':>14s} "
      f"{'γ [J/m²]':>10s} {'fmax':>7s} {'conv':>5s} {'t[s]':>6s}")
print("=" * 110)

for face_key in TARGETS:
    info = INDEX[face_key]
    print()
    face_recs = []
    for term in info["terminations"]:
        cif_path = SLABS / term["cif"]
        atoms = read(cif_path); atoms.calc = calc
        n = len(atoms)
        a, b, _ = atoms.cell.lengths()
        gamma_deg = atoms.cell.angles()[2]
        area = a * b * np.sin(np.radians(gamma_deg))
        e0 = atoms.get_potential_energy()
        t0 = time.time()
        opt = BFGS(atoms, logfile=str(SURF / f"{term['tag']}.log"))
        converged = opt.run(fmax=0.015, steps=500)
        dt = time.time() - t0
        e1 = atoms.get_potential_energy()
        f1 = abs(atoms.get_forces()).max()
        write(SURF / f"{term['tag']}.cif", atoms)
        gamma_J = (e1 - n * e_bulk_per_atom) * EV_TO_J / (2 * area * ANG2_TO_M2)
        face_recs.append({
            "termination_idx": term["termination_idx"],
            "n_atoms": n, "area_A2": float(area),
            "energy_eV": float(e1), "gamma_Jm2": float(gamma_J),
            "fmax": float(f1), "converged": bool(converged), "wall_s": float(dt),
        })
        print(f"{face_key:>8s} {term['termination_idx']:>4d} {n:>5d} {area:>9.2f} "
              f"{e1:>14.4f} {gamma_J:>10.3f} {f1:>7.4f} {str(converged):>5s} {dt:>6.1f}")

    gmin = min(r["gamma_Jm2"] for r in face_recs)
    gmax = max(r["gamma_Jm2"] for r in face_recs)
    results[face_key] = {
        "hkl_nominal": info["hkl_nominal"],
        "terminations": face_recs,
        "gamma_min_Jm2": float(gmin), "gamma_max_Jm2": float(gmax),
        "n_terminations": len(face_recs),
        "gamma_ref_Table1": info["gamma_ref_Table1"],
    }
    print(f"  -> {face_key} γ_min={gmin:.3f}, γ_max={gmax:.3f}")

print()
print("=" * 80)
print(" Final ranking (lowest γ per face from full termination enumeration)")
print("=" * 80)
sorted_pbe = sorted(results.items(), key=lambda kv: kv[1]["gamma_min_Jm2"])
print(f"{'face':>10s} {'γ_min':>8s} {'γ_max':>8s} {'γ_ref':>7s} {'Δ_min':>7s} {'#terms':>6s}")
for k, v in sorted_pbe:
    ref = v["gamma_ref_Table1"]
    delta = (v["gamma_min_Jm2"] - ref) if ref else None
    delta_str = f"{delta:+.2f}" if delta is not None else "—"
    ref_str = f"{ref:.2f}" if ref else "—"
    print(f"{k:>10s} {v['gamma_min_Jm2']:>8.3f} {v['gamma_max_Jm2']:>8.3f} "
          f"{ref_str:>7s} {delta_str:>7s} {v['n_terminations']:>6d}")

(RES / "alpha_surfaces_pmg_PBE.json").write_text(json.dumps(results, indent=2))
print(f"\nSaved: {RES}/alpha_surfaces_pmg_PBE.json")
