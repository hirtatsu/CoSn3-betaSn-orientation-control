"""Si bulk re-relax at fmax=0.001 (manuscript spec).

Supplement S2.1.1 declares "ExpCellFilter ... fmax 0.001 eV/Å" for bulk.
Earlier script 24 used fmax=0.015 → a=5.4652 Å.
This script tightens fmax to 0.001 to align with manuscript spec.

Output:
  Si_alpha_wad/relaxed/Si_bulk_tight.cif  (relaxed structure)
  Si_alpha_wad/results/Si_bulk_tight.json (energy + lattice)
  Si_alpha_wad/logs/Si_bulk_tight.log     (BFGS log)
"""
import json, time
from pathlib import Path
from ase.io import read, write
from ase.optimize import BFGS
from ase.filters import ExpCellFilter
from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator
from pfp_api_client.pfp.estimator import Estimator, EstimatorCalcMode

WORK = Path("/home/jovyan/work_dir/cosn3_betasn_adhesion_pfp/Si_alpha_wad")
STK = WORK / "stacks"
REL = WORK / "relaxed"
RES = WORK / "results"
LOG = WORK / "logs"

calc = ASECalculator(Estimator(calc_mode=EstimatorCalcMode.PBE))

atoms = read(STK / "Si_bulk.cif")
atoms.calc = calc
print(f"Si bulk input: a={atoms.cell.lengths()[0]:.5f} Å, N={len(atoms)}")
e0 = atoms.get_potential_energy()
print(f"  E_initial = {e0:.6f} eV")

flt = ExpCellFilter(atoms)
t0 = time.time()
opt = BFGS(flt, logfile=str(LOG / "Si_bulk_tight.log"))
conv = opt.run(fmax=0.001, steps=500)
dt = time.time() - t0

a, b, c = atoms.cell.lengths()
e1 = atoms.get_potential_energy()
f1 = abs(atoms.get_forces()).max()
mu_Si = e1 / len(atoms)

print(f"  E_final   = {e1:.6f} eV")
print(f"  μ_Si      = {mu_Si:.6f} eV/atom")
print(f"  a, b, c   = {a:.5f}, {b:.5f}, {c:.5f} Å")
print(f"  fmax_final = {f1:.6f} eV/Å (target 0.001)")
print(f"  steps={opt.nsteps}, conv={conv}, t={dt:.1f}s")

write(REL / "Si_bulk_tight.cif", atoms)
out = {
    "fmax_target": 0.001,
    "fmax_final": float(f1),
    "converged": bool(conv),
    "n_steps": int(opt.nsteps),
    "n_atoms": len(atoms),
    "lattice": [float(a), float(b), float(c)],
    "E_eV": float(e1),
    "mu_Si_eV_per_atom": float(mu_Si),
    "wall_seconds": float(dt),
    "vs_loose_run": {
        "loose_a_A": 5.4652,
        "loose_mu_Si_eV": -4.5526,
        "delta_a_A": float(a) - 5.4652,
        "delta_mu_meV": (mu_Si - (-4.5526)) * 1000,
    },
    "vs_manuscript": {
        "manuscript_a_A": 5.469,
        "delta_a_A": float(a) - 5.469,
    },
}
(RES / "Si_bulk_tight.json").write_text(json.dumps(out, indent=2))
print(f"\nSaved: {RES / 'Si_bulk_tight.json'}")
