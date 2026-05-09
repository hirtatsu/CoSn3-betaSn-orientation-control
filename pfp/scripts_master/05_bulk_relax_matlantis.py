"""Bulk relaxation of alpha-CoSn3 and beta-Sn with PFP/PBE on Matlantis.

Run on Matlantis under ~/work_dir/cosn3_betasn_adhesion/.
Outputs:
  bulk_alpha_CoSn3_PBE.json: relaxed lattice + atomic positions + energy
  bulk_beta_Sn_PBE.json:     same
  bulk_alpha_CoSn3_PBE.cif:  relaxed structure (cif)
  bulk_beta_Sn_PBE.cif:      same

Settings (matched to Tatsumi-Ito 2026 manuscript):
  PFP v8, PBE calc mode
  ExpCellFilter (lattice + atomic positions free, isotropic stress allowed)
  fmax = 0.001 eV/A (manuscript Sec 2 says 0.001 for bulk in PFP)
  BFGS optimizer
  max steps = 500

Reference values to reproduce:
  beta-Sn PFP/PBE: a=5.9295, c=3.2008 A (mode-self-consistent)
  alpha-CoSn3:    new territory, no PFP reference
"""
import json, time
from pathlib import Path
from ase.io import read, write
from ase.optimize import BFGS
from ase.filters import ExpCellFilter
from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator
from pfp_api_client.pfp.estimator import Estimator, EstimatorCalcMode

WORK = Path("/home/jovyan/work_dir/cosn3_betasn_adhesion")
WORK.mkdir(parents=True, exist_ok=True)
DATA = WORK / "data"
RES = WORK / "results"
RES.mkdir(exist_ok=True)

calc = ASECalculator(Estimator(calc_mode=EstimatorCalcMode.PBE))

def relax_bulk(name, cif_in):
    print(f"\n[{name}] reading {cif_in}")
    atoms = read(cif_in)
    atoms.calc = calc
    n = len(atoms)
    print(f"  N atoms: {n}")
    print(f"  initial cell (a,b,c): {atoms.cell.lengths()}")
    print(f"  initial cell (al,be,ga): {atoms.cell.angles()}")
    e0 = atoms.get_potential_energy()
    print(f"  initial E: {e0:.6f} eV ({e0/n:.6f} eV/atom)")

    t0 = time.time()
    flt = ExpCellFilter(atoms)
    opt = BFGS(flt, logfile=str(RES / f"{name}_relax.log"))
    opt.run(fmax=0.001, steps=500)
    dt = time.time() - t0

    e1 = atoms.get_potential_energy()
    print(f"  relaxed cell (a,b,c): {atoms.cell.lengths()}")
    print(f"  relaxed cell (al,be,ga): {atoms.cell.angles()}")
    print(f"  final E:   {e1:.6f} eV ({e1/n:.6f} eV/atom)")
    print(f"  dE: {e1-e0:.6f} eV; relax wall {dt:.1f} s")

    cif_out = RES / f"bulk_{name}_PBE.cif"
    write(cif_out, atoms)
    json_out = RES / f"bulk_{name}_PBE.json"
    rec = {
        "name": name,
        "n_atoms": n,
        "cell_lengths": atoms.cell.lengths().tolist(),
        "cell_angles": atoms.cell.angles().tolist(),
        "volume": float(atoms.get_volume()),
        "volume_per_atom": float(atoms.get_volume() / n),
        "energy_eV": float(e1),
        "energy_per_atom": float(e1 / n),
        "wall_seconds": dt,
        "calc_mode": "PBE",
        "fmax_target": 0.001,
        "optimizer": "BFGS+ExpCellFilter",
    }
    with open(json_out, "w") as f:
        json.dump(rec, f, indent=2)
    print(f"  saved: {cif_out}, {json_out}")
    return rec

results = {}
results["alpha_CoSn3"] = relax_bulk("alpha_CoSn3", DATA / "alpha-CoSn3.cif")
results["beta_Sn"]     = relax_bulk("beta_Sn",     DATA / "beta-Sn.cif")

# Comparison summary vs experiment / PFP reference
print("\n" + "=" * 70)
print("Comparison to references")
print("=" * 70)
ref_beta = {"a": 5.9295, "c": 3.2008}  # manuscript PFP/PBE mode-self-consistent
ref_beta_exp = {"a": 5.831, "c": 3.182}
b = results["beta_Sn"]
print(f"beta-Sn PFP/PBE: a={b['cell_lengths'][0]:.4f} (ref {ref_beta['a']:.4f}, "
      f"exp {ref_beta_exp['a']:.4f}), c={b['cell_lengths'][2]:.4f} (ref {ref_beta['c']:.4f}, "
      f"exp {ref_beta_exp['c']:.4f})")
ref_acs_exp = {"a": 16.864, "b": 6.268, "c": 6.270}
a = results["alpha_CoSn3"]
print(f"alpha-CoSn3 PFP/PBE: a={a['cell_lengths'][0]:.4f} (exp {ref_acs_exp['a']:.4f}), "
      f"b={a['cell_lengths'][1]:.4f} (exp {ref_acs_exp['b']:.4f}), "
      f"c={a['cell_lengths'][2]:.4f} (exp {ref_acs_exp['c']:.4f})")

with open(RES / "bulk_summary_PBE.json", "w") as f:
    json.dump(results, f, indent=2)
print(f"\nSummary saved: {RES / 'bulk_summary_PBE.json'}")
