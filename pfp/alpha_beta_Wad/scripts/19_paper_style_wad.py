"""Compute paper-style W_ad following Tatsumi 2026 Acta Mater 304, 121813 Eq. (4):

W_ad,paper = (E_slab,A + E_slab,B − E_interface) / A

where E_slab,A and E_slab,B are at their NATURAL (relaxed) lattices (no strain),
and E_interface is at the constrained (α) lattice.

Differs from our previous "W_ad" by the strain energy in the β slab:
W_ad,paper = W_ad,ours − (strain energy density of β slab) / A

Procedure:
1. For each unique β slab in the 6 pairs, take the strained-β slab CIF.
2. Rescale the in-plane cell from α-strained to β-natural using:
   new_a = old_a / (1 + strain_a/100)
   new_b = old_b / (1 + strain_b/100)
   with maintained fractional coordinates.
3. Relax atomic positions (cell fixed at natural lattice).
4. Get E_β-slab,natural.
5. Compute paper-style W_ad and compare to ours.
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
SLABS = WORK / "slabs"
RELAX = WORK / "relax"
RELAX_NAT = WORK / "relax_beta_natural"; RELAX_NAT.mkdir(exist_ok=True)
RES = WORK / "results"

calc = ASECalculator(Estimator(calc_mode=EstimatorCalcMode.PBE))

idx = json.loads((WORK / "slab_split_index.json").read_text())

# For each pair, take ONE term's beta slab (term0; the strain is the same regardless of term)
unique_pairs = {}
for key, info in idx.items():
    pair = key.rsplit("_term", 1)[0]
    if pair not in unique_pairs:
        unique_pairs[pair] = (key, info)

print("=" * 110)
print(f"{'pair':22s} {'N':>4s} {'A_strained':>11s} {'A_natural':>11s} "
      f"{'εa%':>6s} {'εb%':>6s} {'E_β,strained':>14s} {'E_β,natural':>13s} {'ΔE_strain [eV]':>14s}")
print("=" * 110)

results_natural = {}
for pair, (key, info) in unique_pairs.items():
    cif_in = SLABS / info["beta_slab_cif"]
    atoms_strained = read(cif_in)
    n = len(atoms_strained)
    a_str, b_str, c_str = atoms_strained.cell.lengths()
    gamma_str = atoms_strained.cell.angles()[2]

    # Strain components from index
    eps_a = info["strain_a_pct"] / 100.0
    eps_b = info["strain_b_pct"] / 100.0
    # β natural lattice: a_β = a_α / (1 + ε_a)... no: ε = (β-α)/α => β = α(1+ε)
    # So a_α = a_β/(1+ε), or a_β = a_α(1+ε).
    # Wait: in pymatgen ZSL, u_a = (f_a - s_a)/s_a where s = substrate=α, f = film=β.
    # So f_a = s_a (1 + u_a). f_a is β's lattice IN INTERFACE CELL = s_a (since β strained to α).
    # Actually IN MATCHED CELL, both α and β have same length s_a. The strain reported is
    # how much β had to compress/expand. β natural would be β-without-strain.
    # If ε_a = (f_a-s_a)/s_a is reported negative (e.g., -5.14%), it means β had to
    # contract by 5.14% to match α (β natural was bigger by 5.14%).
    # So β natural a = s_a * (1 + |ε|) when ε is "how much β contracted/expanded to fit α".
    #
    # Sign convention check from script 06:
    # u_a = (f_a - s_a)/s_a * 100, where s_a comes from substrate (α), f_a from film (β).
    # For α(600)/β(100): reported -5.14% means f_a < s_a, i.e., β natural < α lattice.
    # So β had to STRETCH to match α.
    # Hence to recover β natural: β_a = s_a * (1 + ε_a/100) [where ε is in %, signed].
    # i.e., β_a = a_strained × (1 + ε_a/100).
    # Sanity: α(600)/β(100): a_strained = 6.25, ε = -5.14% → β_a = 6.25 × 0.9486 = 5.929 ✓
    factor_a = 1.0 + eps_a
    factor_b = 1.0 + eps_b

    # Build new cell. Keep c (vacuum) the same. Scale a,b vectors.
    cell = atoms_strained.cell.copy()
    cell[0] = cell[0] * factor_a
    cell[1] = cell[1] * factor_b
    # c unchanged

    atoms_natural = atoms_strained.copy()
    # Use scaled positions to preserve relative atomic arrangement
    scaled = atoms_strained.get_scaled_positions(wrap=False)
    atoms_natural.set_cell(cell, scale_atoms=False)
    atoms_natural.set_scaled_positions(scaled)

    a_nat, b_nat, c_nat = atoms_natural.cell.lengths()
    A_natural = a_nat * b_nat * np.sin(np.radians(atoms_natural.cell.angles()[2]))
    A_strained = a_str * b_str * np.sin(np.radians(gamma_str))

    # Relax atoms only (cell fixed at natural lattice)
    atoms_natural.calc = calc
    e0 = atoms_natural.get_potential_energy()
    t0 = time.time()
    opt = BFGS(atoms_natural, logfile=str(RELAX_NAT / f"{pair}_beta_natural.log"))
    converged = opt.run(fmax=0.015, steps=500)
    dt = time.time() - t0
    e1 = atoms_natural.get_potential_energy()
    f1 = abs(atoms_natural.get_forces()).max()
    write(RELAX_NAT / f"{pair}_beta_natural.cif", atoms_natural)

    # E_β,strained from existing relaxed JSON
    e_beta_strained = json.loads((RELAX / f"{key}_beta_PBE.json").read_text())["energy_eV"]
    delta_E_strain = e_beta_strained - e1  # positive = strain raises energy
    delta_per_area = delta_E_strain * EV_TO_J / (A_strained * ANG2_TO_M2)  # J/m^2

    results_natural[pair] = {
        "n_beta": n, "A_strained": float(A_strained), "A_natural": float(A_natural),
        "strain_a_pct": float(eps_a*100), "strain_b_pct": float(eps_b*100),
        "E_beta_strained": float(e_beta_strained),
        "E_beta_natural": float(e1),
        "delta_E_strain_eV": float(delta_E_strain),
        "strain_energy_per_area_Jm2": float(delta_per_area),
        "natural_relax_steps": opt.nsteps, "converged": bool(converged),
        "wall_seconds": float(dt), "fmax_final": float(f1),
    }
    print(f"{pair:22s} {n:>4d} {A_strained:>11.2f} {A_natural:>11.2f} "
          f"{eps_a*100:>+6.2f} {eps_b*100:>+6.2f} {e_beta_strained:>14.4f} {e1:>13.4f} {delta_E_strain:>+14.4f}")

print("=" * 110)

# Now compute paper-style W_ad for all 10 stacks (re-using stacked α data)
INDIV = WORK / "relax"

def load_pair(key):
    info = idx[key]
    ea = json.loads((INDIV / f"{key}_alpha_PBE.json").read_text())["energy_eV"]
    eb_strained = json.loads((INDIV / f"{key}_beta_PBE.json").read_text())["energy_eV"]
    ei = json.loads((INDIV / f"{key}_interface_PBE.json").read_text())["energy_eV"]
    return ea, eb_strained, ei, info

print()
print("=" * 110)
print(" Paper-style W_ad vs ours-style W_ad (PFP/PBE)")
print("=" * 110)
print(f"{'pair':22s} {'term':>4s} {'εa%':>6s} {'εb%':>6s} "
      f"{'W_ad ours':>10s} {'W_ad paper':>11s} {'Δ(strain)':>11s} {'A [Å²]':>8s}")
print("-" * 110)

table = []
for key, info in idx.items():
    pair = key.rsplit("_term", 1)[0]
    term = int(key.rsplit("_term", 1)[1])
    ea, eb_s, ei, _ = load_pair(key)
    eb_n = results_natural[pair]["E_beta_natural"]
    a, b = info["lattice_a"], info["lattice_b"]; gamma = info["gamma"]
    A_int = a * b * np.sin(np.radians(gamma))
    wad_ours = (ea + eb_s - ei) * EV_TO_J / (A_int * ANG2_TO_M2)
    wad_paper = (ea + eb_n - ei) * EV_TO_J / (A_int * ANG2_TO_M2)
    strain_term = wad_paper - wad_ours
    print(f"{pair:22s} {term:>4d} {info['strain_a_pct']:>+6.2f} {info['strain_b_pct']:>+6.2f} "
          f"{wad_ours:>10.3f} {wad_paper:>11.3f} {strain_term:>+11.3f} {A_int:>8.2f}")
    table.append({
        "pair": pair, "term": term, "method": "PBE",
        "Wad_ours": float(wad_ours), "Wad_paper": float(wad_paper),
        "strain_correction": float(strain_term),
        "strain_a_pct": info["strain_a_pct"], "strain_b_pct": info["strain_b_pct"],
        "area_A2": float(A_int),
        "E_alpha_eV": ea, "E_beta_strained_eV": eb_s, "E_beta_natural_eV": eb_n, "E_int_eV": ei,
    })

print("=" * 110)

# Summary by pair (max-over-terms for paper-style)
print()
print("=" * 80)
print(" Pair-level max W_ad (paper-style) ranking")
print("=" * 80)
by_pair = {}
for r in table:
    by_pair.setdefault(r["pair"], []).append(r)

ranked_paper = []
for pair, rs in by_pair.items():
    w_max = max(r["Wad_paper"] for r in rs)
    ranked_paper.append((w_max, pair, rs))
ranked_paper.sort(reverse=True)  # high to low

print(f"{'rank':>4s} {'pair':22s} {'W_ad,paper max':>14s} {'W_ad,ours max':>14s} {'Δ avg':>8s}")
for i, (w, pair, rs) in enumerate(ranked_paper):
    w_ours_max = max(r["Wad_ours"] for r in rs)
    avg_delta = sum(r["strain_correction"] for r in rs) / len(rs)
    print(f"{i+1:>4d} {pair:22s} {w:>14.3f} {w_ours_max:>14.3f} {avg_delta:>+8.3f}")

(RES / "wad_paper_style.json").write_text(json.dumps({
    "natural_beta_slabs": results_natural,
    "table": table,
    "ranking_paper_style_max": [{"rank": i+1, "pair": p, "Wad_paper_max": w}
                                for i, (w, p, _) in enumerate(ranked_paper)],
}, indent=2))
print(f"\nSaved: {RES / 'wad_paper_style.json'}")
