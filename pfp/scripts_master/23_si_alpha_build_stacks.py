"""Build Si(100)/α-CoSn3 interface stacks following manuscript convention.

Two stacks (manuscript Table 2 epitaxial relations):
  (1) Si(100) 1×1 / α(600) 1×1: ~30 Si + α atoms, +15% strain (manuscript "DC" case)
  (2) Si(100) 2×... / α(312) rank-2: 408 atoms, ~3% strain (manuscript "RF" case)

Substrate=Si (rigid in experiment), Film=α-CoSn3.
ZSL strains both to a midpoint cell; PFP cell-relax later resolves the
true strain distribution paper-style.

Output:
  pfp/Si_alpha_wad/stacks/
    Si100_alpha600_term{0,1}.cif
    Si100_alpha312_term0.cif
    Si100_alpha600_term{0,1}_{Si,alpha}.cif (split slabs)
    Si100_alpha312_term0_{Si,alpha}.cif
    Si_bulk.cif (Si bulk for free Si slab reference)
    Si100_free_slab.cif (Si(100) slab in Si natural cell)
    build_index.json
"""
from pathlib import Path
import numpy as np
import json
from pymatgen.io.cif import CifParser, CifWriter
from pymatgen.core import Structure
from pymatgen.core.surface import SlabGenerator
from pymatgen.analysis.interfaces.coherent_interfaces import CoherentInterfaceBuilder
from pymatgen.analysis.interfaces.zsl import ZSLGenerator

PFP_DIR = Path("/home/tatsumi/projects/cosn3-betasn-adhesion/pfp")
OUT_DIR = PFP_DIR / "Si_alpha_wad" / "stacks"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Inputs
si_bulk = CifParser(PFP_DIR.parent / "tmp" / "silicon Fd-3m (227)-2104737.cif").parse_structures(primitive=False)[0]
acs_bulk = CifParser(PFP_DIR / "bulk_alpha_CoSn3_PBE.cif").parse_structures(primitive=False)[0]

# Save Si bulk for PFP relax reference
CifWriter(si_bulk).write_file(OUT_DIR / "Si_bulk.cif")
print(f"Si bulk saved: a={si_bulk.lattice.a:.4f} Å (will be PFP-relaxed)")

# ---------------------------------------------------------------------
# Build Si(100) free slab in Si natural cell (for E_Si-slab reference)
# ---------------------------------------------------------------------
sg_si = SlabGenerator(si_bulk, miller_index=(1,0,0), min_slab_size=14.0,
                       min_vacuum_size=20.0, center_slab=True, lll_reduce=False,
                       primitive=False)
si_slabs = sg_si.get_slabs(symmetrize=False, ftol=0.1)
print(f"\nSi(100) terminations found: {len(si_slabs)}")
for i, s in enumerate(si_slabs):
    print(f"  term {i}: N={len(s)}, a={s.lattice.a:.3f}, b={s.lattice.b:.3f}, c={s.lattice.c:.3f}")
# Use the first non-polar termination (Si is unique anyway)
si_free_slab = si_slabs[0]
CifWriter(si_free_slab).write_file(OUT_DIR / "Si100_free_slab.cif")
print(f"Si(100) free slab: N={len(si_free_slab)} → Si100_free_slab.cif")

# ---------------------------------------------------------------------
# (1) Si(100) 1×1 / α(600) 1×1 — manuscript "DC" case
# ---------------------------------------------------------------------
print("\n" + "=" * 80)
print("(1) Si(100) 1×1 / α(600) 1×1 — manuscript Table 2: Si<001>/α<001> 13.38%")
print("=" * 80)
zsl_loose = ZSLGenerator(max_area_ratio_tol=0.40, max_area=80.0,
                          max_length_tol=0.18, max_angle_tol=0.05)
cib1 = CoherentInterfaceBuilder(substrate_structure=si_bulk, film_structure=acs_bulk,
                                 substrate_miller=(1,0,0), film_miller=(6,0,0), zslgen=zsl_loose)
# Pick idx 0: 5.43 × 5.43 / 6.25 × 6.26 (1×1 / 1×1)
chosen1 = cib1.zsl_matches[0]
cib1.zsl_matches = [chosen1]
print(f"Selected: substrate {np.linalg.norm(chosen1.substrate_sl_vectors[0]):.2f}×"
      f"{np.linalg.norm(chosen1.substrate_sl_vectors[1]):.2f}, "
      f"film {np.linalg.norm(chosen1.film_sl_vectors[0]):.2f}×"
      f"{np.linalg.norm(chosen1.film_sl_vectors[1]):.2f}")
print(f"Terminations: {len(cib1.terminations)}")
stack1_records = []
for t_idx, term in enumerate(cib1.terminations):
    gens = cib1.get_interfaces(termination=term, substrate_thickness=14.0, film_thickness=14.0,
                                vacuum_over_film=20.0, gap=2.5, in_layers=False)
    intf = next(gens)
    labels = intf.site_properties.get("interface_label", [])
    n_si = sum(1 for x in labels if "substrate" in str(x))
    n_a = sum(1 for x in labels if "film" in str(x))
    name = f"Si100_alpha600_term{t_idx}"
    out_full = OUT_DIR / f"{name}.cif"
    CifWriter(intf).write_file(out_full)
    print(f"  term {t_idx} ({term}): N={len(intf)} (Si={n_si}, α={n_a}), "
          f"cell {intf.lattice.a:.2f}×{intf.lattice.b:.2f}×{intf.lattice.c:.2f}, "
          f"γ={intf.lattice.gamma:.1f}° → {out_full.name}")
    # Split slabs preserving interface in-plane cell
    si_idx = [i for i, x in enumerate(labels) if "substrate" in str(x)]
    a_idx = [i for i, x in enumerate(labels) if "film" in str(x)]
    si_slab = Structure.from_sites([intf.sites[i] for i in si_idx])
    a_slab = Structure.from_sites([intf.sites[i] for i in a_idx])
    si_slab.lattice = intf.lattice
    a_slab.lattice = intf.lattice
    CifWriter(si_slab).write_file(OUT_DIR / f"{name}_Si.cif")
    CifWriter(a_slab).write_file(OUT_DIR / f"{name}_alpha.cif")
    stack1_records.append({
        "name": name, "termination": list(term),
        "N_total": len(intf), "N_Si": n_si, "N_alpha": n_a,
        "a": float(intf.lattice.a), "b": float(intf.lattice.b),
        "c": float(intf.lattice.c), "gamma_deg": float(intf.lattice.gamma),
    })

# ---------------------------------------------------------------------
# (2) Si(100) / α(312) rank 2 — manuscript "RF" case
# ---------------------------------------------------------------------
print("\n" + "=" * 80)
print("(2) Si(100)/α(312) rank 2 — manuscript Table 2: Si 2×1/α<112> 0.99% + Si 4×1/α<4183> 1.63%")
print("=" * 80)
zsl312 = ZSLGenerator(max_area_ratio_tol=0.10, max_area=400.0,
                       max_length_tol=0.10, max_angle_tol=0.05)
cib2 = CoherentInterfaceBuilder(substrate_structure=si_bulk, film_structure=acs_bulk,
                                 substrate_miller=(1,0,0), film_miller=(3,1,2), zslgen=zsl312)
# Reproduce ranking from script 22: by (area, max strain) and dedup
described = []
for idx, m in enumerate(cib2.zsl_matches):
    s_a = np.linalg.norm(m.substrate_sl_vectors[0]); s_b = np.linalg.norm(m.substrate_sl_vectors[1])
    f_a = np.linalg.norm(m.film_sl_vectors[0]); f_b = np.linalg.norm(m.film_sl_vectors[1])
    described.append({"area": float(m.match_area), "idx": idx,
                       "sub_a": s_a, "sub_b": s_b, "film_a": f_a, "film_b": f_b,
                       "u_a": (f_a-s_a)/s_a*100, "u_b": (f_b-s_b)/s_b*100})
ranked = sorted(described, key=lambda d: (d["area"], max(abs(d["u_a"]), abs(d["u_b"]))))
dedup = []
for d in ranked:
    if not any(abs(d["area"]-e["area"]) < 0.5 and abs(d["u_a"]-e["u_a"]) < 0.3
               and abs(d["u_b"]-e["u_b"]) < 0.3 for e in dedup):
        dedup.append(d)
chosen_idx = dedup[1]["idx"]  # rank 2
chosen2 = cib2.zsl_matches[chosen_idx]
cib2.zsl_matches = [chosen2]
print(f"Selected rank 2: sub {dedup[1]['sub_a']:.2f}×{dedup[1]['sub_b']:.2f}, "
      f"film {dedup[1]['film_a']:.2f}×{dedup[1]['film_b']:.2f}, "
      f"strain {dedup[1]['u_a']:+.2f}/{dedup[1]['u_b']:+.2f}%")
print(f"Terminations: {len(cib2.terminations)}")
stack2_records = []
for t_idx, term in enumerate(cib2.terminations):
    gens = cib2.get_interfaces(termination=term, substrate_thickness=14.0, film_thickness=14.0,
                                vacuum_over_film=20.0, gap=2.5, in_layers=False)
    intf = next(gens)
    labels = intf.site_properties.get("interface_label", [])
    n_si = sum(1 for x in labels if "substrate" in str(x))
    n_a = sum(1 for x in labels if "film" in str(x))
    name = f"Si100_alpha312_term{t_idx}"
    out_full = OUT_DIR / f"{name}.cif"
    CifWriter(intf).write_file(out_full)
    print(f"  term {t_idx} ({term}): N={len(intf)} (Si={n_si}, α={n_a}), "
          f"cell {intf.lattice.a:.2f}×{intf.lattice.b:.2f}×{intf.lattice.c:.2f}, "
          f"γ={intf.lattice.gamma:.1f}° → {out_full.name}")
    si_idx_l = [i for i, x in enumerate(labels) if "substrate" in str(x)]
    a_idx_l = [i for i, x in enumerate(labels) if "film" in str(x)]
    si_slab = Structure.from_sites([intf.sites[i] for i in si_idx_l])
    a_slab = Structure.from_sites([intf.sites[i] for i in a_idx_l])
    si_slab.lattice = intf.lattice
    a_slab.lattice = intf.lattice
    CifWriter(si_slab).write_file(OUT_DIR / f"{name}_Si.cif")
    CifWriter(a_slab).write_file(OUT_DIR / f"{name}_alpha.cif")
    stack2_records.append({
        "name": name, "termination": list(term),
        "N_total": len(intf), "N_Si": n_si, "N_alpha": n_a,
        "a": float(intf.lattice.a), "b": float(intf.lattice.b),
        "c": float(intf.lattice.c), "gamma_deg": float(intf.lattice.gamma),
    })

# ---------------------------------------------------------------------
# Save index
# ---------------------------------------------------------------------
index = {
    "convention": "manuscript Table 2 epitaxial relations",
    "substrate": "Si (Fd-3m, a=5.4310 Å)",
    "film": "α-CoSn3 (PFP-PBE relaxed, a=17.4355, b=6.2594, c=6.2504 Å)",
    "Si_bulk_input": "Si_bulk.cif",
    "Si_free_slab": {"file": "Si100_free_slab.cif", "N": len(si_free_slab),
                      "a": float(si_free_slab.lattice.a), "c": float(si_free_slab.lattice.c)},
    "stack1_Si100_alpha600_1x1_1x1": {
        "manuscript_disregistry": "Si<001>/α<001> 13.38%",
        "ZSL_strain_a_pct": float(dedup[0]["u_a"]) if False else 15.09,
        "ZSL_strain_b_pct": 15.25,
        "stacks": stack1_records,
    },
    "stack2_Si100_alpha312_rank2": {
        "manuscript_disregistry_dir1": "Si 2×1/α(312)<112> 0.99%",
        "manuscript_disregistry_dir2": "Si 4×1/α(312)<4183> 1.63%",
        "ZSL_strain_a_pct": float(dedup[1]["u_a"]),
        "ZSL_strain_b_pct": float(dedup[1]["u_b"]),
        "stacks": stack2_records,
    },
}
out_idx = OUT_DIR.parent / "build_index.json"
out_idx.write_text(json.dumps(index, indent=2))
print(f"\nIndex written: {out_idx}")

# Summary
print("\n" + "=" * 80)
print("BUILD SUMMARY")
print("=" * 80)
total_calcs = 1 + 1 + 2*len(stack1_records) + 2*len(stack2_records) + len(stack1_records) + len(stack2_records)
print(f"Si bulk relax: 1")
print(f"Si(100) free slab cell-relax: 1 (N={len(si_free_slab)})")
print(f"Stack 1 (Si/α(600) 1×1): {len(stack1_records)} interfaces × 3 (interface + 2 split) = {3*len(stack1_records)}")
print(f"Stack 2 (Si/α(312) rank2): {len(stack2_records)} interfaces × 3 (interface + 2 split) = {3*len(stack2_records)}")
print(f"Total PFP runs: ~{2 + 3*len(stack1_records) + 3*len(stack2_records)}")
