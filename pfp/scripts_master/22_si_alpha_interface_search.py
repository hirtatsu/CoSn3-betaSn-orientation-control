"""Search Si(100)/α-CoSn3(600) and Si(100)/α-CoSn3(312) interface supercells.

Following Tatsumi 2026 Acta Mater convention (paper-style W_ad).
Substrate = Si (rigid, matches experiment: Si bottom, α deposited on top)
Film = α-CoSn3 (strained to Si)

Output:
  pfp/Si_alpha_wad/interface_candidates.json
  pfp/Si_alpha_wad/Si_alpha600_term*.cif
  pfp/Si_alpha_wad/Si_alpha312_term*.cif
  pfp/Si_alpha_wad/Si_alpha600_term*_{Si,alpha,interface}.cif (split slabs)
  pfp/Si_alpha_wad/Si_alpha312_term*_{Si,alpha,interface}.cif

Settings: ZSL max_strain=10%, max_area=400 Å² (matching previous W_ad scripts)
"""
from pathlib import Path
import numpy as np
import json
from pymatgen.io.cif import CifParser
from pymatgen.core import Structure
from pymatgen.analysis.interfaces.coherent_interfaces import CoherentInterfaceBuilder
from pymatgen.analysis.interfaces.zsl import ZSLGenerator

DATA_DIR = Path("/home/tatsumi/projects/cosn3-betasn-adhesion/data")
PFP_DIR = Path("/home/tatsumi/projects/cosn3-betasn-adhesion/pfp")
OUT_DIR = PFP_DIR / "Si_alpha_wad"
OUT_DIR.mkdir(exist_ok=True)
SLABS_OUT = OUT_DIR / "interface_stacks"
SLABS_OUT.mkdir(exist_ok=True)

# Si: use experimental bulk lattice (we'll PFP-relax separately)
si = CifParser(DATA_DIR.parent / "tmp" / "silicon Fd-3m (227)-2104737.cif").parse_structures(primitive=False)[0] \
    if (DATA_DIR.parent / "tmp" / "silicon Fd-3m (227)-2104737.cif").exists() \
    else CifParser(DATA_DIR.parent / "tmp" / "silicon Fd-3m (227)-2104737.cif").parse_structures(primitive=False)[0]
print(f"Si bulk: {si.formula}, N={len(si)}, cell {si.lattice.abc}")

# α-CoSn3: use PFP-relaxed bulk
acs = CifParser(PFP_DIR / "bulk_alpha_CoSn3_PBE.cif").parse_structures(primitive=False)[0]
print(f"α-CoSn3 PFP-relaxed: N={len(acs)}, cell {acs.lattice.abc}")

PAIRS = [
    ("Si100_alpha600",  (1,0,0), (6,0,0), "manuscript: Si<011>/α<011> 13.38%"),
    ("Si100_alpha312",  (1,0,0), (3,1,2), "manuscript: Si<001>/α<112> 0.99%"),
]

zsl = ZSLGenerator(
    max_area_ratio_tol=0.10, max_area=400.0,
    max_length_tol=0.10, max_angle_tol=0.05,
)

print()
print("=" * 110)
print(f"{'pair':22s} {'rank':>4s} {'area':>7s} {'sub axes':>14s} {'film axes':>14s} "
      f"{'εa%':>6s} {'εb%':>6s} {'Δang':>6s}  comment")
print("=" * 110)

candidates = {}
for tag, sub_hkl, film_hkl, note in PAIRS:
    cib = CoherentInterfaceBuilder(
        substrate_structure=si, film_structure=acs,
        substrate_miller=sub_hkl, film_miller=film_hkl, zslgen=zsl,
    )
    matches = cib.zsl_matches
    print(f"\n[{tag}] note: {note}, total ZSL matches: {len(matches)}")
    described = []
    for idx, m in enumerate(matches):
        s_a = np.linalg.norm(m.substrate_sl_vectors[0])
        s_b = np.linalg.norm(m.substrate_sl_vectors[1])
        f_a = np.linalg.norm(m.film_sl_vectors[0])
        f_b = np.linalg.norm(m.film_sl_vectors[1])
        cos_s = np.dot(m.substrate_sl_vectors[0], m.substrate_sl_vectors[1]) / (s_a*s_b)
        cos_f = np.dot(m.film_sl_vectors[0], m.film_sl_vectors[1]) / (f_a*f_b)
        ang_s = np.degrees(np.arccos(np.clip(cos_s, -1, 1)))
        ang_f = np.degrees(np.arccos(np.clip(cos_f, -1, 1)))
        u_a = (f_a - s_a) / s_a * 100
        u_b = (f_b - s_b) / s_b * 100
        described.append({
            "area": float(m.match_area),
            "sub_a": float(s_a), "sub_b": float(s_b), "sub_ang": float(ang_s),
            "film_a": float(f_a), "film_b": float(f_b), "film_ang": float(ang_f),
            "strain_a_pct": float(u_a), "strain_b_pct": float(u_b),
            "ang_mismatch_deg": float(ang_f - ang_s),
            "match_obj_idx": idx,
        })
    # Sort by area then by max strain
    ranked = sorted(described, key=lambda d: (d["area"], max(abs(d["strain_a_pct"]), abs(d["strain_b_pct"]))))
    # Dedup
    dedup = []
    for d in ranked:
        if not any(abs(d["area"]-e["area"]) < 0.5 and abs(d["strain_a_pct"]-e["strain_a_pct"]) < 0.3
                   and abs(d["strain_b_pct"]-e["strain_b_pct"]) < 0.3 for e in dedup):
            dedup.append(d)
    top5 = dedup[:5]
    candidates[tag] = top5

    if not top5:
        print(f"  {tag}: no match within tolerance!")
        continue
    for i, d in enumerate(top5):
        print(f"{tag if i==0 else '':22s} {i+1:>4d} {d['area']:>7.1f} "
              f"{d['sub_a']:>5.2f}x{d['sub_b']:>5.2f}@{d['sub_ang']:>3.0f} "
              f"{d['film_a']:>5.2f}x{d['film_b']:>5.2f}@{d['film_ang']:>3.0f} "
              f"{d['strain_a_pct']:>+6.2f} {d['strain_b_pct']:>+6.2f} "
              f"{d['ang_mismatch_deg']:>+6.2f}  {note if i==0 else ''}")

# Save candidate enumeration
out_json = OUT_DIR / "interface_candidates.json"
out_json.write_text(json.dumps(candidates, indent=2))
print(f"\nSaved candidate enumeration: {out_json}")

# Now try to actually BUILD the rank-1 interfaces and report atom counts
print("\n" + "=" * 110)
print("Building rank-1 interfaces (atom count check; CIF output deferred to user approval)")
print("=" * 110)
print(f"{'pair':22s} {'N_int':>6s} {'N_Si':>5s} {'N_α':>5s} {'a':>6s} {'b':>6s} {'γ':>5s} {'Lz':>6s} {'#term':>5s}")
print("-" * 110)
for tag, sub_hkl, film_hkl, note in PAIRS:
    if not candidates[tag]:
        continue
    cib = CoherentInterfaceBuilder(
        substrate_structure=si, film_structure=acs,
        substrate_miller=sub_hkl, film_miller=film_hkl, zslgen=zsl,
    )
    chosen = cib.zsl_matches[candidates[tag][0]["match_obj_idx"]]
    cib.zsl_matches = [chosen]
    n_terms = len(cib.terminations)
    try:
        gens = cib.get_interfaces(
            termination=cib.terminations[0],
            substrate_thickness=14.0, film_thickness=14.0,
            vacuum_over_film=20.0, gap=2.5, in_layers=False,
        )
        intf = next(gens)
        labels = intf.site_properties.get("interface_label", [])
        n_si = sum(1 for x in labels if "substrate" in str(x))
        n_a = sum(1 for x in labels if "film" in str(x))
        a, b, c = intf.lattice.abc
        gamma = intf.lattice.gamma
        print(f"{tag:22s} {len(intf):>6d} {n_si:>5d} {n_a:>5d} "
              f"{a:>6.2f} {b:>6.2f} {gamma:>5.1f} {c:>6.1f} {n_terms:>5d}")
    except Exception as e:
        print(f"{tag:22s} ERROR: {e}")

