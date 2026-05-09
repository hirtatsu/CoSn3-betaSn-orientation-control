"""Generate alpha-CoSn3 slab CIFs locally with pymatgen SlabGenerator (all terminations).
Output: results/alpha_slabs/<face>_t<n>.cif + index json.
"""
import json
from pathlib import Path
from math import gcd
from pymatgen.io.cif import CifParser
from pymatgen.core.surface import SlabGenerator
from pymatgen.io.ase import AseAtomsAdaptor

RES = Path("/home/tatsumi/projects/cosn3-betasn-adhesion/results")
SLABS = RES / "alpha_slabs"; SLABS.mkdir(exist_ok=True)

bulk = CifParser(RES / "bulk_alpha_CoSn3_PBE.cif").parse_structures(primitive=False)[0]

REF_GAMMA = {(6,0,0): 0.46, (0,1,0): 1.05, (3,2,1): 1.09, (3,0,1): 0.78, (3,1,2): None}
FACES = [(6,0,0), (3,1,2), (0,1,0), (3,2,1), (3,0,1)]

def reduce_hkl(hkl):
    g = gcd(gcd(abs(hkl[0]), abs(hkl[1])), abs(hkl[2])) or 1
    return (hkl[0]//g, hkl[1]//g, hkl[2]//g)

index = {}
print("=" * 80)
print(f"{'face':>10s} {'reduced':>10s} {'#terms':>7s} {'thickness':>10s} {'A [A^2]':>9s} {'N range':>10s}")
print("=" * 80)

for hkl in FACES:
    rhkl = reduce_hkl(hkl)
    sg = SlabGenerator(bulk, miller_index=rhkl,
                       min_slab_size=14.0, min_vacuum_size=15.0,
                       center_slab=True, in_unit_planes=False,
                       primitive=False, max_normal_search=2)
    slabs = sg.get_slabs(symmetrize=False)
    n_atoms_list = [len(s) for s in slabs]
    thickness = max(s.cart_coords[:,2].max() - s.cart_coords[:,2].min() for s in slabs) if slabs else 0
    a, b, _ = slabs[0].lattice.abc if slabs else (0,0,0)
    import numpy as np
    gamma_deg = slabs[0].lattice.gamma if slabs else 90
    area = a*b*np.sin(np.radians(gamma_deg))
    print(f"{str(hkl):>10s} {str(rhkl):>10s} {len(slabs):>7d} {thickness:>10.2f} {area:>9.2f} {min(n_atoms_list)}-{max(n_atoms_list)}")

    rec_list = []
    for ti, slab in enumerate(slabs):
        # Convert to ASE Atoms (preserve cell), write CIF
        ase_at = AseAtomsAdaptor.get_atoms(slab)
        tag = f"alpha_{''.join(map(str,hkl))}_t{ti}"
        cif_path = SLABS / f"{tag}.cif"
        from ase.io import write as ase_write
        ase_write(cif_path, ase_at)
        rec_list.append({
            "tag": tag, "termination_idx": ti, "n_atoms": len(slab),
            "lattice_a": float(slab.lattice.a), "lattice_b": float(slab.lattice.b),
            "lattice_c": float(slab.lattice.c), "gamma_deg": float(slab.lattice.gamma),
            "area_A2": float(area), "cif": cif_path.name,
        })
    index[str(hkl)] = {"hkl_nominal": list(hkl), "hkl_reduced": list(rhkl),
                       "n_terminations": len(slabs),
                       "gamma_ref_Table1": REF_GAMMA.get(hkl),
                       "terminations": rec_list}

(RES / "alpha_slabs_index.json").write_text(json.dumps(index, indent=2))
print(f"\nSaved: {RES / 'alpha_slabs_index.json'}")
print(f"CIFs in: {SLABS}/  ({sum(len(v['terminations']) for v in index.values())} files)")
