# α-CoSn₃ surface energies (PFP/PBE)

Source for **Table 1** (main) and **Table S2** (supplement) MLIP values.

## Procedure

1. Generate slabs with pymatgen `SlabGenerator(symmetrize=False, ftol=0.1)` for full termination enumeration
2. Cell-relax each slab (FrechetCellFilter mask=[T,T,F,F,F,T], fmax = 0.015 eV/Å)
3. Compute γ = (E_slab − N · μ_α) / (2 · A) using Eq. 2 of main text

μ_α = E_bulk(Co₈Sn₂₄) / 32 = `bulks/bulk_alpha_CoSn3_PBE.json` → "energy_eV" / 32

## Files

| Path | Content |
|---|---|
| `slabs_input/alpha_{600,010,301}_t{0,1,2}.cif` | 9 input CIFs (3 face × 3 termination) |
| `slabs_input/alpha_321_t{0,1,2}.cif` | 3 input CIFs for (321) |
| `alpha_slabs_index.json` | termination labels and metadata |
| `slabs_relaxed/alpha_{600,010,301}_t{0,1,2}.cif` | 9 cell-relaxed CIFs |
| **`alpha_surfaces_cellrelax_PBE.json`** ★ | **γ master (12 entries)** |
| `alpha_321_cellrelax/` | (321) PFP relax data (864-atom primitive slabs × 3 terminations) |
| `scripts/16b_local_alpha_slab_gen.py` | Slab generation script |
| `scripts/17_alpha_surface_relax.py` | PFP relax + γ calculation |

## Key numerical values (per Table 1)

| Plane | γ_PFP min (J/m²) | γ_Ma 2020 (J/m²) |
|---|---:|---:|
| (600) | 0.459 | 0.46 |
| (301) | 0.660 | 0.78 |
| (321) | 0.731 | 1.09 |
| (010) | 0.770 | 1.05 |

Min over symmetrically distinct terminations is reported (Table 1, S2).
Max-γ termination is also of interest (interface reconstruction discussion in Supplement).
