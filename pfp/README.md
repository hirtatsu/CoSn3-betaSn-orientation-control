# pfp/ — MLIP (PFP) calculation data

All MLIP / Matlantis results for the paper.

## Subfolders

| Folder | Content | Paper element |
|---|---|---|
| `bulks/` | α-CoSn₃, β-Sn, Si bulk relaxed structures + JSONs | Table S1 |
| `alpha_surface_energies/` | γ for 4 faces × 3 terminations | Table 1, Table S2 |
| `alpha_beta_Wad/` | W_ad for 3 α/β interfaces (interface + α slab + β natural slab) | Fig. 9, Table S5 |
| `Si_alpha_Wad/` | W_ad for Si/α(600) and Si/α(312) | Fig. 7, Table S3 |
| `scripts_master/` | All key reproducibility scripts (numbered execution order) | — |

## Key master JSONs

```bash
# γ master (4 faces × 3 termination, sorted by face)
pfp/alpha_surface_energies/alpha_surfaces_cellrelax_PBE.json

# W_ad α/β master (10 stacks; paper uses 3)
pfp/alpha_beta_Wad/wad_paper_full_consistent.json

# W_ad Si/α master
pfp/Si_alpha_Wad/wad_Si_alpha_PBE.json
```

## Computational procedure

Bulk relaxation (ExpCellFilter, full cell + atoms, fmax = 0.001 eV/Å) →
Slab generation (pymatgen SlabGenerator, full termination enumeration, vacuum ≈ 15 Å) →
Slab cell relaxation (FrechetCellFilter mask=[T,T,F,F,F,T], fmax = 0.015 eV/Å) →
γ or W_ad evaluation per Eq. 2/3 of main text.

Paper-style W_ad recipe (Tatsumi 2026 *Acta Mater.* **304**, 121813 Eq. 4):
- E_int: cell-relaxed interface
- E_α-slab: cell-relaxed α slab (independent, naturally relaxes to α equilibrium)
- E_β-slab: cell-relaxed β slab (β natural cell)
- W_ad = (E_α + E_β − E_int) / A_int
