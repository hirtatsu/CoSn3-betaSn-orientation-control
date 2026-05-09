# Manuscript data → file map

For each paper element, the table below lists the specific files reproducing it.

## Main text

### Table 1 — α-CoSn₃ surface energies (MLIP)
| Plane | γ_PFP (J/m²) | Source |
|---|---:|---|
| (600) | 0.46 | `pfp/alpha_surface_energies/alpha_surfaces_cellrelax_PBE.json` (key: `alpha_600_t*` → min γ) |
| (301) | 0.66 | (same JSON, `alpha_301_t*`) |
| (321) | 0.73 | `pfp/alpha_surface_energies/alpha_321_cellrelax/alpha_321_gamma_PBE.json` |
| (010) | 0.77 | `pfp/alpha_surface_energies/alpha_surfaces_cellrelax_PBE.json` (`alpha_010_t*`) |

Inputs: `pfp/alpha_surface_energies/slabs_input/*.cif` (12 files, 4 face × 3 termination)
Relaxed: `pfp/alpha_surface_energies/slabs_relaxed/*.cif`

### Fig. 6 — Si/α-CoSn₃ lattice matching (Co overlay)
- Bulk lattices for spacing calculation: `pfp/bulks/{alpha_CoSn3,Si}_PBE.cif`
- Disregistry analysis (Si Co): see Table 2 / Table S6

### Table 2 — Disregistry α/β
- Lattice constants from `pfp/bulks/{alpha_CoSn3,beta_Sn}_PBE.cif`
- Computed values (7.48%, 1.48%, 7.48%, 20.16%, 10.12%) computed analytically from these

### Fig. 7 — W_ad Si/α-CoSn₃
| Pair | W_ad (J/m²) | Source |
|---|---:|---|
| Si(100)/α(600) Co-term | **1.06** | `pfp/Si_alpha_Wad/wad_Si_alpha_PBE.json` (key: `Si100_alpha600_term1`) |
| Si(100)/α(312) | **2.57** | (same JSON, `Si100_alpha312_term0`) |

Structures (relaxed): `pfp/Si_alpha_Wad/stacks_relaxed/Si100_alpha600_term1_int.cif`, `Si100_alpha312_term0_int.cif`

### Fig. 9 — W_ad α-CoSn₃/β-Sn
| Pair | W_ad PFP | W_ad DFT | Source |
|---|---:|---:|---|
| α(600)/β(100) | **2.28** | **2.41** | PFP: `pfp/alpha_beta_Wad/wad_paper_full_consistent.json`; DFT: `dft/wad_dft_result.json` (derived from `dft/Done_CoSn3(600)+Sn(100)_Phase1{A,B,C}/`) |
| α(600)/β(001) | **1.89** | — | PFP: same JSON, term1 |
| α(312)/β(100) | **0.82** | — | PFP: same JSON, term0 |

Structures (relaxed): `pfp/alpha_beta_Wad/stacks_relaxed/*_term*_interface_cellrelax_PBE.cif`

## Supplementary text

### Table S1 — Bulk lattice constants
- α-CoSn₃: `pfp/bulks/bulk_alpha_CoSn3_PBE.{cif,json}`, `dft/Done_Bulk_CoSn3/test.cif`
- β-Sn: `pfp/bulks/bulk_beta_Sn_PBE.{cif,json}` (DFT not computed in this work)
- Si: `pfp/bulks/Si_PBE.{cif,json}` (DFT not computed)

### Table S2 — α-CoSn₃ γ comparison
- γ_PFP: `pfp/alpha_surface_energies/alpha_surfaces_cellrelax_PBE.json`
- γ_DFT (600): `dft/Done_Surf_CoSn3(600)/test.ene` + `dft/Done_Bulk_CoSn3/test.ene` → **0.545 J/m²**
- γ_DFT (010): `dft/Done_Surf_CoSn3(010)/test.ene` + `dft/Done_Bulk_CoSn3/test.ene` → **0.878 J/m²**
- γ_DFT (301): see `dft/gamma_DFT.json` (`results[2]`, γ = 0.685 J/m²)
- All γ_DFT collected in `dft/gamma_DFT.json`
- γ_Ref (Ma 2020 [26], Wang 2024 [27]): from cited literature

### Table S3 — Si/α termination dependence W_ad
- All termination data: `pfp/Si_alpha_Wad/wad_Si_alpha_PBE.json` (`Si100_alpha600_term0`, `term1`, `Si100_alpha312_term0`)
- Strain values: `pfp/Si_alpha_Wad/build_index.json`

### Table S4 — α(600)/β(100) PFP vs DFT W_ad
- PFP: `pfp/alpha_beta_Wad/wad_paper_full_consistent.json` (`alpha600_beta100_term1`)
- DFT: `dft/wad_dft_result.json` (W_ad = 2.41 J/m²); raw inputs in `dft/Done_CoSn3(600)+Sn(100)_Phase1{A,B,C}/`
- Interface-relaxation cost (0.609 J/m²) computed from Phase1A (relaxed) vs Phase2A (rigid stack on relaxed-AB geometry)

### Table S5 — α/β termination dependence W_ad
- All values: `pfp/alpha_beta_Wad/wad_paper_full_consistent.json` (filter to 3 paper pairs)

### Fig. S3 — Sn-atom matching diagrams
- Bulk lattices for spacing calculation: `pfp/bulks/{alpha_CoSn3,beta_Sn}_PBE.cif`
- Disregistry analysis (Sn-Sn): see Table S6

### Table S6 — Disregistry Si/α (Co and Sn)
- Lattice constants from `pfp/bulks/{alpha_CoSn3,Si}_PBE.cif` and experimental Si (5.431 Å)
- Disregistry formulas in supplement S3.2

### Fig. S4 — plane-on-plane vs edge-to-edge schematic
- Schematic only, no calculation data

## Equations

### Eq. 2 — γ definition
```
γ = (E_slab − N · μ) / (2 · A)
```
- E_slab from JSON files (e.g., `pfp/alpha_surface_energies/alpha_surfaces_cellrelax_PBE.json`)
- μ from `pfp/bulks/bulk_alpha_CoSn3_PBE.json`

### Eq. 3 — W_ad definition
```
W_ad = (E_A_slab + E_B_slab − E_AB_interface) / A
```
- All energies from `pfp/{alpha,Si}_*_Wad/*.json` master files

## Computational methods reference

For full procedure, see Tatsumi et al. *Acta Mater.* **304**, 121813 (2026) [referenced in main text].
