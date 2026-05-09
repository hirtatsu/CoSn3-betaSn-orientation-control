# W_ad α-CoSn₃ / β-Sn (PFP/PBE)

Source for **Fig. 9** (main) and **Table S5** (supplement).

## 3 paper interface pairs

| Pair | termination | W_ad (J/m²) | Paper role |
|---|---|---:|---|
| α(600) / β(100) | term1 (Co/Sn) | **2.28** | hero (c-axis ∥ substrate) |
| α(600) / β(001) | term1 (Co/Sn) | **1.89** | c-axis ⊥ contrast |
| α(312) / β(100) | term0 (Co₂Sn₃/Sn) | **0.82** | RF (random) control |

(Values used in main text correspond to "energetically favored termination on each side" — see Table S5 for full termination dependence.)

## Files

| Path | Content |
|---|---|
| `stacks_input/{pair}_term*_interface.cif` | 5 interface input CIFs (initial atom-only) |
| `stacks_relaxed/{pair}_term*_interface_cellrelax_PBE.cif` | 5 cell-relaxed interfaces |
| `slabs_split_input/{pair}_term*_alpha.cif`, `*_beta.cif` | 10 split slab inputs |
| `slabs_split_relaxed/{pair}_term*_alpha_cellrelax.cif` | α slabs cell-relaxed (paper-style) |
| `slabs_split_relaxed/{pair}_beta_natural.cif` | β natural-cell slabs (one per pair, terminator-independent) |
| `slab_split_index.json` | termination labels + cell parameters |
| `interface_cellrelax_PBE.json` | E_interface for each stack |
| `alpha_ifc_cellrelax_PBE.json` | E_α-slab for each stack |
| `beta_natural_cellrelax_PBE.json` | E_β-natural for each pair |
| **`wad_paper_full_consistent.json`** ★ | **W_ad master (Eq. 3, paper-style with cell relaxation)** |
| `scripts/19,20,21_*.py` | PFP relax + W_ad calculation |

## Calculation procedure (paper-style)

1. Build interface stack via pymatgen `CoherentInterfaceBuilder` (ZSL match)
2. Cell-relax interface (FrechetCellFilter mask=[T,T,F,F,F,T])
3. Extract α-only and β-only sides → cell-relax independently (each finds its natural equilibrium cell)
4. W_ad = (E_α-slab + E_β-slab − E_int) / A_int

## Two-tier discrimination (paper claim)

- **Tier 1** (c-axis selection on α(600)): β(100) > β(001) by 0.39 J/m² (~17%) — **β c-axis ∥ substrate is favored**
- **Tier 2** (α-face selection on β(100)): α(600) > α(312) by factor ~2.8 — **α(600) is the effective seed plane**
