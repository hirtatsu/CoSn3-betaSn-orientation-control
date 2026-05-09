# W_ad Si(100) / α-CoSn₃ (PFP/PBE)

Source for **Fig. 7** (main) and **Table S3** (supplement).

## Two interface pairs (manuscript convention)

| Pair | Cell | Strain | W_ad (J/m²) | Paper role |
|---|---|---|---:|---|
| Si(100) / α(600) | 1×1 / 1×1 | Si +10%, α −4% | **+1.06** (Co-term) / −0.24 (Sn-term) | DC condition (13.38% mismatch) |
| Si(100) / α(312) | rank-2 supercell | Si −0.8%, α +2.9% | **+2.57** | RF condition (0.99% mismatch) |

Manuscript convention: Si(100)/α(600) uses the 1×1/1×1 stack to retain the 13.38% disregistry geometry of Section 4 (Fig. 6); Si(100)/α(312) preserves the 0.99% Si⟨001⟩ ∥ α⟨11̄2⟩ alignment.

## Files

| Path | Content |
|---|---|
| `stacks_input/Si_bulk.cif`, `Si100_free_slab.cif` | Si bulk + free Si(100) slab inputs |
| `stacks_input/Si100_alpha600_term{0,1}.cif` | 2 input interface CIFs (Sn-term, Co-term) |
| `stacks_input/Si100_alpha312_term0.cif` | 1 input interface CIF (Co₂Sn₃ term) |
| `stacks_input/*_{Si,alpha}.cif` | Split slabs (Si side, α side) |
| `stacks_relaxed/*.cif` | 11 cell-relaxed structures (3 interfaces + 6 split slabs + Si bulk + Si free slab) |
| `Si_bulk_tight.json` | Si bulk PFP relax (fmax = 0.001 eV/Å, a = 5.4653 Å) |
| `build_index.json` | ZSL match metadata + manuscript convention notes |
| **`wad_Si_alpha_PBE.json`** ★ | **W_ad master (3 entries × stages)** |
| `scripts/22,23,24,25_si_alpha_*.py` | Search + build + PFP relax + tight bulk relax |

## W_ad calculation procedure

Same paper-style recipe as α/β (cell-relax interface, then independent cell-relax of split slabs):

```
W_ad = (E_α-slab + E_Si-slab − E_int) × eV→J / (A_int × Å²→m²)
```

## Termination dependence (Si/α(600))

| Termination | E_int (eV) | E_α-slab (eV) | E_Si-slab (eV) | W_ad (J/m²) |
|---|---:|---:|---:|---:|
| Sn-term | -217.36 | -116.58 (low γ at free surface) | -101.32 | -0.24 |
| **Co-term** | **-217.08** | **-113.38 (high γ at free surface)** | **-101.32** | **+1.06** |

→ The Co-terminated α(600) is the **buried-interface-stable termination** (Si₃N₄/TiN-like reconstruction).
The high-surface-γ termination prefers being at the buried interface (Section S2.2.3 of supplement).
