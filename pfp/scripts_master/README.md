# scripts/ — production pipeline

14 scripts numbered in execution order. To reproduce results from scratch, run in order. Some scripts run locally (need `pymatgen` from `/home/tatsumi/test-env/bin/python`); others run on `matlantis` (need PFP).

## Pipeline overview

```
data/*.cif                                                          (input CIFs)
        │
        ▼
   01_load_and_verify.py        [LOCAL] sanity check, reproduce manuscript Table 3 disregistry
        │
        ▼
   05_bulk_relax_matlantis.py   [MATLANTIS] PFP/PBE relax of α-CoSn₃ and β-Sn bulks
        │
        ▼
   06_interface_search_relaxed.py  [LOCAL] pmg ZSL search on PBE-relaxed bulks
        │
        ▼
   07_stack_chosen_relaxed.py   [LOCAL] build 10 initial interface CIFs (chosen ranks per pair)
        │
        ▼
   08b_split_to_slabs.py        [LOCAL] split each interface CIF → α-only and β-only slab CIFs
        │
        ▼ (upload to matlantis)
        ▼
   09_slabs_and_interfaces_pfp.py  [MATLANTIS] initial 30 atom-only PFP/PBE relaxations
        │
        ▼
   11_phaseA_and_PBED3_bulk.py  [MATLANTIS] re-relax 4 unconverged interfaces + PBE+D3 bulk
        │
        ▼
   12_pbed3_stack_split.py      [LOCAL] PBE+D3 stack build for 3 sensitivity-check pairs
        │
        ▼ (upload)
        ▼
   13_pbed3_relax.py            [MATLANTIS] 9 PBE+D3 relaxations (atom-only)
        │
        ▼
   16b_local_alpha_slab_gen.py  [LOCAL] α slab CIFs with full pmg termination enumeration
        │                       (3 faces × 3 terms = 9 slab CIFs for cell-relax)
        ▼ (upload)
        ▼
   17_alpha_surface_relax.py    [MATLANTIS] α slab atom-only relax (initial pass)
        │
        ▼
   19_paper_style_wad.py        [MATLANTIS] β natural-cell slabs (rescale, atom-only)
        │
        ▼
   20_interface_cell_relax.py   [MATLANTIS] cell-relax all 10 interfaces (FrechetCellFilter)
        │
        ▼
   21_slabs_cell_relax.py       [MATLANTIS] ★ MASTER — cell-relax all standalone slabs +
                                            recompute paper-style W_ad
                                            → results/wad_paper_full_consistent.json
```

## Script-by-script summary

| # | Script | Where | Purpose |
|---:|---|---|---|
| 01 | `01_load_and_verify.py` | local | Load CIFs, verify cell, reproduce manuscript Table 3 |
| 05 | `05_bulk_relax_matlantis.py` | matlantis | PFP/PBE bulk relax of both materials |
| 06 | `06_interface_search_relaxed.py` | local | ZSL search with PBE-relaxed bulks (top-5 candidates) |
| 07 | `07_stack_chosen_relaxed.py` | local | Build 10 interface CIFs at chosen ranks |
| 08b | `08b_split_to_slabs.py` | local | Split interface CIFs into α/β slab CIFs |
| 09 | `09_slabs_and_interfaces_pfp.py` | matlantis | 30 atom-only PFP/PBE relaxations |
| 11 | `11_phaseA_and_PBED3_bulk.py` | matlantis | Re-relax 4 unconverged + PBE+D3 bulks |
| 12 | `12_pbed3_stack_split.py` | local | PBE+D3 stack for 3 pairs (sensitivity check) |
| 13 | `13_pbed3_relax.py` | matlantis | 9 PBE+D3 relaxations |
| 16b | `16b_local_alpha_slab_gen.py` | local | α slabs (full term enum) for γ calculation |
| 17 | `17_alpha_surface_relax.py` | matlantis | α slab atom-only relax (input for 21) |
| 19 | `19_paper_style_wad.py` | matlantis | β natural-cell slab init (atom-only, input for 21) |
| 20 | `20_interface_cell_relax.py` | matlantis | Cell-relax 10 interfaces (FrechetCellFilter) |
| **21** | **`21_slabs_cell_relax.py`** | **matlantis** | **★ MASTER — final paper-style W_ad** |

## Where the master output is written

`results/wad_paper_full_consistent.json` contains the final W_ad table, paper-style with full cell relaxation. See `results/README.md` for output index.

## `_archive/`

9 superseded or diagnostic scripts. See `_archive/README.md` for details. Kept for pedagogical value (gotchas, methodological journey).
