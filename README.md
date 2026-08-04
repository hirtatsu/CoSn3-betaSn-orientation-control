# Computational data — Crystallographic orientation control of β-Sn via α-CoSn₃ thin films

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

Companion data and code repository for:

> X. Wang, H. Tatsumi, C.-L. Li, F.-C. Yang, Z. He, I-E. Chen, C. R. Kao,
> L.-C. Chang, J.-W. Lee, H. Nishikawa.
> *Crystallographic orientation control of β-Sn via preferentially oriented
> α-CoSn₃ thin films*. **Manuscript in preparation, 2026**.

This repository contains every input, output, and analysis script required
to reproduce the surface-energy and work-of-adhesion calculations behind the
texture-control argument:

- **PFP/PBE** (Preferred Potential v8 on the Matlantis platform) — primary
  evaluation of γ for four α-CoSn₃ low-index faces and W_ad for three
  α/β interface pairs and two Si/α interface pairs.
- **DFT/PBE** (OpenMX 3.9.9 with the PBE19 norm-conserving pseudopotential
  database, on SQUID @ The University of Osaka) — cross-validation of
  α(600)/β(100) work of adhesion and α-CoSn₃(600),(010),(301) surface
  energies; auxiliary charge-density-difference (CDD) calculation.

## Citation

```
X. Wang, H. Tatsumi, C.-L. Li, F.-C. Yang, Z. He, I-E. Chen, C. R. Kao,
L.-C. Chang, J.-W. Lee, H. Nishikawa.
"Crystallographic orientation control of β-Sn via preferentially oriented
α-CoSn₃ thin films". (2026, manuscript in preparation).
```

A Zenodo DOI will be added on first stable release.

## Directory layout

```
cosn3-betasn-paper-data/
├── README.md                    ← this file (overall navigator)
├── MANUSCRIPT_DATA_INDEX.md     ← which file produces which Table/Figure
├── LICENSE                      MIT (code) + CC-BY-4.0 (data)
├── CITATION.cff
├── .gitignore
│
├── pfp/                         ★ MLIP / Matlantis data
│   ├── README.md
│   ├── bulks/                   α-CoSn₃, β-Sn, Si bulk relax (Table S1)
│   ├── alpha_surface_energies/  4 face × 3 termination γ (Table 1, S2)
│   ├── alpha_beta_Wad/          3 α/β interfaces W_ad (Fig. 9, Table S5)
│   ├── Si_alpha_Wad/            Si/α(600), Si/α(312) W_ad (Fig. 7, Table S3)
│   └── scripts_master/          numbered reproducibility scripts
│
└── dft/                         ★ OpenMX / SQUID data (workdirs preserved)
    ├── README.md
    ├── wad_dft_result.json      W_ad,DFT = 2.412 J/m² (derived)
    ├── gamma_DFT.json           γ_DFT for available surfaces (derived)
    ├── Phase3_CDD/              auxiliary CDD analysis (excluded from git)
    ├── Done_Bulk_CoSn3/         μ_α reference (Table S1, S2)
    ├── Done_Surf_CoSn3(600)/    γ(600)_DFT = 0.545 J/m²  (Table S2)
    ├── Done_Surf_CoSn3(010)/    γ(010)_DFT = 0.878 J/m²  (Table S2)
    ├── Done_Surf_CoSn3(301)/    γ(301)_DFT = 0.787 J/m²  (Table S2)
    ├── Done_CoSn3(600)+Sn(100)_Phase{1A,1B,1C}/  α(600)/β(100) W_ad (Fig. 9)
    ├── Done_CoSn3(600)+Sn(100)_Phase{2A,2B,2C}/  AB/A/B SCF for CDD (auxiliary)
    └── Stop_*/                  cancelled DFT jobs (PFP-replaced; archival only)
```

## Headline results (from the paper)

### α-CoSn₃ surface energies (PFP/PBE; Table 1)

| Plane | γ (J/m²) | Literature DFT (Ma; Wang *et al.*) |
|---|---:|---:|
| (600) | **0.46** | 0.46 |
| (301) | 0.66    | 0.78 |
| (321) | 0.73    | 1.09 |
| (010) | 0.77    | 1.05 |

DFT (this work) cross-check: γ(600) = 0.545 J/m², γ(010) = 0.878 J/m²,
γ(301) = 0.787 J/m².

### Si / α-CoSn₃ adhesion + lattice match (Fig. 6, 7; Table S3, S6)

| Pair | W_ad (J/m²) | Disregistry (in-plane) |
|---|---:|---|
| Si(100) / α-CoSn₃(312) | **2.57** | 0.99 % along Si⟨001⟩ ‖ α⟨11̄2⟩ |
| Si(100) / α-CoSn₃(600), Co-term | 1.06 | 13.38 % along Si⟨001⟩ ‖ α⟨001⟩ |

These numbers explain the experimentally observed (600) → (312) texture
switch as a transition from surface-energy minimization (DC-sputter, high
adatom mobility) to substrate-interface-energy minimization (RF-sputter,
low mobility).

### Works of adhesion at α-CoSn₃ / β-Sn interfaces (PFP/PBE; Fig. 9, Table S4, Table S5)

| Interface (α / β-Sn) | W_ad (J/m²) | DFT |
|---|---:|---:|
| α-CoSn₃(600) / β-Sn(100), c // substrate | **2.28** | 2.41 (≤ 6 % deviation; Table S4) |
| α-CoSn₃(600) / β-Sn(001), c ⊥ substrate  | 1.89 | — |
| α-CoSn₃(312) / β-Sn(100)                 | 0.82 | — |

The high W_ad of α-CoSn₃(600) / β-Sn(100) explains the strong c-axis-parallel
β-Sn templating observed by EBSD on (600)-textured films, and the low W_ad
of α-CoSn₃(312) / β-Sn(100) explains why (312)-textured films fail to
template β-Sn.

## Computational settings (summary)

### MLIP (PFP/PBE on Matlantis)

- Code: PFP v8, `calc_mode = PBE`, ASE Python interface
- Bulk relax: ExpCellFilter (cell + atoms), BFGS, fmax = 0.001 eV/Å
- Slab/interface relax: FrechetCellFilter mask = [T,T,F,F,F,T] (in-plane
  cell + atoms; vacuum dim fixed), fmax = 0.015 eV/Å
- Vacuum ≈ 15 Å normal to slab; symmetric terminations
- Recipe: paper-style W_ad following H. Tatsumi *et al.*, *Acta Mater.*
  **304**, 121813 (2026), Eq. 4 — α slab and β slab are each cell-relaxed
  independently to their natural equilibria, the interface stack is cell-
  relaxed in-plane, and W_ad = (E_α + E_β − E_int) / A_int.

### DFT (OpenMX 3.9.9 / PBE on SQUID @ The University of Osaka)

- GGA-PBE; norm-conserving pseudopotentials from the OpenMX PBE19 database
- Pseudo-atomic-orbital basis: `Co_PBE19S` (Co: Co6.0S-s2p3d2f1) and
  `Sn_PBE19` (Sn: Sn7.0-s2p2d3f1)
- Real-space grid cutoff: 200 Ryd; SCF convergence: 1.0 × 10⁻⁷ Hartree
- Force criterion 1.0 × 10⁻³ Hartree/Bohr (≈ 0.05 eV/Å)
- Geometry optimisation: simultaneous cell + atomic relaxation
  (`MD.Type OptC5` in `in.dat`)
- k-point grid: target spacing Δk ≈ 0.15 rad/Å
- Electronic temperature: 300 K

### Bulk lattice constants (Table S1)

| Phase | PFP a, b, c (Å) | DFT a, b, c (Å) | Exp a, b, c (Å) |
|---|---|---|---|
| α-CoSn₃ | 17.436, 6.259, 6.250 | 17.207, 6.338, 6.349 | 16.864, 6.268, 6.270 |
| β-Sn    | 5.929, 5.929, 3.201  | 5.930, 5.930, 3.201  | 5.831, 5.831, 3.182 |
| Si      | 5.465, 5.465, 5.465  | —                    | 5.431, 5.431, 5.431 |

## How to navigate

1. Open `MANUSCRIPT_DATA_INDEX.md` — maps every paper element (Fig./Table/Eq.)
   to specific files.
2. Each subfolder has its own `README.md` with details on inputs/outputs/scripts.
3. Master numerical values are in JSON files; structures are in CIF.

## Reproducing the work

> **Note on paths.** Some scripts under `pfp/scripts_master/` and
> `pfp/*/scripts/` retain hardcoded absolute paths from the Matlantis Jupyter
> environment (`/home/jovyan/work_dir/...`) and the analysis workstation
> (`/home/tatsumi/projects/...`). They are kept verbatim as a record of how
> the data were produced. To re-run on your own machine, edit the `WORK = …`
> / `PFP_DIR = …` lines at the top of each script. The numerical results
> they produced are already cached in JSON/CIF form under `pfp/*/`.

### MLIP (Matlantis)
- `pfp/scripts_master/05_bulk_relax_matlantis.py` — bulk relax (Table S1)
- `pfp/scripts_master/16b_local_alpha_slab_gen.py` → `17_alpha_surface_relax.py`
  — α-CoSn₃ slabs and γ (Table 1, S2)
- `pfp/scripts_master/19_paper_style_wad.py`, `20_interface_cell_relax.py`,
  `21_slabs_cell_relax.py` — α/β W_ad (Fig. 9, Table S5)
- `pfp/scripts_master/22…25_si_alpha_*.py` — Si/α W_ad (Fig. 7, Table S3)
- Requires Matlantis (Preferred Networks) account.

### DFT (OpenMX 3.9.9 on SQUID)
- Per-job inputs: `dft/<run>/in.dat` (and `dft/<run>/job_SQUID_vec.sh` job
  scripts, modify queue parameters for your site)
- Pseudopotential database: OpenMX 3.9.9 `DFT_DATA19` (PBE19; not redistributed)
- Reference structures: `pfp/bulks/{alpha_CoSn3,beta_Sn,Si}_PBE.cif`

## Excluded heavy artefacts (preserved offline for future Zenodo archive)

To keep this repository under GitHub's recommended size, the following are
excluded by `.gitignore` and preserved offline for the future Zenodo deposit:

- `**/*_rst/` — OpenMX restart-state directories (~2.7 GB total)
- `**/*.cube` — charge-density / wavefunction cubes for Phase 2A/B/C and
  Phase3_CDD (~1.4 GB total)
- `**/openmx_vec` — OpenMX MPI executable that ended up in each run dir
- `**/log.txt`, `*.xsf` — very large stdout / charge-difference dumps
- `manuscript/` — manuscript drafts (will appear via the journal once published)

The numerical results extracted from these (γ, W_ad, μ_α) are kept as
`*.json` and per-run summary `.dat`/`.out` files.

## Related repositories

- [`hirtatsu/beta-Sn-DFT-PFP-MEAM`](https://github.com/hirtatsu/beta-Sn-DFT-PFP-MEAM)
  — companion DFT/PFP/MEAM benchmark for β-Sn elastic constants and surface
  energies (the protocol underlying the W_ad recipe used here).
- [`hirtatsu/beta-Sn-foundation-MLIP`](https://github.com/hirtatsu/beta-Sn-foundation-MLIP)
  — follow-up benchmark of three universal **foundation MLIPs** (MACE-MPA-0,
  ORB v3, SevenNet-Omni) on β-Sn using the same protocol.

## License

Code: MIT.
Data (CIFs, JSON results, raw OpenMX outputs): CC-BY-4.0.

## Contact

For questions: tatsumi.jwri@osaka-u.ac.jp
