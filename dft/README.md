# dft/ — DFT (OpenMX/PBE) calculation data on SQUID

All DFT calculation directories preserved as-is from the SQUID supercomputer working directory.
Each subfolder corresponds to one SQUID job. Naming prefix indicates status:

- **`Done_*`** — calculation completed successfully (geometry/SCF converged)
- **`Stop_*`** — calculation cancelled or stopped before completion (kept as historical record)

Top-level files:
- **`wad_dft_result.json`** — derived analysis: W_ad,DFT for α(600)/β(100) interface (= 2.412 J/m²)
- **`gamma_DFT.json`** — derived analysis: γ_DFT values for the surfaces evaluated (γ(600) = 0.545, γ(010) = 0.878, γ(301) = 0.685 J/m²)
- **`Phase3_CDD/`** — derived charge-density-difference results (auxiliary, NOT used in the main paper after revision; only the JSON summary, the z-profile PNG, and the Fortran helpers are tracked in git — see "What is in git vs offline" below)

## What is in git vs offline

Each `Done_*` / `Stop_*` directory has been preserved on disk in full
(in.dat, all `test.*` outputs, restart `*_rst/`, electron-density `.cube`
files, the `openmx_vec` MPI binary, and the multi-GB `log.txt`). To keep
this repository at a reasonable size, only the small **textual results**
needed to reproduce the surface energy and work-of-adhesion numbers are
tracked in git:

| File / pattern | In this git repo? | Where to find the rest |
|---|:---:|---|
| `in.dat`, `in.dat#`, `job_SQUID_vec.sh*` | ✓ | — |
| `test.cif` (final relaxed structure)    | ✓ | — |
| `test.ene` (Utot per MD step)           | ✓ | — |
| `test.md`, `test.md2`, `test.xyz`, `test.bulk.xyz` (trajectory) | ✓ | — |
| `test.out` (final summary)              | ✓ | — |
| `test.DFTSCF`, `test.SD`, `test.MC`, `test.EV` (small SCF aux.)  | ✓ | — |
| `log.txt` (full OpenMX stdout, multi-GB)| ✗ | offline / Zenodo |
| `openmx_vec` (~85 MB MPI binary)        | ✗ | offline / Zenodo |
| `test_rst/` (restart files, hundreds of MB) | ✗ | offline / Zenodo |
| `test.*.cube` (electron density / KS potential, ~20 MB each)   | ✗ | offline / Zenodo |
| `Phase3_CDD/*.cube`, `Phase3_CDD/*.xsf` | ✗ | offline / Zenodo |

The numerical values that the paper actually quotes are recovered from the
small text files alone — `test.ene`, `test.md`, and `test.out` are
sufficient to reconstruct E_total, the relaxed cell, the in-plane area,
and the converged Utot. The cube/restart artefacts are only required for
re-running OpenMX or for re-doing the CDD post-processing, and will be
included in the future Zenodo deposit.

## Directory inventory

### Done_* (used in paper / converged)

| Folder | SQUID job | Atoms | Paper element | Notes |
|---|---|---:|---|---|
| `Done_Bulk_CoSn3/` | Bulk_CoSn3 | 32 | **Table S1**, μ_α reference | α-CoSn₃ bulk full-cell relax (OptC5). Utot = −2562.0571508 Ha → μ_α = −80.0643 Ha/atom |
| `Done_Surf_CoSn3(600)/` | Surf_CoSn3(600) | 32 | **Table S2** (γ_DFT = 0.545 J/m²) | (600) low-γ termination slab, OptC5 |
| `Done_Surf_CoSn3(010)/` | Surf_CoSn3(010) | 96 | **Table S2** (γ_DFT = 0.878 J/m²) | (010) low-γ termination slab, OptC5 (610 MD steps) |
| `Done_CoSn3(600)+Sn(100)_Phase1A/` | Phase1A | 56 | **Fig. 9, Table S4** (W_ad,DFT = 2.41 J/m²) | AB interface OptC5 (cell + atoms relaxed) |
| `Done_CoSn3(600)+Sn(100)_Phase1B/` | Phase1B | 32 | Same | α slab from interface, OptC5 |
| `Done_CoSn3(600)+Sn(100)_Phase1C/` | Phase1C | 24 | Same | β slab from interface, OptC5 |
| `Done_CoSn3(600)+Sn(100)_Phase2A/` | Phase2A | 56 | **Auxiliary** (CDD prep, not in main paper) | AB single-point SCF on relaxed-interface geometry; produces `AB.cube` (offline only) |
| `Done_CoSn3(600)+Sn(100)_Phase2B/` | Phase2B | 32 | **Auxiliary** | A (α only) single-point SCF in interface supercell; produces `A.cube` (offline only) |
| `Done_CoSn3(600)+Sn(100)_Phase2C/` | Phase2C | 24 | **Auxiliary** | B (β only) single-point SCF in interface supercell; produces `B.cube` (offline only) |
| `Phase3_CDD/` | (post-processing) | — | **Auxiliary** | dAB.cube = ρ_AB − (ρ_A + ρ_B) computed from Phase2 cubes (offline only); see Phase3 section below |

### Stop_* (cancelled, not used in paper)

These were initially submitted to SQUID but cancelled when the project pivoted to PFP-primary methodology.
Kept as historical record / partial data.

| Folder | SQUID job | Status when cancelled | Reason cancelled |
|---|---|---|---|
| `Stop_CoSn3(312)+Sn(100)_Phase1A/` | Phase1A | Submitted only (not started) | PFP代替済 (PFP value 0.82 J/m² adequate) |
| `Stop_CoSn3(312)+Sn(100)_Phase1B/` | Phase1B | Mid-run, partial relaxation data | Same |
| `Stop_CoSn3(312)+Sn(100)_Phase1C/` | Phase1C | Submitted only | Same |
| `Stop_CoSn3(600)+Sn(001)_Phase1A/` | Phase1A | Mid-run, stuck in steepest-descent fallback | PFP代替済 (PFP value 1.89 J/m² adequate) |
| `Stop_CoSn3(600)+Sn(001)_Phase1B/` | Phase1B | Mid-run | Same |
| `Stop_CoSn3(600)+Sn(001)_Phase1C/` | Phase1C | Submitted only | Same |

## Computational settings (all DFT runs, identical except where noted)

| Item | Value |
|---|---|
| Code | OpenMX 3.9.9 (PBE19 norm-conserving database) |
| XC functional | GGA-PBE |
| Electronic temperature | 300 K |
| Spin polarization | On (`scf.SpinPolarization` in `in.dat`; Co is magnetic) |
| Pseudo-atomic-orbital basis | Co_PBE19S (Co: Co6.0S-s2p3d2f1), Sn_PBE19 (Sn: Sn7.0-s2p2d3f1) |
| Initial spin seed (Co UP/DOWN) | 8.5 / 6.5 (μ ≈ 2) |
| Real-space grid cutoff | 200 Ryd (per paper §S1.2) |
| SCF criterion | 1.0 × 10⁻⁷ Hartree |
| Force criterion | 1.0 × 10⁻³ Hartree/Bohr (≈ 0.05 eV/Å) |
| Geometry optimisation | simultaneous cell + atomic relaxation (`MD.Type OptC5`) |
| k-grid | target spacing Δk ≈ 0.15 rad/Å |
| Walltime / parallelisation | 120 h (system max); 64 VEs × 5 cores, OMP 2 → 320 MPI |

## File contents (per directory)

A typical OpenMX run directory holds the files below. The "git" column
shows whether the file is included in this public repository.

| File | git | Content |
|---|:---:|---|
| `in.dat` | ✓ | OpenMX input file (final version used) |
| `in.dat#` | ✓ | Auto-generated backup of in.dat |
| `job_SQUID_vec.sh` | ✓ | PBS submit script for SQUID Vector engine |
| `job_SQUID_vec.sh.o######` | ✓ | Standard output from SQUID job |
| `job_SQUID_vec.sh.e######` | ✓ | Standard error from SQUID job |
| `openmx_vec` | ✗ | OpenMX MPI binary (~85 MB) — kept offline |
| `log.txt` | ✗ | Full OpenMX stdout (multi-GB) — kept offline |
| **`test.cif`** | ✓ | **Final relaxed structure (CIF)** |
| **`test.ene`** | ✓ | **Energy trajectory per MD step** (Utot in column 14) |
| **`test.md`** | ✓ | **Geometry trajectory per MD step** (cell + atom positions) |
| `test.md2` | ✓ | Alternate trajectory format |
| `test.out` | ✓ | Summary output (numeric force/energy/cell at end) |
| `test.xyz`, `test.bulk.xyz` | ✓ | XYZ format snapshots |
| `test.DFTSCF`, `test.SD`, `test.MC`, `test.EV` | ✓ | Small SCF / MD auxiliary files |
| `test.tden.cube` | ✗ | Total electron density (ρ↑ + ρ↓), Gaussian cube — offline |
| `test.dden.cube` | ✗ | Spin density (ρ↑ − ρ↓) — offline |
| `test.den0.cube`, `test.den1.cube` | ✗ | Spin-up / spin-down densities separately — offline |
| `test.sden.cube` | ✗ | Spin density (alternate format) — offline |
| `test.v0.cube`, `test.v1.cube` | ✗ | Kohn–Sham potential (spin-up / spin-down) — offline |
| `test.vhart.cube` | ✗ | Hartree potential — offline |
| `test_rst/` | ✗ | Restart files (OpenMX wave-function snapshots) — offline |
| `AB.cube`, `A.cube`, `B.cube` | ✗ | (Phase2 only) named CDD-input cubes copied from `test.tden.cube` — offline |

## Phase3 CDD files

`Phase3_CDD/` contains:

| File | git | Content |
|---|:---:|---|
| `A.cube`, `B.cube`, `AB.cube` | ✗ | Copied from Phase2A/B/C (offline) |
| `A_B.cube` | ✗ | Sum (A + B) computed by `add_gcube` (offline) |
| `dAB.cube` | ✗ | Charge density difference, AB − (A + B), by `diff_gcube` (offline) |
| `dAB.xsf` | ✗ | Same Δρ in XCrySDen format for VESTA (offline) |
| `dAB_z_profile.png` | ✓ | Planar-averaged Δρ(z) profile |
| `cdd_summary.json` | ✓ | Numerical summary (planar-averaged Δρ peaks, integrated charge transfer) |
| `add_gcube`, `cube2xsf`, `diff_gcube` | ✓ | Fortran helper executables |

**Important**: Fig. 9(b) and (c) showing CDD were removed from the paper during revision
(the metallic α-CoSn₃/β-Sn interface produces only modest charge redistribution; W_ad alone
is sufficient to make the bonding argument). Phase3 data retained here as an auxiliary
record. The cube/xsf files needed to regenerate the figures will be in the Zenodo deposit.

## How W_ad,DFT is computed

```
W_ad,DFT = (E_α-slab + E_β-slab − E_AB) × eV→J / (A_int × Å²→m²)

E_α-slab  ←  Done_CoSn3(600)+Sn(100)_Phase1B/test.ene  (last column 14, last row)
E_β-slab  ←  Done_CoSn3(600)+Sn(100)_Phase1C/test.ene
E_AB      ←  Done_CoSn3(600)+Sn(100)_Phase1A/test.ene
A_int     ←  Done_CoSn3(600)+Sn(100)_Phase1A/test.md (last cell)
```

See `wad_dft_result.json` for the final numerical result and intermediate values.

## How γ_DFT is computed

```
γ_DFT(face) = (E_slab − N · μ_α) × eV→J / (2 · A · Å²→m²)

E_slab  ←  Done_Surf_CoSn3(600)/test.ene
N       =  number of atoms in slab (32 for (600))
μ_α     =  −80.0643 Ha/atom (= Done_Bulk_CoSn3/test.ene last Utot / 32)
A       =  in-plane area from Done_Surf_CoSn3(600)/test.md (last cell vectors)
```
