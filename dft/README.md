# dft/ — DFT (OpenMX/PBE) calculation data on SQUID

All DFT calculation directories preserved as-is from the SQUID supercomputer working directory.
Each subfolder corresponds to one SQUID job. Naming prefix indicates status:

- **`Done_*`** — calculation completed successfully (geometry/SCF converged)
- **`Stop_*`** — calculation cancelled or stopped before completion (kept as historical record)

Top-level files:
- **`wad_dft_result.json`** — derived analysis: W_ad,DFT for α(600)/β(100) interface (= 2.412 J/m²)
- **`gamma_DFT.json`** — derived analysis: γ_DFT values for available surfaces (γ(600) = 0.545, γ(010) = 0.878 J/m²; γ(301) provisional pending convergence)
- **`Phase3_CDD/`** — derived charge density difference cubes (auxiliary, NOT in main paper after revision)

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
| `Done_CoSn3(600)+Sn(100)_Phase2A/` | Phase2A | 56 | **Auxiliary** (CDD prep, not in main paper) | AB single-point SCF on relaxed-interface geometry → AB.cube |
| `Done_CoSn3(600)+Sn(100)_Phase2B/` | Phase2B | 32 | **Auxiliary** | A (α only) single-point SCF in interface supercell → A.cube |
| `Done_CoSn3(600)+Sn(100)_Phase2C/` | Phase2C | 24 | **Auxiliary** | B (β only) single-point SCF in interface supercell → B.cube |
| `Phase3_CDD/` | (post-processing) | — | **Auxiliary** | dAB.cube = ρ_AB − (ρ_A + ρ_B) computed from Phase2 cubes |

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
| Code | OpenMX 3.9.9 |
| XC functional | GGA-PBE |
| Spin polarization | On |
| Pseudopotentials | Co_PBE19S (15 valence, Soft), Sn_PBE19 (14 valence) |
| PAO basis | Co6.0S-s2p3d2f1, Sn7.0-s2p2d3f1 |
| Initial spin (Co UP/DOWN) | 8.5 / 6.5 (μ = 2 seed) |
| Real-space grid cutoff | 200 Hartree |
| SCF criterion | 1.0 × 10⁻⁷ Hartree |
| Force criterion | 1.0 × 10⁻³ Hartree/Bohr (≈ 0.05 eV/Å) |
| MD type | OptC5 (cell + atoms simultaneous) |
| k-grid | target Δk ≈ 0.15 rad/Å |
| Walltime | 120 h (system max), 64 VEs × 5 cores / OMP 2 = 320 MPI |

## File contents (per directory)

Typical OpenMX run directory contains:

| File | Content |
|---|---|
| `in.dat` | OpenMX input file (final version used) |
| `in.dat#` | Auto-generated backup of in.dat |
| `job_SQUID_vec.sh` | PBS submit script for SQUID Vector engine |
| `job_SQUID_vec.sh.o######` | Standard output from SQUID job |
| `job_SQUID_vec.sh.e######` | Standard error from SQUID job |
| `openmx_vec` | OpenMX executable (compiled for SQUID Vector engine, ~80 MB) |
| **`log.txt`** | **Full OpenMX log (large; latest force/energy at end)** |
| **`test.cif`** | **Final relaxed structure (CIF)** |
| **`test.ene`** | **Energy trajectory per MD step** (Utot in column 14) |
| **`test.md`** | **Geometry trajectory per MD step** (cell + atom positions) |
| `test.md2` | Alternate trajectory format |
| `test.out` | Summary output (numeric force/energy/cell at end) |
| `test.xyz`, `test.bulk.xyz` | XYZ format snapshots |
| `test.tden.cube` | Total electron density (ρ↑ + ρ↓), Gaussian cube format |
| `test.dden.cube` | Spin density (ρ↑ − ρ↓) |
| `test.den0.cube`, `test.den1.cube` | Spin-up and spin-down densities separately |
| `test.sden.cube` | Spin density (alternate format) |
| `test.v0.cube`, `test.v1.cube` | Kohn-Sham potential (spin-up, spin-down) |
| `test.vhart.cube` | Hartree potential |
| `test_rst/` | Restart files (OpenMX wave function snapshots, used for resume runs) |
| `AB.cube`, `A.cube`, `B.cube` | (Phase2 only) named CDD-input cubes copied from test.tden.cube |

## Phase3 CDD files

`Phase3_CDD/` contains:
- `A.cube`, `B.cube`, `AB.cube` (copied from Phase2A/B/C)
- `A_B.cube` = sum (A + B) computed by `add_gcube`
- **`dAB.cube`** = AB − (A + B), the charge density difference, by `diff_gcube`
- `dAB.xsf` = same Δρ in XCrySDen format (for VESTA)
- `dAB_z_profile.png` = planar-averaged Δρ(z) profile
- `cdd_summary.json` = numerical summary
- `add_gcube`, `cube2xsf`, `diff_gcube` = helper executables (Fortran)

**Important**: Fig. 9(b) and (c) showing CDD were removed from the paper during revision
(the metallic α-CoSn₃/β-Sn interface produces only modest charge redistribution; W_ad alone
sufficient to make the bonding argument). Phase3 data retained as auxiliary record.

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
