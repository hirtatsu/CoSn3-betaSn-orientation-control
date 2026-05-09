# pfp/bulks/ — PFP/PBE relaxed bulks

Reference structures for chemical potentials μ used in γ and W_ad calculations.

## Files

| Phase | CIF | JSON |
|---|---|---|
| α-CoSn₃ (Co₈Sn₂₄, 32-atom conventional cell) | `bulk_alpha_CoSn3_PBE.cif` | `bulk_alpha_CoSn3_PBE.json` |
| β-Sn (4-atom tetragonal) | `bulk_beta_Sn_PBE.cif` | `bulk_beta_Sn_PBE.json` |
| Si (8-atom Fd-3m) | `Si_PBE.cif` | `Si_PBE.json` |

## Lattice constants (PFP/PBE)

| Phase | a (Å) | b (Å) | c (Å) | μ (eV/atom) |
|---|---:|---:|---:|---:|
| α-CoSn₃ | 17.4355 | 6.2594 | 6.2504 | −3.71287 |
| β-Sn | 5.9292 | 5.9292 | 3.2008 | −3.11613 |
| Si | 5.4653 | 5.4653 | 5.4653 | −4.55263 |

## Settings

- ExpCellFilter (full cell + atoms)
- BFGS optimizer
- fmax = 0.001 eV/Å (tight)
- PFP v8 PBE mode

## Summary
`bulk_summary_PBE.json` collects α-CoSn₃ and β-Sn lattice/μ in one file for convenience (Si is in its own JSON).
