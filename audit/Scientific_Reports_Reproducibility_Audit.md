# Scientific Reports: bibliography and numerical reproducibility audit
Date: 15 November 2026

## Materials inspected
- Uploaded 26-page manuscript PDF and LaTeX source
- `mechanistic_dose_response_unified_framework.bib`
- `dose-response_update_band(1).py`
- `simulation_summary_with_band.xlsx` (four worksheets)

## Bibliography: structural verification
- 129 BibTeX entries; 60 distinct citation keys used in active manuscript text.
- All 60 cited keys are defined. No missing citation keys or duplicate DOI strings were detected in the BibTeX source.
- The uploaded PDF has 60 numbered bibliography items; thus it is not suffering from unresolved bibliography placeholders.
- 69 uncited entries remain in the database; this is not an error.
- **Limit:** bibliographic author/title/DOI metadata were not individually verified against publisher/Crossref records. This is a structural check, not full bibliographic fact-checking.

## Numerical workbook: full internal consistency check
- Workbook: 12 dose rows, 100 Monte Carlo repetitions, 360,000 host-level observations and 4,800 trajectory time rows.
- Each of 1,200 dose-by-repetition groups contains 300 hosts; no missing worksheet cells.
- All host outcomes exactly equal `int(max_pathogen_burden >= 1e6)`.
- `P0 = eta * dose` agrees within numerical floating-point rounding.
- All 1,200 run-level proportions agree with host-level outcome averages (maximum discrepancy 5.6e-17).
- Summary means, medians, 2.5th and 97.5th percentiles agree with calculations directly from all 100 runs (maximum discrepancy 5.6e-17).
- Twelve displayed example trajectories (4 at each of 3 doses), 400 sampled times each, agree exactly with their corresponding host maxima.
- At D=1e8, the median estimated positive-response proportion is 0.685, with percentile limits 0.630 and 0.726667.

## Independent ODE replay (first Monte Carlo repetition)
- Reproduced 3,600 ODE integrations independently using SciPy RK45, the original seed, dose order, parameter draws, tolerances and time grid.
- All 3,600 integrations reported success.
- Maximum absolute difference in peak pathogen burden from workbook: 1.21e-7; maximum relative difference: 7.0e-13.
- Parameter draws agree to floating-point precision.
- **Limit:** this is an independent replay of 1 of 100 Monte Carlo repetitions, not all 360,000 integrations.

## Remaining pre-submission tasks
1. Archive the final script and full output workbook together with an environment specification, e.g. Python/NumPy/SciPy/pandas versions, and an explicit solver-success check in the released script.
2. Consider continuous event detection or demonstrate that the 400-time-point threshold-crossing approximation is insensitive to time-grid refinement.
3. Keep the numerical endpoint labeled finite-time pathogen-burden threshold crossing, not validated infection or permanent establishment.
4. Keep the shaded band labeled a Monte Carlo percentile interval for finite-host estimates; do not describe it as a biological predictive interval.
5. Verify publication metadata for the 60 cited works against DOI/publisher records before final submission.
6. Replace prospective Data/Code Availability wording with the actual deposited repository DOI or supplementary-file identifiers once those files are deposited.

## Conclusion
Internal numerical calculations are verified, and an independent replay of the first 3,600 trajectories agrees with the workbook. Citation keys resolve structurally and the uploaded PDF includes a numbered bibliography. Full independent replay, time-grid convergence, publisher-level bibliographic metadata verification, and public archival are still outstanding.
