# Scientific Reports: time-grid and reference metadata audit update

## Numerical time-grid sensitivity
The original model was independently integrated with SciPy RK45 (rtol=1e-6, atol=1e-9) using dense output. For each of 12 dose levels, the first 10 virtual hosts from the first Monte Carlo repetition were tested (120 trajectories total). The maximum pathogen burden evaluated on 400 equally spaced time points was compared against 1600 points over t=0..30.

- Successful integrations: 120/120.
- Threshold-crossing classifications changed: 0/120.
- Disagreements with saved classifications: 0/120.
- Maximum relative change in sampled peak pathogen burden: 2.771745e-5 (0.00277%).

**Scope:** This is a stratified 120-host numerical sensitivity check, not a complete 360,000-trajectory convergence study. The full-run attempt exceeded execution time and was not completed. The test supports the existing sampling grid for the checked subset, but does not rule out rare between-grid crossings elsewhere.

## Bibliographic metadata checked against publisher/PubMed records
Verified and corrected in the accompanying BibTeX file:

1. Abe et al. 2020, DOI 10.1128/AEM.01551-20: original title was inaccurate. Correct title: *Modeling Invasion of Campylobacter jejuni into Human Small Intestinal Epithelial-Like Cells by Bayesian Inference*. https://pubmed.ncbi.nlm.nih.gov/33067190/
2. Abe et al. 2021, DOI 10.1128/AEM.01299-21: original title was inaccurate. Correct title: *A New Dose-Response Model for Estimating the Infection Probability of Campylobacter jejuni Based on the Key Events Dose-Response Framework*. https://pubmed.ncbi.nlm.nih.gov/34347512/
3. Sego et al. 2020, DOI 10.1371/journal.pcbi.1008451: original title was abbreviated/inaccurate. Correct full title: *A modular framework for multiscale, multicellular, spatiotemporal modeling of acute primary viral infection and immune response in epithelial tissues and its application to drug therapy timing and effectiveness*. https://pmc.ncbi.nlm.nih.gov/articles/PMC7785254/
4. The manuscript's `Schmidt2013QMRA` was cited as a 2013 Wiley-Blackwell book by Schmidt, Pintar and Fazil. This book citation could not be substantiated. The key has been replaced with a verified, relevant 2013 Risk Analysis article by Schmidt, Pintar, Fazil and Topp, DOI 10.1111/risa.12006. https://pubmed.ncbi.nlm.nih.gov/23311599/

The corrected bibliography preserves citation keys, so manuscript LaTeX citations remain syntactically intact. **Not all 60 cited references have undergone publisher-level verification.** Other records still require review before claiming full bibliographic verification.

## Recommendations
- Replace the manuscript's original .bib file with the corrected file (rename to the expected basename).
- Recompile LaTeX bibliography to propagate corrected titles and the replacement Schmidt reference.
- Add the sensitivity test description to numerical methods or supplementary methods if desired.
- Check the remaining bibliographic records against publisher metadata and archive the final code/data with a persistent location.
