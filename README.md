# Mechanistic dose–response modeling in infectious diseases

Research materials accompanying the manuscript **Mechanistic dose–response modeling in infectious diseases: a unified framework across biological scales** (Scientific Reports submission preparation).

## Repository layout

- `manuscript/Scientific_Reports_Manuscript.tex` — manuscript source.
- `manuscript/mechanistic_dose_response_unified_framework.bib` — corrected bibliography (same basename expected by LaTeX).
- `manuscript/*.png` — two referenced numerical figures.
- `manuscript/Scientific_Reports_preview.pdf` — **previously compiled preview**, not rebuilt with the corrected bibliography. Do not treat its references as final.
- `scripts/dose_response_with_band.py` — principal simulation script for the interval figure and 360,000-host workbook.
- `scripts/dose_response_original.py` and `scripts/dose_response_updated.py` — earlier simulation versions, retained for provenance.
- `scripts/verify_replay.py` — independent replay of first 3,600 host integrations.
- `scripts/time_grid_sensitivity.py` — 400 vs 1,600 time-point check on 120 sampled hosts.
- `data/` — saved Excel outputs; the `simulation_summary_with_band.xlsx` workbook is the main result.
- `results/` — numerical sensitivity summary.
- `audit/` — scope, findings and limitations of completed checks.

## Python environment

Use Python 3.10+ (recommended). From the repository root:

```bash
python -m pip install -r requirements.txt
python scripts/verify_replay.py
python scripts/time_grid_sensitivity.py
```

The verification scripts use repository-relative paths. To rerun the **main simulation**, change into `data/` first so its relative output files are written there, e.g. `cd data && python ../scripts/dose_response_with_band.py`. This can take considerable time and overwrite saved outputs. Preserve the original workbook before rerunning.

## Compile LaTeX

The manuscript references the two PNGs and the bibliography by basename, all stored in `manuscript/`. Compile from that directory using a TeX installation with BibTeX:

```bash
pdflatex Scientific_Reports_Manuscript.tex
bibtex Scientific_Reports_Manuscript
pdflatex Scientific_Reports_Manuscript.tex
pdflatex Scientific_Reports_Manuscript.tex
```

Alternatively upload the contents of `manuscript/` to Overleaf.

## Verification status and caveats

- Saved output contains 12 doses × 100 Monte Carlo repetitions × 300 hosts = 360,000 host records.
- Saved host-level outcomes reproduce the workbook summaries.
- First repetition (3,600 trajectories) was independently replayed; see audit.
- Time-grid sensitivity compared 400 vs 1,600 evaluation points for **120 trajectories**, with no changed classifications in that subset. It is not a full 360,000-trajectory convergence study.
- Four bibliographic entries were corrected or replaced; **the remaining references were not all independently verified against publisher metadata**.
- The preview PDF predates those bibliography corrections; regenerate it for submission.
- The model's positive outcome is finite-time pathogen-burden threshold crossing, not clinically validated infection.
- Before publishing this repository, review author contact details, third-party figure permissions, and any journal policies. Add a LICENSE only after deciding on reuse permissions.
