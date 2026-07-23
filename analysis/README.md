# DStressR Analysis Area

This folder is intentionally outside the CRAN package. `.Rbuildignore`
excludes `analysis/`, so the package submitted to CRAN stays lightweight and
self-contained while the repository can still keep reproducibility and
manuscript workflows nearby.

## Layout

- `ecoli_promoter_screen/`: public *E. coli* promoter-compound screen
  reproducibility checks and manuscript application analyses. These scripts use
  the dataset shipped with DStressR plus public author score/Z-score data.
- `campylobacter_manuscript/`: the previous Campylobacter comparison and figure
  workflow. These scripts depend on local package-generated TSV outputs and
  external data that are not shipped with the package.
- `outputs/`: generated local tables and figures. This directory is ignored by
  Git and is not part of the package.

Analysis scripts may consume package outputs, but should not reimplement
estimators, p-value calculations, empirical-Bayes moderation, replicate
aggregation, or multiple-testing correction. If a statistic belongs to the
package API, add it under `R/` and test it under `tests/`.
