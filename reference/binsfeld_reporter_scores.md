# Binsfeld et al. reporter scores and Z-scores

Long-form version of the PLOS Biology S4 Data supplement from Binsfeld
et al. (2025). These values reproduce the authors' score/Z-score
hit-calling workflow and can be compared with DStressR model-based calls
from
[binsfeld_reporter_auc](https://muellsen.github.io/DStressR/reference/binsfeld_reporter_auc.md).

## Usage

``` r
binsfeld_reporter_scores
```

## Format

A data frame with one row per well, strain, statistic, promoter, and
replicate. Columns are `well`, `drug`, `compound`,
`concentration_ug_ml`, `strain`, `statistic`, `promoter`, `replicate`,
and `value`.

## Source

Binsfeld et al. (2025), PLOS Biology, S4 Data.
