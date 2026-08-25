# Empirical p-values from replicate-averaged perturbation effects and DMSO nulls

Computes empirical p-values by first averaging replicate-level values
for each reporter-perturbation-stratum combination, then comparing each
averaged non-control perturbation value to the corresponding
distribution of averaged DMSO control values from the same reporter and
stratum, such as library plate. Optionally, it also computes Monte Carlo
permutation p-values by repeatedly drawing replicate-sized DMSO sets
from the same matched null stratum.

## Usage

``` r
empirical_replicate_pvalues(
  table,
  value,
  reporter = "reporter",
  perturbation = "perturbation",
  control,
  replicate = NULL,
  strata = NULL,
  min_replicates = 2,
  min_null = 5,
  permutation = FALSE,
  B = 1000,
  seed = NULL,
  alternative = c("two.sided", "greater", "less"),
  padj_method = "BH"
)
```

## Arguments

- table:

  Replicate-level data frame.

- value:

  Numeric value column to test, for example a DStressR adjusted effect
  such as `destress_eb_effect_centered`.

- reporter, perturbation:

  Columns identifying reporters and perturbations.

- control:

  Character vector of control perturbation IDs, usually DMSO wells.

- replicate:

  Optional replicate column. The function does not require replicate
  labels for averaging, but the argument documents the intended
  replicate-level input and is checked when supplied.

- strata:

  Optional columns defining matched null strata. Use
  `strata = "libplate"` to compare perturbations only to DMSO wells from
  the same library plate.

- min_replicates:

  Minimum finite replicate values required for an averaged perturbation
  or DMSO control value.

- min_null:

  Minimum number of averaged DMSO controls required to compute an
  empirical p-value.

- permutation:

  If `TRUE`, compute an additional permutation p-value by sampling
  replicate-sized DMSO sets within each matched null stratum.

- B:

  Number of permutation draws when `permutation = TRUE`.

- seed:

  Optional random seed for reproducible permutation p-values.

- alternative:

  One of `"two.sided"`, `"greater"`, or `"less"`.

- padj_method:

  Multiple-testing correction method passed to
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html), applied
  within reporter.

## Value

A data frame with one row per non-control reporter-perturbation-stratum
average, including empirical p-values and reporter-wise adjusted
p-values.
