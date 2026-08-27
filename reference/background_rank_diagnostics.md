# Diagnose low-rank background structure

Computes singular values of a reporter-by-perturbation effect matrix and
compares them with a permutation null. The default permutation shuffles
reporter labels within each perturbation, preserving the
perturbation-wise marginal distribution while breaking shared
reporter-loading structure.

## Usage

``` r
background_rank_diagnostics(
  table,
  effect = "total_effect",
  reporter = "reporter",
  perturbation = "perturbation",
  rank_max = 10,
  permutations = 100,
  threshold = 0.99,
  impute = "column_mean",
  seed = NULL
)
```

## Arguments

- table:

  Data frame with reporter, perturbation, and effect columns.

- effect:

  Numeric effect column to decompose, usually `total_effect` or
  `background_adjusted_effect`.

- reporter, perturbation:

  Column names identifying reporters and perturbations.

- rank_max:

  Maximum component index to report.

- permutations:

  Number of null permutations. Use `0` to skip the null.

- threshold:

  Permutation reference quantile used for automatic rank selection and
  reported as `null_threshold`.

- impute:

  Method used to fill missing matrix entries before singular-value
  decomposition.

- seed:

  Optional random seed for reproducible permutations.

## Value

A data frame with observed singular values, variance fractions, and
optional permutation summaries.
