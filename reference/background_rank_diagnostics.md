# Diagnose low-rank background structure

Computes singular values of a promoter-by-compound effect matrix and
compares them with a permutation null. The default permutation shuffles
promoter labels within each compound, preserving the compound-wise
marginal distribution while breaking shared promoter-loading structure.

## Usage

``` r
background_rank_diagnostics(
  table,
  effect = "specific_effect",
  promoter = "promoter",
  compound = "compound",
  rank_max = 10,
  permutations = 100,
  seed = NULL
)
```

## Arguments

- table:

  Data frame with promoter, compound, and effect columns.

- effect:

  Numeric effect column to decompose, usually `specific_effect`.

- promoter, compound:

  Column names identifying promoters and compounds.

- rank_max:

  Maximum component index to report.

- permutations:

  Number of null permutations. Use `0` to skip the null.

- seed:

  Optional random seed for reproducible permutations.

## Value

A data frame with observed singular values, variance fractions, and
optional permutation summaries.
