# Simulate a chemical-genomics screen

Simulate a chemical-genomics screen

## Usage

``` r
simulate_screen(
  n_promoters = 12,
  n_compounds = 24,
  n_replicates = 2,
  sigma = 0.15,
  seed = NULL
)
```

## Arguments

- n_promoters, n_compounds:

  Dimensions excluding DMSO.

- n_replicates:

  Number of technical replicates.

- sigma:

  Observation noise standard deviation.

- seed:

  Optional random seed.

## Value

A data frame suitable for
[`prepare_assay()`](https://bio-datascience.github.io/DStressR/reference/prepare_assay.md).
