# Extract model results

Extract model results

## Usage

``` r
results(fit, perturbations = NULL, reporters = NULL, verbose = FALSE)
```

## Arguments

- fit:

  A `destress_fit` or `destress_summary_fit` object.

- perturbations:

  Optional perturbation subset.

- reporters:

  Optional reporter subset.

- verbose:

  If `TRUE`, print stage-level progress messages while extracting model
  results.

## Value

A data frame with total and reporter-specific effects.
