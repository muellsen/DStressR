# Fit a three-part empirical-null mixture to promoter-compound effects

Fits, separately for each promoter, a three-component Student-t mixture
to adjusted promoter-compound effects. The ordered components are
interpreted as repressed, null, and activated effects. This second-stage
model is intended for empirical-null calibration after the first-stage
DStressR model has already adjusted growth, technical factors,
compound-wide effects, and promoter-specific variance.

## Usage

``` r
fit_effect_mixture(
  table,
  value = "specific_effect",
  promoter = "promoter",
  df = 4,
  max_iter = 2000,
  tol = 1e-06,
  min_scale = 1e-04,
  min_prior = 1e-04,
  padj_method = "BH"
)
```

## Arguments

- table:

  A data frame with one row per promoter-compound pair.

- value:

  Numeric effect column, usually `specific_effect` or a centered
  DStressR EB effect column.

- promoter:

  Column identifying promoters.

- df:

  Degrees of freedom for each Student-t component. Smaller values give
  heavier tails.

- max_iter:

  Maximum EM iterations per promoter.

- tol:

  Relative log-likelihood convergence tolerance.

- min_scale:

  Lower bound for component scale.

- min_prior:

  Lower bound for component mixing proportions.

- padj_method:

  Multiple-testing correction method passed to
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html), applied
  within promoter to empirical-null p-values.

## Value

The input table with posterior probabilities, local FDR, posterior
non-null probability, local-FDR q-values within promoter, empirical-null
p-values, promoter-wise adjusted p-values, and posterior class appended.
The promoter-level fitted parameters are available as
`attr(result, "mixture_summary")`.
