# Fit a three-part empirical-null mixture to reporter-perturbation effects

Fits, separately for each reporter, a three-component Student-t mixture
to adjusted reporter-perturbation effects. The ordered components are
interpreted as repressed, null, and activated effects. This second-stage
model is intended for empirical-null calibration after the first-stage
DStressR model has already adjusted growth, technical factors,
perturbation-wide effects, and reporter-specific variance.

## Usage

``` r
fit_effect_mixture(
  table,
  value = "specific_effect",
  reporter = "reporter",
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

  A data frame with one row per reporter-perturbation pair.

- value:

  Numeric effect column, usually `specific_effect` or a centered
  DStressR EB effect column.

- reporter:

  Column identifying reporters.

- df:

  Degrees of freedom for each Student-t component. Smaller values give
  heavier tails.

- max_iter:

  Maximum EM iterations per reporter.

- tol:

  Relative log-likelihood convergence tolerance.

- min_scale:

  Lower bound for component scale.

- min_prior:

  Lower bound for component mixing proportions.

- padj_method:

  Multiple-testing correction method passed to
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html), applied
  within reporter to empirical-null p-values.

## Value

The input table with posterior probabilities, local FDR, empirical-null
p-values, reporter-wise adjusted p-values, and posterior class appended.
The reporter-level fitted parameters are available as
`attr(result, "mixture_summary")`.
