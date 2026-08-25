# Perturbation-level response diagnostics

Summarizes a reporter-by-perturbation effect table into
perturbation-level diagnostics. The primary use is to compare the
absolute mean response of a perturbation with the variance of its
reporter-specific responses, analogous in spirit to mean-variance
diagnostic plots used for count-data workflows.

## Usage

``` r
perturbation_diagnostics(
  table,
  mean_effect = "total_effect",
  variance_effect = mean_effect,
  reporter = "reporter",
  perturbation = "perturbation",
  perturbation_label = perturbation,
  min_reporters = 2,
  trend_span = 0.45
)
```

## Arguments

- table:

  A data frame with one row per reporter-perturbation pair.

- mean_effect:

  Numeric column used to compute the perturbation-level mean effect. For
  DStressR results, `total_effect` or `rank_adjusted_total_effect` are
  typical choices.

- variance_effect:

  Numeric column whose cross-reporter variance is computed for each
  perturbation. Defaults to `mean_effect`.

- reporter, perturbation:

  Columns identifying reporters and perturbations.

- perturbation_label:

  Optional human-readable perturbation label column. Defaults to
  `perturbation`.

- min_reporters:

  Minimum number of finite reporter-level effects required for a
  perturbation.

- trend_span:

  Span used for the loess trend of variance over
  `rank(abs(mean effect))`.

## Value

A data frame with one row per perturbation.

## Details

The function is deliberately generic: columns are referred to as
reporters and perturbations, although DStressR result tables usually
contain reporters and perturbations. Large absolute mean effects with
low variance indicate coherent, broad responses, whereas unusually large
variance conditional on the absolute mean effect indicates heterogeneous
reporter behavior.
