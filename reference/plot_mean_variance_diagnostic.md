# Plot perturbation-level mean-variance diagnostics

Draws a DESeq-style diagnostic plot in which perturbations are ordered
by their absolute mean response and the y-axis shows the variance of
reporter-level effects. A non-parametric trend can be overlaid, and the
most heterogeneous perturbations can be labelled.

## Usage

``` r
plot_mean_variance_diagnostic(
  table = NULL,
  diagnostics = NULL,
  mean_effect = "total_effect",
  variance_effect = mean_effect,
  reporter = "reporter",
  perturbation = "perturbation",
  perturbation_label = perturbation,
  min_reporters = 2,
  trend_span = 0.45,
  add_trend = TRUE,
  label_by = c("residual", "variance", "abs_mean", "none"),
  top_n = 8,
  point_size = 1.8,
  color_limits = NULL,
  color_quantile = 0.92,
  title = NULL,
  subtitle = NULL,
  xlab = NULL,
  ylab = NULL,
  legend_position = c(0.035, 0.965),
  legend_justification = c(0, 1)
)
```

## Arguments

- table:

  A data frame with one row per reporter-perturbation pair.

- diagnostics:

  Optional output from
  [`perturbation_diagnostics()`](https://muellsen.github.io/DStressR/reference/perturbation_diagnostics.md).
  If supplied, `table` is ignored.

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

- add_trend:

  If `TRUE`, overlay the loess trend.

- label_by:

  Which perturbations to label. `"residual"` labels perturbations with
  the largest variance residual above the trend, `"variance"` labels the
  largest variances, `"abs_mean"` labels the largest absolute mean
  effects, and `"none"` suppresses labels.

- top_n:

  Number of perturbations to label.

- point_size:

  Point size.

- color_limits:

  Optional numeric vector of length two for the color scale. If omitted,
  symmetric robust limits are computed from `color_quantile`.

- color_quantile:

  Quantile of the absolute mean effect used for robust color limits when
  `color_limits` is omitted.

- title, subtitle, xlab, ylab:

  Plot labels.

- legend_position, legend_justification:

  Numeric vectors passed to
  [`ggplot2::theme()`](https://ggplot2.tidyverse.org/reference/theme.html)
  to place the legend inside the plotting region.

## Value

A `ggplot` object. The diagnostic table is available as
`attr(plot, "diagnostics")`.
