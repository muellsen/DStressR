# Heatmap of a DStressR reporter-by-perturbation response matrix

Creates a standard heatmap for normalized reporter-perturbation
responses. The default `value` is `specific_effect`, matching
[`results()`](https://muellsen.github.io/DStressR/reference/results.md),
but workflow tables can use columns such as
`destress_eb_effect_centered`.

## Usage

``` r
plot_response_heatmap(
  table,
  value = "specific_effect",
  reporter = "reporter",
  perturbation = "perturbation",
  perturbation_label = perturbation,
  show_perturbation_ids = TRUE,
  top_n_perturbations = 160,
  reporter_order = NULL,
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  clip_quantile = 0.98,
  color_limit = NULL,
  show_perturbation_labels = NULL,
  top_perturbation_labels = 40,
  perturbation_label_score = NULL,
  perturbation_label_min_gap = NULL,
  perturbation_label_angle = 45,
  title = "DStressR reporter-by-perturbation matrix",
  subtitle = NULL,
  xlab = "Perturbations",
  ylab = "Reporters",
  legend_title = value,
  low = "#2166AC",
  mid = "white",
  high = "#B2182B"
)
```

## Arguments

- table:

  A data frame with one row per reporter-perturbation pair.

- value:

  Numeric response/effect column to show in the heatmap.

- reporter, perturbation:

  Columns identifying reporters and perturbations.

- perturbation_label:

  Optional human-readable perturbation-name column. Defaults to
  `perturbation`.

- show_perturbation_ids:

  If `TRUE`, append perturbation IDs in square brackets to perturbation
  labels.

- top_n_perturbations:

  If finite, show only the top perturbations by mean absolute response.
  Use `Inf` to show all perturbations.

- reporter_order:

  Optional global reporter order used when `cluster_rows = FALSE`. If
  omitted, the option `DStressR.reporter_order` is used when set;
  otherwise known DStressR paper reporters are shown in their manuscript
  order and remaining reporters are sorted alphabetically.

- cluster_rows, cluster_cols:

  If `TRUE`, hierarchically cluster reporters and/or perturbations.

- clip_quantile:

  Quantile of absolute response values used to clip the color scale. Set
  to `1` to use the observed maximum.

- color_limit:

  Optional positive color-scale limit. If supplied, values are clipped
  to `[-color_limit, color_limit]`; otherwise the limit is computed from
  `clip_quantile`.

- show_perturbation_labels:

  If `TRUE`, draw all x-axis perturbation labels. If `FALSE`, suppress
  x-axis perturbation labels. The default labels the
  `top_perturbation_labels` perturbations with largest absolute column
  sums, or all perturbations when fewer are plotted.

- top_perturbation_labels:

  Number of highest-signal perturbations to label when
  `show_perturbation_labels = NULL`.

- perturbation_label_score:

  Optional numeric score used to choose the top-labelled perturbations
  when `show_perturbation_labels = NULL`. If named, values are matched
  to perturbation labels; otherwise the order must match the displayed
  matrix columns.

- perturbation_label_min_gap:

  Minimum number of matrix columns between automatically selected
  labels. The default chooses a gap from the displayed matrix size and
  `top_perturbation_labels`.

- perturbation_label_angle:

  Angle used for visible perturbation labels.

- title, subtitle, xlab, ylab:

  Plot labels.

- legend_title:

  Colorbar title. Defaults to the selected `value` column.

- low, mid, high:

  Colors for negative, zero, and positive responses.

## Value

A `ggplot` object. The plotted matrix is available as
`attr(plot, "response_matrix")`.
