# Heatmap of significant DStressR hits

Shows a reporter-by-perturbation effect matrix with significant pairs
highlighted by color and non-significant pairs shown as a light
background. This plot is a compact companion to
[`plot_response_heatmap()`](https://muellsen.github.io/DStressR/reference/plot_response_heatmap.md)
for inspecting the discovered hit structure.

## Usage

``` r
plot_hit_heatmap(
  table,
  effect = "specific_effect",
  padj = "specific_padj_by_reporter",
  reporter = "reporter",
  perturbation = "perturbation",
  perturbation_label = perturbation,
  show_perturbation_ids = TRUE,
  top_n_perturbations = 160,
  fdr = 0.05,
  lfc = 0,
  drop_empty_perturbations = TRUE,
  reporter_order = NULL,
  order_rows = c("global", "cluster", "input", "frequency"),
  order_cols = c("frequency", "cluster", "input"),
  clip_quantile = 0.98,
  color_limit = NULL,
  show_perturbation_labels = NULL,
  top_perturbation_labels = 40,
  perturbation_label_score = NULL,
  perturbation_label_min_gap = NULL,
  perturbation_label_angle = 45,
  title = "DStressR significant-hit matrix",
  subtitle = NULL,
  xlab = "Perturbations",
  ylab = "Reporters",
  legend_title = effect,
  low = "#2166AC",
  mid = "white",
  high = "#B2182B"
)
```

## Arguments

- table:

  A data frame with one row per reporter-perturbation pair.

- effect:

  Effect-size column shown by color for significant hits.

- padj:

  Adjusted p-value column used for hit calls.

- reporter, perturbation:

  Columns identifying reporters and perturbations.

- perturbation_label:

  Optional human-readable perturbation-name column. Defaults to
  `perturbation`.

- show_perturbation_ids:

  If `TRUE`, append perturbation IDs in square brackets to perturbation
  labels.

- top_n_perturbations:

  If finite, show only perturbations with the strongest hit evidence,
  ranked by hit count and effect size. Use `Inf` to show all
  perturbations.

- fdr:

  FDR threshold for hit highlighting.

- lfc:

  Minimum absolute effect size for hit highlighting.

- drop_empty_perturbations:

  If `TRUE`, remove perturbations with no significant hits from the
  displayed matrix.

- reporter_order:

  Optional global reporter order used when `order_rows = "global"`. If
  omitted, the option `DStressR.reporter_order` is used when set;
  otherwise known DStressR paper reporters are shown in their manuscript
  order and remaining reporters are sorted alphabetically.

- order_rows, order_cols:

  Ordering strategy for reporters and perturbations. Use `"global"` for
  the package-wide reporter order, `"input"` to preserve the
  input/factor order, `"frequency"` to order by number of hits, or
  `"cluster"` for hierarchical clustering of the hit matrix.

- clip_quantile:

  Quantile of absolute significant effects used to clip the color scale.
  Set to `1` to use the observed maximum.

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

  Colorbar title. Defaults to the selected `effect` column.

- low, mid, high:

  Colors for negative, zero, and positive hit effects.

## Value

A `ggplot` object with attributes `hit_matrix`, `hit_summary`,
`plotted_pairs`, and `color_limit`.
