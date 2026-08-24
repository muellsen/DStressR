# Heatmap of significant DStressR hits

Shows a promoter-by-compound effect matrix with significant pairs
highlighted by color and non-significant pairs shown as a light
background. This plot is a compact companion to
[`plot_response_heatmap()`](https://muellsen.github.io/DStressR/reference/plot_response_heatmap.md)
for inspecting the discovered hit structure.

## Usage

``` r
plot_hit_heatmap(
  table,
  effect = "specific_effect",
  padj = "specific_padj_by_promoter",
  promoter = "promoter",
  compound = "compound",
  compound_label = compound,
  show_compound_ids = TRUE,
  top_n_compounds = 160,
  fdr = 0.05,
  lfc = 0,
  drop_empty_compounds = TRUE,
  promoter_order = NULL,
  order_rows = c("global", "cluster", "input", "frequency"),
  order_cols = c("frequency", "cluster", "input"),
  clip_quantile = 0.98,
  color_limit = NULL,
  show_compound_labels = NULL,
  top_compound_labels = 40,
  compound_label_score = NULL,
  compound_label_min_gap = NULL,
  compound_label_angle = 45,
  title = "DStressR significant-hit matrix",
  subtitle = NULL,
  xlab = "Compounds",
  ylab = "Promoters",
  legend_title = effect,
  low = "#2166AC",
  mid = "white",
  high = "#B2182B"
)
```

## Arguments

- table:

  A data frame with one row per promoter-compound pair.

- effect:

  Effect-size column shown by color for significant hits.

- padj:

  Adjusted p-value column used for hit calls.

- promoter, compound:

  Columns identifying promoters and compounds.

- compound_label:

  Optional human-readable compound-name column. Defaults to `compound`.

- show_compound_ids:

  If `TRUE`, append compound IDs in square brackets to compound labels.

- top_n_compounds:

  If finite, show only compounds with the strongest hit evidence, ranked
  by hit count and effect size. Use `Inf` to show all compounds.

- fdr:

  FDR threshold for hit highlighting.

- lfc:

  Minimum absolute effect size for hit highlighting.

- drop_empty_compounds:

  If `TRUE`, remove compounds with no significant hits from the
  displayed matrix.

- promoter_order:

  Optional global promoter order used when `order_rows = "global"`. If
  omitted, the option `DStressR.promoter_order` is used when set;
  otherwise known DStressR paper promoters are shown in their manuscript
  order and remaining promoters are sorted alphabetically.

- order_rows, order_cols:

  Ordering strategy for promoters and compounds. Use `"global"` for the
  package-wide promoter order, `"input"` to preserve the input/factor
  order, `"frequency"` to order by number of hits, or `"cluster"` for
  hierarchical clustering of the hit matrix.

- clip_quantile:

  Quantile of absolute significant effects used to clip the color scale.
  Set to `1` to use the observed maximum.

- color_limit:

  Optional positive color-scale limit. If supplied, values are clipped
  to `[-color_limit, color_limit]`; otherwise the limit is computed from
  `clip_quantile`.

- show_compound_labels:

  If `TRUE`, draw all x-axis compound labels. If `FALSE`, suppress
  x-axis compound labels. The default labels the `top_compound_labels`
  compounds with largest absolute column sums, or all compounds when
  fewer are plotted.

- top_compound_labels:

  Number of highest-signal compounds to label when
  `show_compound_labels = NULL`.

- compound_label_score:

  Optional numeric score used to choose the top-labelled compounds when
  `show_compound_labels = NULL`. If named, values are matched to
  compound labels; otherwise the order must match the displayed matrix
  columns.

- compound_label_min_gap:

  Minimum number of matrix columns between automatically selected
  labels. The default chooses a gap from the displayed matrix size and
  `top_compound_labels`.

- compound_label_angle:

  Angle used for visible compound labels.

- title, subtitle, xlab, ylab:

  Plot labels.

- legend_title:

  Colorbar title. Defaults to the selected `effect` column.

- low, mid, high:

  Colors for negative, zero, and positive hit effects.

## Value

A `ggplot` object with attributes `hit_matrix`, `hit_summary`,
`plotted_pairs`, and `color_limit`.
