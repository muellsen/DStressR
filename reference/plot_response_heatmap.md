# Heatmap of a DStressR promoter-by-compound response matrix

Creates a standard heatmap for normalized promoter-compound responses.
The default `value` is `specific_effect`, matching
[`results()`](https://muellsen.github.io/DStressR/reference/results.md),
but workflow tables can use columns such as
`destress_eb_effect_centered`.

## Usage

``` r
plot_response_heatmap(
  table,
  value = "specific_effect",
  promoter = "promoter",
  compound = "compound",
  compound_label = compound,
  show_compound_ids = TRUE,
  top_n_compounds = 160,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  clip_quantile = 0.98,
  title = "DStressR promoter-by-compound matrix",
  subtitle = NULL,
  xlab = "Compounds",
  ylab = "Promoters",
  low = "#2166AC",
  mid = "white",
  high = "#B2182B"
)
```

## Arguments

- table:

  A data frame with one row per promoter-compound pair.

- value:

  Numeric response/effect column to show in the heatmap.

- promoter, compound:

  Columns identifying promoters and compounds.

- compound_label:

  Optional human-readable compound-name column. Defaults to `compound`.

- show_compound_ids:

  If `TRUE`, append compound IDs in square brackets to compound labels.

- top_n_compounds:

  If finite, show only the top compounds by mean absolute response. Use
  `Inf` to show all compounds.

- cluster_rows, cluster_cols:

  If `TRUE`, hierarchically cluster promoters and/or compounds.

- clip_quantile:

  Quantile of absolute response values used to clip the color scale. Set
  to `1` to use the observed maximum.

- title, subtitle, xlab, ylab:

  Plot labels.

- low, mid, high:

  Colors for negative, zero, and positive responses.

## Value

A `ggplot` object. The plotted matrix is available as
`attr(plot, "response_matrix")`.
