# Clustered block map of a DStressR promoter-by-compound response matrix

Hierarchically clusters promoters and compounds, cuts the dendrograms
into interpretable groups, and plots the mean response for each
promoter-cluster by compound-cluster block. This is useful as a compact
overview when the full compound library is too large for individual
compound labels.

## Usage

``` r
plot_response_cluster_blocks(
  table,
  value = "specific_effect",
  promoter = "promoter",
  compound = "compound",
  compound_label = compound,
  show_compound_ids = TRUE,
  n_promoter_clusters = 6,
  n_compound_clusters = 14,
  missing_value = 0,
  clip_quantile = 0.98,
  show_counts = TRUE,
  title = "DStressR clustered response map",
  subtitle = NULL,
  xlab = "Compound clusters",
  ylab = "Promoter clusters",
  low = "#2166AC",
  mid = "white",
  high = "#B2182B"
)
```

## Arguments

- table:

  A data frame with one row per promoter-compound pair.

- value:

  Numeric response/effect column to summarize.

- promoter, compound:

  Columns identifying promoters and compounds.

- compound_label:

  Optional human-readable compound-name column. Defaults to `compound`.

- show_compound_ids:

  If `TRUE`, append compound IDs in square brackets to compound labels
  before clustering.

- n_promoter_clusters, n_compound_clusters:

  Number of dendrogram clusters to use for promoters and compounds.

- missing_value:

  Value used only for clustering missing matrix entries. Block summaries
  are still computed from observed finite values.

- clip_quantile:

  Quantile of absolute block means used to clip the color scale. Set to
  `1` to use the observed maximum.

- show_counts:

  If `TRUE`, annotate each tile with the number of compounds in that
  compound cluster.

- title, subtitle, xlab, ylab:

  Plot labels.

- low, mid, high:

  Colors for negative, zero, and positive responses.

## Value

A `ggplot` object with attributes `response_matrix`,
`promoter_clusters`, `compound_clusters`, `block_summary`, `row_hclust`,
and `col_hclust`.
