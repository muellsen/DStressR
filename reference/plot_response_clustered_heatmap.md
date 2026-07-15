# Clustered heatmap with promoter and compound dendrograms

Draws a clustered promoter-by-compound response heatmap with
hierarchical trees on both axes. Unlike
[`plot_response_cluster_blocks`](https://muellsen.github.io/DStressR/reference/plot_response_cluster_blocks.md),
this keeps the individual matrix cells visible and uses the dendrograms
to reveal structure without collapsing the data into coarse blocks.

## Usage

``` r
plot_response_clustered_heatmap(
  table,
  value = "specific_effect",
  promoter = "promoter",
  compound = "compound",
  compound_label = compound,
  show_compound_ids = TRUE,
  top_n_compounds = 400,
  n_promoter_clusters = 6,
  n_compound_clusters = 14,
  missing_value = 0,
  clip_quantile = 0.98,
  file = NULL,
  width = 14,
  height = 8,
  res = 300,
  title = "DStressR clustered response heatmap",
  subtitle = NULL,
  show_rownames = TRUE,
  show_colnames = FALSE,
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

  If `TRUE`, append compound IDs in square brackets to compound labels
  before clustering.

- top_n_compounds:

  If finite, show only the top compounds by mean absolute response. Use
  `Inf` to show all compounds.

- n_promoter_clusters, n_compound_clusters:

  Number of dendrogram clusters returned in the cluster assignment
  tables.

- missing_value:

  Value used only for clustering missing matrix entries. Heatmap cells
  with missing values are left missing.

- clip_quantile:

  Quantile of absolute response values used to clip the color scale. Set
  to `1` to use the observed maximum.

- file:

  Optional output file. Supports `.png` and `.pdf`. If `NULL`, the plot
  is drawn on the active graphics device.

- width, height:

  Plot size in inches when `file` is supplied.

- res:

  PNG resolution in dots per inch.

- title, subtitle:

  Plot title and subtitle.

- show_rownames, show_colnames:

  Whether to draw row and column labels.

- low, mid, high:

  Colors for negative, zero, and positive responses.

## Value

Invisibly returns a list containing the response matrix, clustering
objects, ordered matrix, cluster assignments, and color limit.
