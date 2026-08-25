# Clustered block map of a DStressR reporter-by-perturbation response matrix

Hierarchically clusters reporters and perturbations, cuts the
dendrograms into interpretable groups, and plots the mean response for
each reporter-cluster by perturbation-cluster block. This is useful as a
compact overview when the full perturbation library is too large for
individual perturbation labels.

## Usage

``` r
plot_response_cluster_blocks(
  table,
  value = "specific_effect",
  reporter = "reporter",
  perturbation = "perturbation",
  perturbation_label = perturbation,
  show_perturbation_ids = TRUE,
  n_reporter_clusters = 6,
  n_perturbation_clusters = 14,
  missing_value = 0,
  clip_quantile = 0.98,
  show_counts = TRUE,
  title = "DStressR clustered response map",
  subtitle = NULL,
  xlab = "Perturbation clusters",
  ylab = "Reporter clusters",
  low = "#2166AC",
  mid = "white",
  high = "#B2182B"
)
```

## Arguments

- table:

  A data frame with one row per reporter-perturbation pair.

- value:

  Numeric response/effect column to summarize.

- reporter, perturbation:

  Columns identifying reporters and perturbations.

- perturbation_label:

  Optional human-readable perturbation-name column. Defaults to
  `perturbation`.

- show_perturbation_ids:

  If `TRUE`, append perturbation IDs in square brackets to perturbation
  labels before clustering.

- n_reporter_clusters, n_perturbation_clusters:

  Number of dendrogram clusters to use for reporters and perturbations.

- missing_value:

  Value used only for clustering missing matrix entries. Block summaries
  are still computed from observed finite values.

- clip_quantile:

  Quantile of absolute block means used to clip the color scale. Set to
  `1` to use the observed maximum.

- show_counts:

  If `TRUE`, annotate each tile with the number of perturbations in that
  perturbation cluster.

- title, subtitle, xlab, ylab:

  Plot labels.

- low, mid, high:

  Colors for negative, zero, and positive responses.

## Value

A `ggplot` object with attributes `response_matrix`,
`reporter_clusters`, `perturbation_clusters`, `block_summary`,
`row_hclust`, and `col_hclust`.
