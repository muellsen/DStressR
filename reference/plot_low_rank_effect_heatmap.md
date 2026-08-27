# Plot a low-rank effect decomposition

Shows the observed effect matrix, the estimated low-rank component, and
the rank-adjusted residual matrix. Rows are ordered by the first
reporter score within each group, making broad background axes visible
even when reporter labels are too dense to display.

## Usage

``` r
plot_low_rank_effect_heatmap(
  decomposition,
  matrices = c("effect", "low_rank_effect", "rank_adjusted_effect"),
  group = NULL,
  show_reporter_labels = FALSE,
  perturbation_order = NULL,
  clip_quantile = 0.985,
  color_limit = NULL,
  title = NULL,
  subtitle = NULL,
  xlab = "Perturbation",
  ylab = "Reporters ordered by rank-1 score",
  legend_title = "Effect"
)
```

## Arguments

- decomposition:

  Output from
  [`low_rank_effect_decomposition()`](https://muellsen.github.io/DStressR/reference/low_rank_effect_decomposition.md).

- matrices:

  Which matrices to display.

- group:

  Optional grouping columns used for faceting.

- show_reporter_labels:

  If `TRUE`, draw reporter labels.

- perturbation_order:

  Optional perturbation order.

- clip_quantile:

  Quantile of absolute effects used for color clipping.

- color_limit:

  Optional positive color limit.

- title, subtitle, xlab, ylab, legend_title:

  Plot labels.

## Value

A `ggplot` object. The plotted data are available as
`attr(plot, "plot_data")`.
