# Plot diagnostic effect-tail histograms

Shows effect distributions and highlights observations in the diagnostic
tails. The plot is intended for exploratory effect matrices, especially
when formal replicate-based inference is unavailable.

## Usage

``` r
plot_effect_tail_histogram(
  table,
  effect = "effect",
  tail_probability = "tail_probability",
  facet = NULL,
  tail_threshold = 0.05,
  bins = 40,
  title = NULL,
  subtitle = NULL,
  xlab = NULL,
  ylab = "Effects"
)
```

## Arguments

- table:

  A data frame, usually returned by
  [`effect_tail_scores()`](https://muellsen.github.io/DStressR/reference/effect_tail_scores.md).

- effect:

  Numeric effect column to plot.

- tail_probability:

  Diagnostic tail-probability column.

- facet:

  Optional columns used for faceting.

- tail_threshold:

  Tail-probability threshold used for highlighting.

- bins:

  Number of histogram bins.

- title, subtitle, xlab, ylab:

  Plot labels.

## Value

A `ggplot` object.
