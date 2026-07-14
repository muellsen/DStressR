# Histogram of DStressR promoter-compound effects

Shows the empirical distribution of normalized promoter-compound
effects, either over all matrix entries or faceted by promoter.

## Usage

``` r
plot_effect_histogram(
  table,
  value = "specific_effect",
  promoter = "promoter",
  by = c("all", "promoter"),
  bins = 80,
  xlim = NULL,
  scales = "fixed",
  title = NULL,
  subtitle = NULL,
  xlab = NULL,
  ylab = "Count",
  fill = "#4E79A7",
  border = "white"
)
```

## Arguments

- table:

  A data frame with one row per promoter-compound pair.

- value:

  Numeric effect column to plot.

- promoter:

  Column identifying promoters, used when `by = "promoter"`.

- by:

  Plot one pooled histogram (`"all"`) or promoter-faceted histograms
  (`"promoter"`).

- bins:

  Number of histogram bins.

- xlim:

  Optional two-element x-axis limit.

- scales:

  Facet scale behavior for `by = "promoter"`.

- title, subtitle, xlab, ylab:

  Plot labels.

- fill, border:

  Histogram fill and border colors.

## Value

A `ggplot` object.
