# Histogram of DStressR reporter-perturbation effects

Shows the empirical distribution of normalized reporter-perturbation
effects, either over all matrix entries or faceted by reporter.

## Usage

``` r
plot_effect_histogram(
  table,
  value = "specific_effect",
  reporter = "reporter",
  by = c("all", "reporter"),
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

  A data frame with one row per reporter-perturbation pair.

- value:

  Numeric effect column to plot.

- reporter:

  Column identifying reporters, used when `by = "reporter"`.

- by:

  Plot one pooled histogram (`"all"`) or reporter-faceted histograms
  (`"reporter"`).

- bins:

  Number of histogram bins.

- xlim:

  Optional two-element x-axis limit.

- scales:

  Facet scale behavior for `by = "reporter"`.

- title, subtitle, xlab, ylab:

  Plot labels.

- fill, border:

  Histogram fill and border colors.

## Value

A `ggplot` object.
