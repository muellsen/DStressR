# Decompose a reporter-by-perturbation effect matrix

Estimates a low-rank background component from an effect matrix and
returns the observed effect, low-rank component, and rank-adjusted
residual in long format. If grouping columns are supplied, the
decomposition is performed independently within each group.

## Usage

``` r
low_rank_effect_decomposition(
  table,
  effect = "total_effect",
  reporter = "reporter",
  perturbation = "perturbation",
  group = NULL,
  rank = 1,
  impute = c("column_mean", "global_mean", "zero"),
  reporter_label = reporter,
  perturbation_label = perturbation
)
```

## Arguments

- table:

  A data frame with one row per reporter-perturbation pair.

- effect:

  Numeric effect column to decompose.

- reporter, perturbation:

  Columns identifying reporters and perturbations.

- group:

  Optional grouping columns. A separate low-rank decomposition is fitted
  within each group.

- rank:

  Non-negative rank of the background component.

- impute:

  Method used to fill missing entries before singular-value
  decomposition. The default `"column_mean"` replaces missing entries by
  the observed mean of the corresponding perturbation column.

- reporter_label, perturbation_label:

  Optional display-label columns.

## Value

A data frame containing the original effect, the low-rank effect, the
rank-adjusted effect, and rank-1 reporter/perturbation scores when
`rank >= 1`.
