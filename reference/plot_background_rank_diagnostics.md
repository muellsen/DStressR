# Plot low-rank background diagnostics

Draws observed singular values together with permutation-null summaries
from
[`background_rank_diagnostics()`](https://muellsen.github.io/DStressR/reference/background_rank_diagnostics.md).
This is a diagnostic for choosing whether a low-rank background
component is visible in the effect matrix.

## Usage

``` r
plot_background_rank_diagnostics(
  table = NULL,
  diagnostics = NULL,
  effect = "total_effect",
  reporter = "reporter",
  perturbation = "perturbation",
  group = NULL,
  rank_max = 10,
  permutations = 100,
  threshold = 0.99,
  impute = c("column_mean", "global_mean", "zero"),
  seed = NULL,
  title = NULL,
  subtitle = NULL,
  xlab = "Component",
  ylab = "Singular value"
)
```

## Arguments

- table:

  Optional effect table. Ignored when `diagnostics` is supplied.

- diagnostics:

  Optional data frame returned by
  [`background_rank_diagnostics()`](https://muellsen.github.io/DStressR/reference/background_rank_diagnostics.md).

- effect:

  Numeric effect column to decompose, usually `total_effect` or
  `background_adjusted_effect`.

- reporter, perturbation:

  Column names identifying reporters and perturbations.

- group:

  Optional grouping columns. A separate diagnostic is computed and
  faceted for each group.

- rank_max:

  Maximum component index to report.

- permutations:

  Number of null permutations. Use `0` to skip the null.

- threshold:

  Permutation reference quantile.

- impute:

  Method used to fill missing entries before singular-value
  decomposition.

- seed:

  Optional random seed for reproducible permutations.

- title, subtitle, xlab, ylab:

  Plot labels.

## Value

A `ggplot` object. The diagnostic table is available as
`attr(plot, "diagnostics")`.
