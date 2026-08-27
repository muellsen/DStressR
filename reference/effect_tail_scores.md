# Robust diagnostic tail scores for effect distributions

Computes robust Gaussian tail probabilities for effect values,
optionally within strata. These values are diagnostics for unreplicated
or exploratory effect matrices and should not be interpreted as formal
p-values.

## Usage

``` r
effect_tail_scores(
  table,
  effect = "effect",
  group = NULL,
  min_n = 5,
  alternative = c("two_sided", "greater", "less")
)
```

## Arguments

- table:

  A data frame containing an effect column.

- effect:

  Numeric effect column.

- group:

  Optional grouping columns. Centers and scales are estimated separately
  within each group.

- min_n:

  Minimum number of finite effects required within a group.

- alternative:

  Tail alternative: `"two_sided"`, `"greater"`, or `"less"`.

## Value

The input table with diagnostic center, scale, z-score, tail
probability, and negative log10 tail score columns appended.
