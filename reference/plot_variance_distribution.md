# Plot diagnostic distributions for perturbation-level variances

Shows the empirical density of log10 variances together with fitted
diagnostic distributions from
[`fit_variance_distribution()`](https://muellsen.github.io/DStressR/reference/fit_variance_distribution.md).

## Usage

``` r
plot_variance_distribution(
  diagnostics,
  fit = NULL,
  variance = "effect_variance",
  distributions = c("beta_prime", "inverse_gamma", "log_normal", "log_t"),
  adjust = 1.15,
  grid_n = 512,
  title = NULL,
  subtitle = NULL,
  xlab = NULL,
  ylab = "Density"
)
```

## Arguments

- diagnostics:

  A data frame, usually from
  [`perturbation_diagnostics()`](https://muellsen.github.io/DStressR/reference/perturbation_diagnostics.md).

- fit:

  Optional output from
  [`fit_variance_distribution()`](https://muellsen.github.io/DStressR/reference/fit_variance_distribution.md).
  If omitted, distributions are fitted before plotting.

- variance:

  Numeric column containing positive variance estimates.

- distributions:

  Character vector passed to
  [`fit_variance_distribution()`](https://muellsen.github.io/DStressR/reference/fit_variance_distribution.md)
  when `fit` is omitted.

- adjust:

  Bandwidth adjustment for the empirical density.

- grid_n:

  Number of grid points for fitted curves.

- title, subtitle, xlab, ylab:

  Plot labels.

## Value

A `ggplot` object. The fit table is available as `attr(plot, "fit")`.
