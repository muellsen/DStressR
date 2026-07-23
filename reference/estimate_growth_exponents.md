# Estimate promoter-specific growth normalization exponents

Estimates how luminescence scales with growth in negative-control wells,
optionally adjusting for technical factors such as batch, plate, or
replicate:

## Usage

``` r
estimate_growth_exponents(
  data,
  promoter = "promoter",
  compound = "compound",
  lux = "lux",
  growth = "growth",
  covariates = NULL,
  numeric_covariates = NULL,
  controls = "DMSO",
  pseudocount = 1e-08,
  min_control_n = 8,
  shrink = TRUE,
  alpha_bounds = c(-2, 3)
)
```

## Arguments

- data:

  A data frame with one row per well.

- promoter, compound, lux, growth:

  Column names.

- covariates:

  Optional technical covariate column names to include as additive
  adjustment terms when estimating growth slopes. Only covariates with
  more than one observed value in the relevant control subset are used.

- numeric_covariates:

  Optional subset of `covariates` that should enter the growth-slope
  model as numeric variables. Other covariates are converted to factors.

- controls:

  Control values in `compound`, usually DMSO wells.

- pseudocount:

  Added before log2 transformation.

- min_control_n:

  Minimum control wells needed for a promoter-specific raw slope.
  Promoters with fewer controls use the global slope.

- shrink:

  If `TRUE`, shrink promoter-specific slopes toward the global control
  slope.

- alpha_bounds:

  Optional numeric length-2 bounds for the final exponent. Use `NULL`
  for no clipping.

## Value

A data frame with raw promoter intercepts and raw and shrunken growth
exponents per promoter.

## Details

\$\$\log_2(LUX_i) = a_g + \alpha_g \log_2(growth_i) + X_i\theta +
e_i\$\$

Raw promoter-specific slopes are then shrunk toward a global
control-well slope using an empirical-Bayes normal prior. The shrunken
\\\alpha_g\\ values can be used in
[`prepare_assay()`](https://muellsen.github.io/DStressR/reference/prepare_assay.md)
to compute:

\$\$y_i = \log_2(LUX_i) - \alpha_g \log_2(growth_i)\$\$
