# Prepare a chemical-genomics assay table

Computes a growth-adjusted log2 promoter-activity response from
luminescence and growth measurements. By default, promoter-specific
growth exponents are estimated from control wells with available
technical-factor adjustment and shrunk toward a global control-well
slope. Set `growth_exponent = 1` to reproduce the current workflow's
log2(LUX / OD) response.

## Usage

``` r
prepare_assay(
  data,
  promoter = "promoter",
  compound = "compound",
  control = "DMSO",
  lux = "lux",
  growth = "growth",
  growth_exponent = "estimate",
  control_values = control,
  response = NULL,
  batch = NULL,
  plate = NULL,
  replicate = NULL,
  pseudocount = 1e-08
)
```

## Arguments

- data:

  A data frame with one row per promoter-compound-replicate well.

- promoter, compound:

  Column names identifying promoter and compound.

- control:

  Label in `compound` for the negative control, usually DMSO.

- lux, growth:

  Column names for luminescence and growth summaries.

- growth_exponent:

  Fixed coefficient for growth normalization, a named vector keyed by
  promoter, or `"estimate"` to estimate promoter-specific exponents from
  controls.

- control_values:

  Values in `compound` used as controls for growth exponent estimation.
  Defaults to `control`.

- response:

  Optional existing response column. If supplied, `lux` and `growth` are
  not used to compute the response.

- batch, plate, replicate:

  Optional technical-factor column names. When
  `growth_exponent = "estimate"`, these columns are also used as
  covariates while estimating promoter-specific growth exponents.

- pseudocount:

  Added before log2 transformation.

## Value

A `destress_assay` data frame.
