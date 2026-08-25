# Prepare a chemical-genomics assay table

Computes a growth-adjusted log2 reporter-activity response from
luminescence and growth measurements. By default, reporter-specific
growth exponents are estimated from control wells with available
technical-factor adjustment and shrunk toward a global control-well
slope. Set `growth_exponent = 1` to reproduce the current workflow's
log2(LUX / OD) response.

## Usage

``` r
prepare_assay(
  data,
  reporter = "reporter",
  perturbation = "perturbation",
  control = "DMSO",
  lux = "lux",
  growth = "growth",
  growth_exponent = "estimate",
  control_values = control,
  response = NULL,
  batch = NULL,
  plate = NULL,
  replicate = NULL,
  growth_covariates = NULL,
  numeric_covariates = NULL,
  background_reporter = NULL,
  background_method = c("none", "subtract", "lm", "huber"),
  background_by = NULL,
  pseudocount = 1e-08
)
```

## Arguments

- data:

  A data frame with one row per reporter-perturbation-replicate well.

- reporter, perturbation:

  Column names identifying reporter and perturbation.

- control:

  Label or labels in `perturbation` for the negative control, usually
  DMSO. If several labels are supplied, they are pooled into the first
  label as the model reference while the original input column is kept
  unchanged.

- lux, growth:

  Column names for luminescence and growth summaries.

- growth_exponent:

  Fixed coefficient for growth normalization, a named vector keyed by
  reporter, or `"estimate"` to estimate reporter-specific exponents from
  controls.

- control_values:

  Values in `perturbation` used as controls for growth exponent
  estimation. Defaults to `control`.

- response:

  Optional existing response column. If supplied, `lux` and `growth` are
  not used to compute the response.

- batch, plate, replicate:

  Optional technical-factor column names. When
  `growth_exponent = "estimate"`, these columns are also used as
  covariates while estimating reporter-specific growth exponents unless
  `growth_covariates` is supplied.

- growth_covariates:

  Optional technical covariate column names used only while estimating
  reporter-specific growth exponents from control wells. If `NULL`,
  DStressR uses the supplied `batch`, `plate`, and `replicate` columns
  for backwards compatibility.

- numeric_covariates:

  Optional subset of technical covariate column names that should remain
  numeric in model matrices. Other optional covariates are converted to
  factors.

- background_reporter:

  Optional background reporter used as a background reference, e.g. an
  Empty Vector Control. When supplied, the default background method is
  `"huber"`. The background reporter is matched to other reporters by
  `background_by`, the response is calibrated, and the background
  reporter is excluded from model-based testing.

- background_method:

  One of `"none"`, `"subtract"`, `"lm"`, or `"huber"`. If omitted,
  DStressR uses `"none"` when no `background_reporter` is supplied and
  `"huber"` when one is supplied. `"subtract"` subtracts the matched
  background response. `"lm"` and `"huber"` replace each non-background
  reporter response by residuals from a reporter-wise calibration
  against the matched background.

- background_by:

  Columns used to match each observation to the background reporter.
  Defaults to `perturbation` plus the supplied technical columns.

- pseudocount:

  Added before log2 transformation.

## Value

A `destress_assay` data frame.
