# Add DGrowthR-derived growth metrics to a DStressR assay table

DStressR is designed to work hand in hand with DGrowthR: DGrowthR models
the optical-density growth curves, and DStressR can use one chosen
DGrowthR growth metric as the growth covariate in
[`prepare_assay()`](https://muellsen.github.io/DStressR/reference/prepare_assay.md).

## Usage

``` r
add_dgrowthr_growth(
  data,
  object,
  by = "curve_id",
  model_covariate = by,
  growth_metric = "max_growth",
  output = "growth",
  keep_dgrowthr_columns = FALSE
)
```

## Arguments

- data:

  A data frame with one row per DStressR assay observation.

- object:

  A DGrowthR object after
  [`DGrowthR::estimate_growth_parameters()`](https://bio-datascience.github.io/DGrowthR/reference/estimate_growth_parameters.html).

- by:

  Column in `data` and `object@metadata` identifying the growth curve
  for each assay row. DGrowthR calls this `curve_id`.

- model_covariate:

  Metadata column used as DGrowthR's
  `estimate_growth_parameters(model_covariate = ...)`. Defaults to `by`,
  corresponding to one GP fit per growth curve.

- growth_metric:

  Column in `object@growth_parameters` to use as the growth measurement.
  Common choices are `max_growth`, `AUC`, `OD_16`, and `AUC_16`.

- output:

  Column name to create in `data`.

- keep_dgrowthr_columns:

  If `TRUE`, also keep the DGrowthR join key and the unrenamed metric
  column when possible.

## Value

`data` with an additional numeric growth column named by `output`.

## Details

This function is deliberately opt-in. The current hit-determination
scripts use the original exported growth summaries unless the analyst
explicitly calls this helper and passes the resulting column to
[`prepare_assay()`](https://muellsen.github.io/DStressR/reference/prepare_assay.md).
This keeps the present analysis reproducible while making the DGrowthR
bridge available for future comparisons.

This helper joins `object@growth_parameters` from a fitted DGrowthR
object onto a promoter/luminescence assay table. Run
[`DGrowthR::estimate_growth_parameters()`](https://bio-datascience.github.io/DGrowthR/reference/estimate_growth_parameters.html)
first, optionally with `od_auc_at_t` if you want time-specific columns
such as `OD_16` or `AUC_16`.
