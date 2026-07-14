# Fit DStressR with staged statistical options

`fit_destress()` is the main DStressR entry point. By default it fits
the model-based analysis, but it can also run named compatibility
presets for the legacy median-polish and Empty Vector Control analyses.

## Usage

``` r
fit_destress(
  assay,
  technical = NULL,
  empirical_bayes = TRUE,
  normalization = NULL,
  testing = NULL,
  aggregation = NULL,
  adjustment = NULL,
  preset = NULL,
  ...
)
```

## Arguments

- assay:

  A `destress_assay` produced by
  [`prepare_assay()`](https://bio-datascience.github.io/DStressR/reference/prepare_assay.md)
  or a raw assay data frame for `normalization = "linear_model"`, or a
  long expression table for the compatibility presets.

- technical:

  Character vector of batch, plate, replicate, or other technical-factor
  columns to include.

- empirical_bayes:

  If `TRUE`, lightly shrinks standard errors toward a common prior
  variance. This maps to `testing = "moderated_t"` for the model path;
  `FALSE` maps to `testing = "student_t"`.

- normalization:

  One of `"linear_model"`, `"median_polish"`, or `"empty_vector"`.
  `"model"` and `"evc"` are accepted aliases.

- testing:

  One of `"student_t"`, `"moderated_t"`, or `"gaussian_z"`.

- aggregation:

  One of `"none"` or `"max_p"`.

- adjustment:

  One of `"global"`, `"by_promoter"`, or `"none"`.

- preset:

  Optional named preset: `"model"`, `"median_polish_legacy"`, or
  `"empty_vector_control"`. Common aliases such as `"median_polish"` and
  `"evc"` are accepted.

- ...:

  For `normalization = "linear_model"` with a raw data frame, arguments
  passed to
  [`prepare_assay()`](https://bio-datascience.github.io/DStressR/reference/prepare_assay.md),
  including `growth_exponent`. For compatibility presets, arguments
  passed to the selected engine.

## Value

A fitted DStressR object. The model path returns a `destress_fit`;
compatibility presets return their corresponding legacy result objects.

## Details

The staged options make the major statistical choices explicit:
normalization, test statistic and p-value calculation, replicate
aggregation, and p-value adjustment. Only implemented combinations are
accepted. For the model-based path, growth-response normalization is
performed upstream by
[`prepare_assay()`](https://bio-datascience.github.io/DStressR/reference/prepare_assay.md),
where `growth_exponent` can be fixed, estimated, or supplied as
promoter-specific values.
