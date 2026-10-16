# Fit DStressR with staged statistical options

`fit_destress()` is the main DStressR entry point. By default it fits
the model-based analysis, but it can also run named presets for the
screen-specific median-polish and empty-vector-control analyses.

## Usage

``` r
fit_destress(
  assay,
  technical = NULL,
  empirical_bayes = TRUE,
  empty_vector_reporter = NULL,
  background_rank = 0,
  background_rank_max = 5,
  background_rank_permutations = 100,
  background_rank_threshold = 0.99,
  background_rank_seed = NULL,
  background_impute = c("column_mean", "global_mean", "zero"),
  normalization = NULL,
  testing = NULL,
  aggregation = NULL,
  adjustment = NULL,
  interaction = FALSE,
  preset = NULL,
  verbose = FALSE,
  ...
)
```

## Arguments

- assay:

  A `destress_assay` produced by
  [`prepare_assay()`](https://muellsen.github.io/DStressR/reference/prepare_assay.md)
  or a raw assay data frame for `normalization = "linear_model"`, or a
  long expression table for the compatibility presets.

- technical:

  Character vector of batch, plate, replicate, or other technical-factor
  columns to include.

- empirical_bayes:

  If `TRUE`, moderates standard errors toward an empirical prior
  variance with prior degrees of freedom estimated from the observed
  variance distribution, using each contrast's degrees of freedom.
  Reporter-wise specific effects combine raw total-effect variances with
  Welch-Satterthwaite degrees of freedom and are moderated once. This
  maps to `testing = "moderated_t"` for the model path; `FALSE` maps to
  `testing = "student_t"`.

- empty_vector_reporter:

  Optional reporter/control strain used as an empty-vector reporter in
  the model-based path. When supplied, its reference-relative
  perturbation effect is subtracted from every reporter's
  reference-relative perturbation effect before reporter-library
  centering. Its raw contrast variance is added to the reporter's raw
  contrast variance, with Welch-Satterthwaite degrees of freedom; the
  combined contrast variance is moderated once when requested. This
  effect-level adjustment is an alternative to response-level background
  calibration in
  [`prepare_assay()`](https://muellsen.github.io/DStressR/reference/prepare_assay.md).
  The two routes cannot be combined. Supported for
  `interaction = FALSE`.

- background_rank:

  Non-negative integer. The default `0` removes no latent background.
  Values `1` or `2` additionally subtract a low-rank background term
  from the reference-relative total-effect matrix before testing
  rank-adjusted total and reporter-specific residual effects. Use
  `"auto"` to choose the rank by permutation parallel analysis.

- background_rank_max:

  Maximum rank considered when `background_rank = "auto"`.

- background_rank_permutations:

  Number of permutations used when `background_rank = "auto"`.

- background_rank_threshold:

  Permutation reference quantile used when `background_rank = "auto"`.

- background_rank_seed:

  Optional random seed used when `background_rank = "auto"`.

- background_impute:

  Method used to fill missing entries of the effect matrix before
  singular-value decomposition. The default `"column_mean"` replaces
  missing entries by the observed mean of the corresponding perturbation
  column. `"global_mean"` and `"zero"` are also available.

- normalization:

  One of `"linear_model"`, `"median_polish"`, or `"empty_vector"`.
  `"model"` and `"evc"` are accepted aliases.

- testing:

  One of `"student_t"`, `"moderated_t"`, or `"gaussian_z"`.

- aggregation:

  One of `"none"` or `"max_p"`.

- adjustment:

  One of `"global"`, `"by_reporter"`, or `"none"`.

- interaction:

  If `FALSE`, fit one Gaussian linear model per reporter with the
  control perturbation as reference and the supplied technical
  covariates as design terms. The latter is the scalable path for
  reporter-specific perturbation effects. If `TRUE`, fit the full
  reporter-by-perturbation interaction model.

- preset:

  Optional named preset: `"model"`, `"median_polish"`, or
  `"empty_vector_control"`. The model-based workflow remains the default
  when omitted. `"median_polish_legacy"` is accepted as an alias for
  `"median_polish"`, and `"evc"` as an alias for
  `"empty_vector_control"`.

- verbose:

  If `TRUE`, print stage-level progress messages during model fitting.
  This is useful for large screens and full interaction models.

- ...:

  For `normalization = "linear_model"` with a raw data frame, arguments
  passed to
  [`prepare_assay()`](https://muellsen.github.io/DStressR/reference/prepare_assay.md),
  including `growth_exponent`. For compatibility presets, arguments
  passed to the selected engine.

## Value

A fitted DStressR object. The model path returns a `destress_fit`;
screen-specific presets return their corresponding workflow result
objects.

## Details

The staged options make the major statistical choices explicit:
normalization, test statistic and p-value calculation, replicate
aggregation, and p-value adjustment. Only implemented combinations are
accepted. For the model-based path, growth-response normalization is
performed upstream by
[`prepare_assay()`](https://muellsen.github.io/DStressR/reference/prepare_assay.md),
where `growth_exponent` can be fixed, estimated, or supplied as
reporter-specific values.
