# Fit DStressR from supplied response estimates and variances

`fit_destress_summary()` analyzes data where each reporter-perturbation
combination is already summarized by a response estimate and an
estimated sampling variance. This is useful for public datasets that
release processed promoter activities with uncertainty estimates, but
not replicate-level wells. The resulting p-values are summary-level
Gaussian approximations rather than replicate-level moderated DStressR
tests.

## Usage

``` r
fit_destress_summary(
  data,
  reporter = "reporter",
  perturbation = "perturbation",
  response = "response",
  response_variance = "response_variance",
  response_se = NULL,
  reference_perturbation,
  global_method = c("weighted_mean", "mean"),
  background_rank = 0,
  background_rank_max = 5,
  background_rank_permutations = 100,
  background_rank_threshold = 0.99,
  background_rank_seed = NULL,
  background_impute = c("column_mean", "global_mean", "zero"),
  adjustment = c("global", "by_reporter", "none")
)
```

## Arguments

- data:

  A data frame with one row per reporter and perturbation.

- reporter, perturbation:

  Columns identifying reporters and perturbations.

- response:

  Numeric response estimate column.

- response_variance:

  Numeric response-variance column. This is the primary uncertainty
  input.

- response_se:

  Optional response-standard-error column. Used only when
  `response_variance` is `NULL`, in which case it is squared internally.

- reference_perturbation:

  Reference perturbation/condition.

- global_method:

  How to estimate perturbation-level global effects. `"weighted_mean"`
  uses inverse-variance weights; `"mean"` uses equal weights.

- background_rank:

  Non-negative integer or `"auto"` for permutation-based rank selection.

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

- adjustment:

  One of `"global"`, `"by_reporter"`, or `"none"`.

## Value

A `destress_summary_fit` object. Use
[`results()`](https://muellsen.github.io/DStressR/reference/results.md)
to extract the effect table.
