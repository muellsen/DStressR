# Extract model results

Extract model results

## Usage

``` r
results(fit, perturbations = NULL, reporters = NULL, verbose = FALSE)
```

## Arguments

- fit:

  A `destress_fit` or `destress_summary_fit` object.

- perturbations:

  Optional perturbation subset.

- reporters:

  Optional reporter subset.

- verbose:

  If `TRUE`, print stage-level progress messages while extracting model
  results.

## Value

A data frame with total and reporter-relative effects. Reporter-wise
results also include `total_var_raw` and `relative_var` (before
moderation), `residual_df`, `relative_df_raw` (Welch-Satterthwaite
effective degrees of freedom), and `relative_df_for_test` (including
prior degrees of freedom when moderation is enabled). Full
interaction-model results include `total_var_raw`, `residual_df`, and
`total_df_for_test`. Their `relative_var` uses the centered contrast and
full coefficient covariance; `relative_df_raw` and
`relative_df_for_test` record its ordinary and moderated degrees of
freedom. Both targets honor `empirical_bayes`.
