# Reproduce the legacy median-polish workflow

This function implements the original median-polish hit-determination
workflow used for the Campylobacter reporter-library screen. It starts
from a long expression table, centers each
reporter-library-plate-replicate group by its DMSO wells, applies
[`stats::medpolish()`](https://rdrr.io/r/stats/medpolish.html) to the
resulting reporter-libplate-replicate by perturbation matrix, and
computes z-test p-values from the polished DMSO residual distribution.

## Usage

``` r
fit_median_polish(
  data,
  reporter = "reporter",
  perturbation = "srn_code",
  libplate = "libplate",
  replicate = "replicate",
  response = "log2.auc.16hmeasured.normed",
  control,
  exclude = character(),
  fdr = 0.05,
  normality = FALSE,
  normality_methods = c("shapiro", "lilliefors"),
  maxiter = 1000,
  eps = 1e-08
)
```

## Arguments

- data:

  Long expression table with one row per reporter-perturbation-replicate
  observation.

- reporter, perturbation, libplate, replicate:

  Column names identifying the reporter, perturbation/library well,
  library plate, and replicate.

- response:

  Column containing the already growth-normalized log2 response, for
  example `log2.auc.16hmeasured.normed`.

- control:

  Character vector of perturbation/library-well IDs used as DMSO
  controls.

- exclude:

  Character vector of perturbation/library-well IDs to remove before
  median polishing and hit calling, for example noisy DMSO wells.

- fdr:

  FDR threshold used to assign the `hit` class in the pair-level table.

- normality:

  If `TRUE`, test pre-polish DMSO-centered fold changes within each
  reporter-library-plate-replicate group.

- normality_methods:

  Character vector containing `"shapiro"` and/or `"lilliefors"`. The
  Lilliefors test requires the suggested `nortest` package.

- maxiter, eps:

  Passed to
  [`stats::medpolish()`](https://rdrr.io/r/stats/medpolish.html).

## Value

A list of class `destress_median_polish` with `replicate_results`,
`pair_results`, `polished_matrix`, `medpolish`, and optional
`normality_results` components.

## Details

The reporter-perturbation hit table follows the original conservative
replicate aggregation: DMSO and excluded control wells are removed, the
largest replicate-level p-value is retained for each
reporter-perturbation pair, and p-values are BH-adjusted within
reporter.
