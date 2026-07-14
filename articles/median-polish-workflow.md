# Original median-polish workflow

`DStressR` includes a legacy-compatible median-polish path so that the
original workflow can be reproduced before moving to the model-based
DStressR analyses. This is useful for regression testing, method
comparison, and explaining which parts of the new framework are
replacements rather than silent changes.

The named workflow entry point is
`fit_workflow(..., workflow = "median_polish")`. The lower-level
[`fit_median_polish()`](https://bio-datascience.github.io/DStressR/reference/fit_median_polish.md)
function remains available for existing scripts.

## Required input files

For the original Campylobacter promoter-library workflow, the
median-polish baseline starts from two exported files.

### `expression_values.tsv.gz`

This is a long table with one row per measured
promoter-library-well-replicate observation. Its shape is approximately:

``` text
n_promoters x (n_compound_wells + n_dmso_wells + n_noisy_dmso_wells) x n_replicates
```

The required columns are:

- promoter identifier, for example `promoter`
- compound/library-well identifier, for example `srn_code`
- library plate, for example `libplate`
- replicate, for example `replicate`
- growth-normalized log2 luminescence response, originally
  `log2.auc.16hmeasured.normed`

If `srn_code` is not already present,
[`read_campylobacter_expression()`](https://bio-datascience.github.io/DStressR/reference/read_campylobacter_expression.md)
can reconstruct it from `libplate` and `well`.

### `LibMap.txt`

This is a library annotation table with one row per compound/library
well. The required columns are:

- `Library plate`
- `Well`
- `ProductName`
- `Catalog Number`

The original workflow identifies ordinary DMSO controls and noisy DMSO
wells from this table. DStressR uses the same compound key:

``` r

srn_code = paste0("lp", `Library plate`, "_", Well)
```

## Reproduce the original workflow

``` r

library(DStressR)

expression_df <- read_campylobacter_expression(
  expression_file = "expression_values.tsv.gz",
  libmap_file = "LibMap.txt"
)

libmap <- read.delim("LibMap.txt", check.names = FALSE)
libmap$srn_code <- paste0("lp", libmap[["Library plate"]], "_", libmap[["Well"]])

dmso_srn_codes <- libmap$srn_code[libmap$ProductName == "DMSO"]
dmso_noisy_srn_codes <- libmap$srn_code[libmap$ProductName == "DMSO noisy"]

legacy <- fit_workflow(
  expression_df,
  workflow = "median_polish",
  promoter = "promoter",
  compound = "srn_code",
  libplate = "libplate",
  replicate = "replicate",
  response = "log2.auc.16hmeasured.normed",
  control = dmso_srn_codes,
  exclude = dmso_noisy_srn_codes,
  normality = TRUE,
  maxiter = 1000,
  eps = 1e-8
)
```

The optional `normality_results` table corresponds to the original DMSO
normality check: Shapiro-Wilk and Lilliefors p-values are computed from
pre-polish DMSO-centered fold changes within each
promoter-libplate-replicate group, then BH-adjusted across groups.

The replicate-level table corresponds to the original
`expression_df.pvalues.tsv.gz` output:

``` r

replicate_pvalues <- legacy$replicate_results

head(replicate_pvalues)
```

The pair-level table corresponds to the original hit-calling table: DMSO
and excluded controls are removed, the largest replicate-level p-value
is retained for each promoter-compound pair, and BH correction is
applied within promoter.

``` r

hit_table <- legacy$pair_results

head(hit_table)
table(hit_table$hit)
```

## What the function does

`fit_workflow(..., workflow = "median_polish")` follows the original
hit-determination steps:

1.  Build
    `promoter_libplate_replicate = paste(promoter, libplate, replicate, sep = "_")`.
2.  Compute the DMSO mean normalized response within each
    promoter-libplate-replicate group.
3.  Compute `log2FC = normalized.lux - avg.dmso.expression`.
4.  Form a matrix with promoter-libplate-replicate groups in rows and
    `srn_code` compound/library wells in columns.
5.  Apply
    `stats::medpolish(..., na.rm = TRUE, maxiter = 1000, eps = 1e-8)`.
6.  Use the median-polish residuals as `log2FC.polished`.
7.  Estimate the DMSO residual mean and standard deviation within each
    promoter-libplate-replicate group.
8.  Compute
    `zscore = (log2FC.polished - dmso.avg_dmsoFC) / dmso.stdv_dmsoFC`.
9.  Compute two-sided Gaussian p-values.
10. For hit calling, remove DMSO/noisy-DMSO wells, keep the largest
    p-value across replicates for each promoter-compound pair, and apply
    BH correction within promoter.

If `normality = TRUE`, the function also runs the workflow’s DMSO
normality screen before median polishing. The Lilliefors test uses
[`nortest::lillie.test()`](https://rdrr.io/pkg/nortest/man/lillie.test.html);
install `nortest` or set `normality_methods = "shapiro"` when only the
Shapiro-Wilk check is needed.

## Minimal executable example

This toy example is deliberately tiny; it only shows the mechanics of
the legacy path. Real screens need multiple DMSO wells per
promoter-libplate-replicate group to estimate the DMSO null
distribution.

``` r

toy <- expand.grid(
  promoter = c("P1", "P2"),
  libplate = "lp1",
  replicate = c("r1", "r2"),
  srn_code = c("DMSO1", "DMSO2", "C1", "C2"),
  stringsAsFactors = FALSE
)

toy$log2.auc.16hmeasured.normed <- c(
  10.0, 10.2, 11.4, 9.4,
  10.1, 10.3, 11.5, 9.5,
  12.0, 12.2, 13.6, 11.6,
  12.1, 12.3, 13.7, 11.7
)

legacy <- fit_workflow(
  toy,
  workflow = "median_polish",
  control = c("DMSO1", "DMSO2")
)

legacy$replicate_results
#>    promoter_libplate_replicate promoter libplate replicate srn_code
#> 1                    P1_lp1_r1       P1      lp1        r1       C1
#> 2                    P1_lp1_r1       P1      lp1        r1       C2
#> 3                    P1_lp1_r1       P1      lp1        r1    DMSO1
#> 4                    P1_lp1_r1       P1      lp1        r1    DMSO2
#> 5                    P1_lp1_r2       P1      lp1        r2       C1
#> 6                    P1_lp1_r2       P1      lp1        r2       C2
#> 7                    P1_lp1_r2       P1      lp1        r2    DMSO1
#> 8                    P1_lp1_r2       P1      lp1        r2    DMSO2
#> 9                    P2_lp1_r1       P2      lp1        r1       C1
#> 10                   P2_lp1_r1       P2      lp1        r1       C2
#> 11                   P2_lp1_r1       P2      lp1        r1    DMSO1
#> 12                   P2_lp1_r1       P2      lp1        r1    DMSO2
#> 13                   P2_lp1_r2       P2      lp1        r2       C1
#> 14                   P2_lp1_r2       P2      lp1        r2       C2
#> 15                   P2_lp1_r2       P2      lp1        r2    DMSO1
#> 16                   P2_lp1_r2       P2      lp1        r2    DMSO2
#>    log2FC.polished dmso.avg_dmsoFC dmso.stdv_dmsoFC n        zscore    pvalue
#> 1            -0.05            0.05     3.140185e-16 2 -3.184526e+14 0.0000000
#> 2            -0.05            0.05     3.140185e-16 2 -3.184526e+14 0.0000000
#> 3             0.05            0.05     3.140185e-16 2  7.071068e-01 0.4795001
#> 4             0.05            0.05     3.140185e-16 2 -7.071068e-01 0.4795001
#> 5             0.05           -0.05     3.140185e-16 2  3.184526e+14 0.0000000
#> 6             0.05           -0.05     3.140185e-16 2  3.184526e+14 0.0000000
#> 7            -0.05           -0.05     3.140185e-16 2  7.071068e-01 0.4795001
#> 8            -0.05           -0.05     3.140185e-16 2 -7.071068e-01 0.4795001
#> 9            -0.05            0.05     9.420555e-16 2 -1.061509e+14 0.0000000
#> 10           -0.05            0.05     9.420555e-16 2 -1.061509e+14 0.0000000
#> 11            0.05            0.05     9.420555e-16 2 -7.071068e-01 0.4795001
#> 12            0.05            0.05     9.420555e-16 2  7.071068e-01 0.4795001
#> 13            0.05           -0.05     3.140185e-16 2  3.184526e+14 0.0000000
#> 14            0.05           -0.05     3.140185e-16 2  3.184526e+14 0.0000000
#> 15           -0.05           -0.05     3.140185e-16 2  7.071068e-01 0.4795001
#> 16           -0.05           -0.05     3.140185e-16 2 -7.071068e-01 0.4795001
legacy$pair_results
#>    promoter_libplate_replicate promoter libplate replicate srn_code
#> 1                    P1_lp1_r1       P1      lp1        r1       C1
#> 2                    P1_lp1_r1       P1      lp1        r1       C2
#> 9                    P2_lp1_r1       P2      lp1        r1       C1
#> 10                   P2_lp1_r1       P2      lp1        r1       C2
#>    log2FC.polished dmso.avg_dmsoFC dmso.stdv_dmsoFC n        zscore pvalue
#> 1            -0.05            0.05     3.140185e-16 2 -3.184526e+14      0
#> 2            -0.05            0.05     3.140185e-16 2 -3.184526e+14      0
#> 9            -0.05            0.05     9.420555e-16 2 -1.061509e+14      0
#> 10           -0.05            0.05     9.420555e-16 2 -1.061509e+14      0
#>    pvalue.adj           hit
#> 1           0 Downregulated
#> 2           0 Downregulated
#> 9           0 Downregulated
#> 10          0 Downregulated
```

## Compare against DStressR model-based inference

The median-polish path is intended as a reproducible baseline. The
model-based DStressR analysis starts from the same long table but uses
[`prepare_assay()`](https://bio-datascience.github.io/DStressR/reference/prepare_assay.md)
and `fit_workflow(..., workflow = "model")`:

``` r

assay <- prepare_assay(
  expression_df,
  promoter = "promoter",
  compound = "srn_code",
  control = "DMSO",
  lux = "LUX.AUC_16",
  growth = "od_16h.measured",
  batch = "batch",
  plate = "libplate",
  replicate = "replicate"
)

fit <- fit_workflow(
  assay,
  workflow = "model",
  technical = c("batch", "libplate", "replicate"),
  empirical_bayes = TRUE
)

model_results <- results(fit)
```

This gives a direct apples-to-apples setup: `workflow = "median_polish"`
recovers the legacy workflow, while `workflow = "model"` provides the
model-based replacement.
