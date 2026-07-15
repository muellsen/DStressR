# Salmonella Empty Vector workflow

The Salmonella StressRegNet workflow uses Empty Vector Control (EVC)
reporter strains to estimate compound-specific background Lux signal. In
the original analysis, the `PEVC3` reporter is used as the reference.
For each compound, the two PEVC3 replicate values are averaged, and this
compound-specific EVC average is subtracted from every reporter strain.

This is different from the Campylobacter median-polish baseline. The
Salmonella baseline is:

``` text
growth-normalized Lux
  -> promoter-replicate median centering
  -> subtract compound-specific PEVC3 average
  -> estimate DMSO null per promoter-replicate
  -> z-test p-values
  -> retain largest replicate p-value per promoter-compound
  -> BH adjustment within promoter
```

DStressR implements this compatibility path with
`fit_workflow(..., workflow = "empty_vector_control")`. The lower-level
[`fit_empty_vector_control()`](https://muellsen.github.io/DStressR/reference/fit_empty_vector_control.md)
function remains available for existing scripts.

## Input files

The workflow starts from the processed Salmonella Lux-estimation output:

### `lux_auc_filtered_median.tsv.gz`

This is a long table with one row per promoter-library-well-replicate
observation. Its rows are approximately:

``` text
n_promoters x n_library_wells x n_replicates
```

after growth filtering and reporter-quality filtering. The required
columns are:

- `promoter`
- `srn_code`
- `replicate`
- `log2.lux.normed.centered`

The response column `log2.lux.normed.centered` is already:

1.  Lux AUC until 10 hours,
2.  divided by the OD reached at 10 hours,
3.  log2-transformed,
4.  median-centered within promoter-replicate.

### `LibMap.tsv.gz`

This file annotates library wells and identifies DMSO controls. Required
columns are:

- `Library plate`
- `New well`
- `Catalog Number`
- `ProductName`

The compound key used by the workflow is:

``` r

srn_code = paste(sub("LibPlate", "lp", `Library plate`), `New well`, sep = "_")
```

## Reproduce the original Salmonella workflow

``` r

expression_df <- read.delim(gzfile("lux_auc_filtered_median.tsv.gz"),
                            check.names = FALSE)
libmap <- read.delim(gzfile("LibMap.tsv.gz"), check.names = FALSE)

libmap$libplate <- sub("LibPlate", "lp", libmap[["Library plate"]])
libmap$srn_code <- paste(libmap$libplate, libmap[["New well"]], sep = "_")

dmso_srn_codes <- libmap$srn_code[libmap[["Catalog Number"]] == "DMSO"]
dmso_noisy_srn_codes <- libmap$srn_code[libmap[["Catalog Number"]] == "DMSO noisy"]

evc <- fit_workflow(
  expression_df,
  workflow = "empty_vector_control",
  promoter = "promoter",
  compound = "srn_code",
  replicate = "replicate",
  response = "log2.lux.normed.centered",
  empty_vector_promoter = "PEVC3",
  control = dmso_srn_codes,
  exclude = dmso_noisy_srn_codes,
  remove_promoters = "PmgrR"
)
```

The replicate-level table contains the EVC-normalized fold change, DMSO
null parameters, z-scores, and raw p-values:

``` r

replicate_pvalues <- evc$replicate_results

head(replicate_pvalues)
```

The pair-level table reproduces the original conservative aggregation:
for each promoter-compound pair, retain the largest replicate-level
p-value and apply BH adjustment within promoter.

``` r

hit_table <- evc$pair_results

head(hit_table)
table(hit_table$hit)
```

## What the function does

`fit_workflow(..., workflow = "empty_vector_control")` follows the
original Salmonella hit-determination logic:

1.  Remove noisy DMSO wells and failed reporters, for example `PmgrR`.
2.  Keep the already centered response column, usually
    `log2.lux.normed.centered`.
3.  For each `srn_code`, compute the average PEVC3 response across PEVC3
    replicates.
4.  Compute
    `log.evcfc = expression_value - average_PEVC3_expression_for_srn_code`.
5.  Remove the PEVC3 reference rows from the tested reporter table.
6.  Gather DMSO `log.evcfc` values within each promoter-replicate group.
7.  Estimate `dmso.mean` and `dmso.stdv`.
8.  Compute `zscore = (log.evcfc - dmso.mean) / dmso.stdv`.
9.  Compute two-sided Gaussian p-values.
10. Retain the largest replicate-level p-value for each
    promoter-compound pair.
11. Apply BH correction within promoter and label hits.

## Minimal executable example

``` r

toy <- expand.grid(
  promoter = c("PEVC3", "P1", "P2"),
  replicate = c("r1", "r2"),
  srn_code = c("DMSO1", "DMSO2", "C1"),
  stringsAsFactors = FALSE
)

toy$value <- NA_real_
toy$value[toy$promoter == "PEVC3" & toy$srn_code == "DMSO1"] <- c(1.0, 1.2)
toy$value[toy$promoter == "PEVC3" & toy$srn_code == "DMSO2"] <- c(1.1, 1.3)
toy$value[toy$promoter == "PEVC3" & toy$srn_code == "C1"] <- c(2.0, 2.2)

toy$value[toy$promoter == "P1" & toy$srn_code == "DMSO1"] <- c(1.5, 1.7)
toy$value[toy$promoter == "P1" & toy$srn_code == "DMSO2"] <- c(1.6, 1.8)
toy$value[toy$promoter == "P1" & toy$srn_code == "C1"] <- c(4.5, 4.7)

toy$value[toy$promoter == "P2" & toy$srn_code == "DMSO1"] <- c(0.8, 1.0)
toy$value[toy$promoter == "P2" & toy$srn_code == "DMSO2"] <- c(0.9, 1.1)
toy$value[toy$promoter == "P2" & toy$srn_code == "C1"] <- c(1.5, 1.7)

evc <- fit_workflow(
  toy,
  workflow = "empty_vector_control",
  response = "value",
  control = c("DMSO1", "DMSO2")
)

evc$replicate_results
#>    promoter_replicate promoter replicate srn_code log.evcfc empty_vector_mean
#> 1               P1_r1       P1        r1    DMSO1       0.4               1.1
#> 2               P1_r1       P1        r1    DMSO2       0.4               1.2
#> 3               P1_r1       P1        r1       C1       2.4               2.1
#> 4               P2_r1       P2        r1    DMSO1      -0.3               1.1
#> 5               P2_r1       P2        r1    DMSO2      -0.3               1.2
#> 6               P2_r1       P2        r1       C1      -0.6               2.1
#> 7               P1_r2       P1        r2    DMSO1       0.6               1.1
#> 8               P1_r2       P1        r2    DMSO2       0.6               1.2
#> 9               P1_r2       P1        r2       C1       2.6               2.1
#> 10              P2_r2       P2        r2    DMSO1      -0.1               1.1
#> 11              P2_r2       P2        r2    DMSO2      -0.1               1.2
#> 12              P2_r2       P2        r2       C1      -0.4               2.1
#>    empty_vector_n dmso.mean    dmso.stdv n        zscore    pvalue
#> 1               2       0.4 0.000000e+00 2           NaN       NaN
#> 2               2       0.4 0.000000e+00 2           NaN       NaN
#> 3               2       0.4 0.000000e+00 2           Inf 0.0000000
#> 4               2      -0.3 7.850462e-17 2  7.071068e-01 0.4795001
#> 5               2      -0.3 7.850462e-17 2 -7.071068e-01 0.4795001
#> 6               2      -0.3 7.850462e-17 2 -3.821431e+15 0.0000000
#> 7               2       0.6 0.000000e+00 2           NaN       NaN
#> 8               2       0.6 0.000000e+00 2           NaN       NaN
#> 9               2       0.6 0.000000e+00 2           Inf 0.0000000
#> 10              2      -0.1 0.000000e+00 2           NaN       NaN
#> 11              2      -0.1 0.000000e+00 2           NaN       NaN
#> 12              2      -0.1 0.000000e+00 2          -Inf 0.0000000
evc$pair_results
#>   promoter_replicate promoter replicate srn_code log.evcfc empty_vector_mean
#> 3              P1_r1       P1        r1       C1       2.4               2.1
#> 6              P2_r1       P2        r1       C1      -0.6               2.1
#> 4              P2_r1       P2        r1    DMSO1      -0.3               1.1
#> 5              P2_r1       P2        r1    DMSO2      -0.3               1.2
#>   empty_vector_n dmso.mean    dmso.stdv n        zscore    pvalue pvalue.adj
#> 3              2       0.4 0.000000e+00 2           Inf 0.0000000  0.0000000
#> 6              2      -0.3 7.850462e-17 2 -3.821431e+15 0.0000000  0.0000000
#> 4              2      -0.3 7.850462e-17 2  7.071068e-01 0.4795001  0.4795001
#> 5              2      -0.3 7.850462e-17 2 -7.071068e-01 0.4795001  0.4795001
#>             hit
#> 3   Upregulated
#> 6 Downregulated
#> 4        Not DE
#> 5        Not DE
```

## Local workflow equivalence check

Using the local Salmonella workflow files,
`fit_workflow(..., workflow = "empty_vector_control")` recovers the
original `04-hit_determination/hit_table.tsv.gz` output numerically:

``` text
pair rows: 68238
max absolute delta in log.evcfc: 4.44e-16
max absolute delta in zscore:    3.55e-15
max absolute delta in pvalue:    1.22e-15
max absolute delta in BH pvalue: 2.11e-15
hit label mismatches: 0
```

Those values are floating-point roundoff, so this is effectively a
one-to-one reproduction of the original Salmonella workflow.
