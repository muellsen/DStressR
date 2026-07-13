![](assets/Logo-DStressR.svg)

# DStressR: Differential stress-response modeling for chemical genomics screens

This repository hosts the `DStressR` R package and companion analysis workflow
for high-throughput bacterial chemical genomics screens. `DStressR` is designed
to go hand in hand with
[`DGrowthR`](https://bio-datascience.github.io/DGrowthR/): `DGrowthR` models
bacterial growth curves, while `DStressR` models promoter-activity responses
after accounting for growth, compound-wide effects, technical covariates, and
promoter-specific uncertainty.

![](assets/dstressr_panel.svg)

> [!TIP]
> `DStressR` is intended as the promoter-response counterpart to
> [`DGrowthR`](https://bio-datascience.github.io/DGrowthR/). The current default
> hit-determination workflow is based on the original exported luminescence and
> growth summaries; DGrowthR-derived growth parameters can be handed over
> explicitly for sensitivity analyses.

## Installation Guide

To install the `DStressR` R package directly from this repository, first clone
the repository and enter the cloned folder. Then execute the following commands
in R.

1. Ensure that you have the `devtools` package installed. If not, you can
   install it using the following command:

```r
# Install devtools
install.packages("devtools")

# Load the library
library(devtools)
```

2. Use the `install` function to install the `DStressR` package:

```r
install()
```

For permutation-based empirical p-values, install the
[`permApprox`](https://github.com/stefpeschel/permApprox) package:

```r
remotes::install_github("stefpeschel/permApprox")
```

## Get started

A worked example is available in the
[`Get started with DStressR`](https://bio-datascience.github.io/DStressR/articles/DStressR.html)
article. It walks through a complete simulated screen, from assay preparation to
model fitting, hit calling, volcano plots, and response heatmaps.

For compatibility with the original analysis, see the
[`Original median-polish workflow`](https://bio-datascience.github.io/DStressR/articles/median-polish-workflow.html)
article. Running `fit_median_polish()` on the original exported expression
table, DMSO well IDs, and noisy-DMSO exclusions reproduces the median-polish
normalization and p-value workflow used as the baseline analysis.

## How to use `DStressR`

The typical DStressR analysis starts from promoter-level luminescence and growth
summaries, together with compound, promoter, replicate, plate, and batch
metadata.

```r
library(DStressR)

assay <- prepare_assay(
  expression_df,
  promoter = "promoter",
  compound = "srn_code",
  control = "DMSO",
  lux = "LUX.AUC_16",
  growth = "od_16h.measured",
  batch = "batch",
  replicate = "replicate"
)

fit <- fit_destress(
  assay,
  technical = c("batch", "replicate"),
  empirical_bayes = TRUE
)

tab <- results(fit)
hits <- call_hits(tab, fdr = 0.05, effect = "specific_effect")
```

The fitted model separates two related quantities:

- `total_effect`: DMSO-relative promoter response for a compound.
- `specific_effect`: promoter-specific response after subtracting the
  compound-wide effect shared across promoters.

This distinction is important for compounds that globally perturb growth,
luminescence, metabolism, or assay chemistry.

## Standard plots

`DStressR` includes standard output plots for common screening summaries,
including volcano plots, promoter-compound response heatmaps, clustered
heatmaps, effect histograms, and empirical-Bayes variance diagnostics.

```r
plot_volcano(
  tab,
  effect = "specific_effect",
  padj = "specific_padj",
  top_n = 12,
  top_promoters = 6
)

plot_response_heatmap(
  tab,
  value = "specific_effect"
)
```

## Required input shape

DStressR expects a long promoter-compound table with one row per measured
promoter-compound-replicate observation. In a complete rectangular screen, the
number of rows is approximately:

```text
n_promoters x (n_compound_wells + n_control_wells) x n_replicates
```

Additional rows can occur when the same design is repeated across batches,
library plates, measurement plates, or experimental days.

At minimum, the expression table must contain columns that identify:

- promoter or reporter construct, for example `promoter`
- compound or library well, for example `compound`, `srn_code`, or
  `libplate` plus `well`
- negative-control compound, usually `DMSO`
- luminescence summary, for example `LUX.AUC_16`
- growth summary, for example `od_16h.measured`
- replicate and technical covariates, for example `replicate`, `batch`,
  `plate`, or `libplate`

These columns are mapped explicitly in `prepare_assay()`, so projects can use
their own column names:

```r
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
```

## Campylobacter workflow input files

For the original Campylobacter promoter-library workflow, DStressR includes the
helper `read_campylobacter_expression()`. It joins two exported files:

1. `expression_values.tsv.gz`

   A long measurement table with one row per promoter-library-well-replicate
   observation. It should contain either:

   - `srn_code`, a unique library-well identifier, or
   - both `libplate` and `well`, from which `srn_code` is reconstructed as
     `paste(libplate, well, sep = "_")`.

   It should also contain the promoter identifier, luminescence summary, growth
   summary, and technical covariates used downstream in `prepare_assay()`.

2. `LibMap.txt`

   A library annotation table with one row per compound/library well. Required
   columns are:

   - `Library plate`
   - `Well`
   - `ProductName`
   - `Catalog Number`

   The helper converts these to a compound key
   `srn_code = paste0("lp", Library plate, "_", Well)` and joins
   `ProductName` and `Catalog Number` onto the expression table.

```r
expression_df <- read_campylobacter_expression(
  expression_file = "expression_values.tsv.gz",
  libmap_file = "LibMap.txt"
)
```

The returned object has the same number of rows as `expression_values.tsv.gz`,
with compound annotations added. This joined table can then be passed directly
to `prepare_assay()`.

To reproduce the original median-polish workflow, provide the DMSO library-well
IDs and optional noisy-DMSO well IDs from `LibMap.txt`:

```r
libmap <- read.delim("LibMap.txt", check.names = FALSE)
libmap$srn_code <- paste0("lp", libmap[["Library plate"]], "_", libmap[["Well"]])

dmso_srn_codes <- libmap$srn_code[libmap$ProductName == "DMSO"]
dmso_noisy_srn_codes <- libmap$srn_code[libmap$ProductName == "DMSO noisy"]

legacy <- fit_median_polish(
  expression_df,
  response = "log2.auc.16hmeasured.normed",
  control = dmso_srn_codes,
  exclude = dmso_noisy_srn_codes
)

replicate_pvalues <- legacy$replicate_results
hit_table <- legacy$pair_results
```

## Optional DGrowthR handoff

If growth curves have already been modeled with DGrowthR, DStressR can use a
chosen DGrowthR growth parameter as the growth column for a sensitivity
analysis:

```r
expression_df2 <- add_dgrowthr_growth(
  expression_df,
  object = dgrowthr_fit,
  by = "curve_id",
  model_covariate = "curve_id",
  growth_metric = "OD_16",
  output = "dgrowthr_od16"
)

assay <- prepare_assay(
  expression_df2,
  promoter = "promoter",
  compound = "srn_code",
  control = "DMSO",
  lux = "LUX.AUC_16",
  growth = "dgrowthr_od16"
)
```

## Analysis workflow

The `analysis/` folder contains scripts used during model development and
benchmarking against the original median-polish workflow, including p-value
comparisons, empirical-Bayes diagnostics, response matrices, clustered
heatmaps, network summaries, and empirical replicate/permutation tests.
