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
  value = "specific_effect",
  padj = "specific_padj"
)
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

For the current Campylobacter workflow, `read_campylobacter_expression()` reads
the exported `expression_values.tsv.gz` and `LibMap.txt` files into the expected
shape.

## Analysis workflow

The `analysis/` folder contains scripts used during model development and
benchmarking against the original median-polish workflow, including p-value
comparisons, empirical-Bayes diagnostics, response matrices, clustered
heatmaps, network summaries, and empirical replicate/permutation tests.
