# E. coli promoter-compound screen

This article describes the public *E. coli* reporter screen used as the
package application for DStressR. The data are from Binsfeld et
al. (2025), PLOS Biology, and are distributed with the package as
compact, documented R data objects. The AUC-level screen is available as
`binsfeld_reporter_auc`; the matching author score and Z-score table is
available as `binsfeld_reporter_scores`.

The source files are the PLOS Biology supplementary data and the
associated Zenodo archive:

- article: <https://doi.org/10.1371/journal.pbio.3003260>
- archive: <https://doi.org/10.5281/zenodo.15600688>

## Public data objects

``` r

binsfeld_reporter_auc <- DStressR::binsfeld_reporter_auc
binsfeld_reporter_scores <- DStressR::binsfeld_reporter_scores

dim(binsfeld_reporter_auc)
#> [1] 24576    13
table(binsfeld_reporter_auc$strain)
#> 
#> DmarA  Drob DsoxS    WT 
#>  6144  6144  6144  6144
length(unique(binsfeld_reporter_auc$promoter))
#> [1] 8
length(unique(binsfeld_reporter_auc$compound))
#> [1] 95
```

The AUC table is a long table with one row per strain, promoter,
compound or control, replicate, and serial-dilution observation. The
package-level `compound` column keeps the original compound labels,
except that the water control wells are collapsed to `Water`. The
original label remains in `drug`. The source `concentration_index`
increases with dilution; `dose_level` reverses this coding so that
larger values correspond to higher concentration.

``` r

names(binsfeld_reporter_auc)
#>  [1] "strain"              "promoter"            "replicate"          
#>  [4] "well"                "drug"                "compound"           
#>  [7] "concentration_index" "dose_level"          "concentration_ug_ml"
#> [10] "od_auc"              "lux_auc"             "od_auc_per_lux_auc" 
#> [13] "removed"
```

The score table stores the public author scores and Z-scores in long
form. It is used by the reproducibility workflow to reconstruct the
reference hit rule.

``` r

names(binsfeld_reporter_scores)
#> [1] "well"                "drug"                "compound"           
#> [4] "concentration_ug_ml" "strain"              "statistic"          
#> [7] "promoter"            "replicate"           "value"
```

## DStressR analysis

For the package application, the wild-type reporter data are filtered to
rows that passed the original quality-control flag. Water is used as the
reference condition. The Empty Vector Control reporter is supplied to
[`prepare_assay()`](https://muellsen.github.io/DStressR/reference/prepare_assay.md)
as a matched background reporter, so the default response model performs
Huber-calibrated background adjustment before model-based testing.

``` r

wt_auc <- subset(
  binsfeld_reporter_auc,
  strain == "WT" & removed == "No"
)

assay <- prepare_assay(
  wt_auc,
  promoter = "promoter",
  compound = "compound",
  control = "Water",
  lux = "lux_auc",
  growth = "od_auc",
  growth_exponent = "estimate",
  batch = "dose_level",
  replicate = "replicate",
  background_promoter = "EVC",
  background_by = c("compound", "dose_level", "replicate")
)

fit <- fit_destress(
  assay,
  preset = "model",
  technical = c("replicate", "dose_level"),
  empirical_bayes = TRUE,
  adjustment = "by_promoter",
  interaction = FALSE
)

hits <- call_hits(
  results(fit),
  fdr = 0.05,
  effect = "specific_effect"
)
```

The full comparison workflow is kept outside the package vignette build
because it creates manuscript-scale figures and tables. From the
repository root, run:

``` sh
Rscript analysis/ecoli_promoter_screen/run_evc_calibrated_analysis.R
```

The workflow compares three hit sets:

- the reconstructed Binsfeld-style Wilcoxon/Z-score reference rule
- DStressR with the default modeled response
- DStressR with the default modeled response and EVC-Huber calibration

The current wild-type comparison gives:

``` text
Reference WT hits: 53
DStressR default modeled-response WT hits: 77
DStressR EVC-Huber WT hits: 92
All three analyses: 35
Reference-only hits: 16
DStressR modeled-only hits: 8
DStressR EVC-Huber-only hits: 21
Modeled + EVC-Huber only: 34
Union significant by at least one primary analysis: 116
```

Additional scripts in `analysis/ecoli_promoter_screen/` recreate the
two-method report, volcano plots, heatmaps, Venn diagrams, and
response-modeling figures used during manuscript development. Generated
outputs are written under `analysis/outputs/` and are intentionally
ignored by Git.
