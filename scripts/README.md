# Local Model Export Scripts

This folder contains small, data-free scripts that document the exact package
calls used to create local analysis outputs. The scripts are meant to be
versioned on GitHub; the proprietary screening data and generated local output
tables are not.

## Campylobacter Default DStressR Model

`export_campy_default_model_template.R` is the public template for the
Campylobacter manuscript analysis. It expects the proprietary expression table
to be available locally through `DSTRESSR_DATA_ROOT` or the local path
conventions in `analysis/_helpers.R`.

The script fits the current default model-based analysis:

```r
fit_default <- fit_destress(
  assay,
  technical = c("libplate", "replicate"),
  empirical_bayes = TRUE,
  interaction = FALSE,
  adjustment = "global",
  background_rank = 0
)
```

It writes the promoter-specific centered result columns from `results()` to:

```text
analysis/outputs/package_results/destress_moderated_pair_results.tsv
```

This table is the input used by the downstream manuscript figure scripts. The
ordinary unmoderated model remains available as a sensitivity analysis, but the
moderated model is the default used for Campylobacter hit discovery.
