# Campylobacter Manuscript Analysis

This is the previous downstream manuscript comparison workflow. It is kept
separate from the CRAN package and from the public Binsfeld example data.

## Input Contract

Package-generated pair-level outputs belong in:

```text
analysis/outputs/package_results/
```

Each method must provide one TSV with one row per promoter-compound pair and
the columns:

```text
promoter  compound  effect  pvalue  padj_global  padj_by_promoter
```

The current comparison scripts look for:

```text
median_polish_pair_results.tsv
destress_standard_pair_results.tsv
destress_moderated_pair_results.tsv
growth_parameters.tsv
legacy_ratio_response.tsv
modeled_response.tsv
```

For model-based DStressR outputs, `effect`, `pvalue`, and adjusted p-values
should refer to the promoter-specific effect after subtracting the compound's
global effect across promoters (`specific_*` columns from `results()`), not the
raw promoter-wise compound-vs-control total effect (`total_*`).

Set `DSTRESSR_COMPARISON_ADJUSTMENT=global` or
`DSTRESSR_COMPARISON_ADJUSTMENT=by_promoter` before running comparison scripts
to choose the adjusted p-value family used for hit calls. The default is
`global`.

## Commands

Run the manuscript figure set from the repository root:

```sh
Rscript analysis/campylobacter_manuscript/run_figures.R
```

Export the clean all-pairs table:

```sh
Rscript analysis/campylobacter_manuscript/export_default_median_pair_results.R
```

Generated PNG, PDF, matrix, and summary files are written under
`analysis/outputs/` and stay local.
