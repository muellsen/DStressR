# DStressR Analysis Scripts

This folder is a downstream comparison layer. It must not implement
estimators, p-value calculations, variance moderation, replicate aggregation,
or multiple-testing correction.

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

These files should be exported by package code, for example from
`fit_workflow(..., workflow = "median_polish")$pair_results` or
`results(fit_workflow(..., workflow = "model", ...))`, after any selected
package-level p-value adjustment.

For model-based DStressR outputs, `effect`, `pvalue`, and adjusted p-values
should refer to the promoter-specific effect after subtracting the compound's
global effect across promoters (`specific_*` columns from `results()`), not the
raw promoter-wise compound-vs-control total effect (`total_*`).

`padj_global` is the primary DStressR-style full-family Benjamini-Hochberg
adjustment across all promoter-compound pairs. `padj_by_promoter` is retained
for legacy median-polish comparability.

Set `DSTRESSR_COMPARISON_ADJUSTMENT=global` or
`DSTRESSR_COMPARISON_ADJUSTMENT=by_promoter` before running comparison scripts
to choose the adjusted p-value family used for hit calls. The default is
`global`.

## Scripts

- `run_campy_manuscript_figures.R`: regenerates the Campylobacter manuscript
  figure set from local package-generated outputs.
- `plot_growth_parameter_estimates.R`: plots the response-scale parameters
  returned by the package (`a_raw`, `alpha_raw`, and `alpha_shrunk`).
- `plot_response_heatmap_comparison.R`: plots the fixed-ratio response,
  modeled response, and modeled-minus-fixed-ratio difference heatmaps.
- `plot_significance_summary.R`: creates the DStressR standard-model and
  median-polish max-p p-value/table panels, plus the volcano/Venn summary.
- `plot_moderated_hit_bipartite_heatmap.R`: plots the DStressR standard-model
  hit/effect heatmap against the median-polish max-p calls in the same
  compound ordering. The file name is historical; the method used by the
  current script is documented in the plot title and output tables.
- `plot_moderated_hit_network_no_pcmeA_unique.R`: plots the DStressR standard
  hit network after removing compounds whose significant hit is unique to
  `PcmeA`. The historical manuscript output names are still emitted as
  compatibility aliases.
- `compare_pair_level_pvalues.R`: merges package pair-level outputs and writes
  comparison summaries, hit membership, and p-value scatter plots.
- `plot_bh_hit_venn.R`: plots overlap of package-level hit calls.
- `export_rejected_pair_list.R`: exports the union of differential
  promoter-compound pairs with method-specific values.
- `plot_bh_hit_network.R`: draws a network view of the differential pairs.

If a script needs a statistic that is not already in the package output, add it
to the package output first rather than computing it here.

## Campylobacter Manuscript Figures

Run from the repository root:

```sh
Rscript analysis/run_campy_manuscript_figures.R
```

The command assumes that the proprietary package outputs listed above are
already present under `analysis/outputs/package_results/`. Generated PNG, PDF,
matrix, and summary files are written under `analysis/outputs/comparisons/`.
That directory is ignored by Git and should stay local.

Current manuscript figure mapping:

| Manuscript panel | Analysis script | Primary output |
| --- | --- | --- |
| Growth parameters | `plot_growth_parameter_estimates.R` | `comparisons/growth_parameters/growth_parameter_estimates.png` |
| Modeled response heatmap | `plot_response_heatmap_comparison.R` | `comparisons/response_heatmaps/modeled_response_heatmap.png` |
| Fixed-ratio response heatmap | `plot_response_heatmap_comparison.R` | `comparisons/response_heatmaps/legacy_ratio_response_heatmap.png` |
| Modeled-minus-fixed-ratio heatmap | `plot_response_heatmap_comparison.R` | `comparisons/response_heatmaps/modeled_minus_ratio_response_heatmap.png` |
| DStressR p-value/table panel | `plot_significance_summary.R` | `comparisons/significance_summary/dstressr_standard_model_pvalue_histograms_with_hit_table.png` |
| Median-polish p-value/table panel | `plot_significance_summary.R` | `comparisons/significance_summary/median_polish_model_pvalue_histograms_with_hit_table.png` |
| Volcano/Venn summary | `plot_significance_summary.R` | `comparisons/significance_summary/volcano_venn_significance_summary.png` |
| Hit/effect heatmap | `plot_moderated_hit_bipartite_heatmap.R` | `comparisons/hit_bipartite_heatmap/dstressr_standard_order_hit_bipartite_heatmap.png` |
| Hit network | `plot_moderated_hit_network_no_pcmeA_unique.R` | `comparisons/hit_network/dstressr_standard_network_no_pcmeA_unique.png` |
