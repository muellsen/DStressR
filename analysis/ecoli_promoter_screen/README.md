# E. coli Promoter-Compound Screen Analysis

This folder contains the public *E. coli* promoter-compound screen analyses
used for the DStressR package application. The data are derived from Binsfeld
et al. (2025), PLOS Biology, and the installed DStressR package ships compact
documented data objects prepared from the public supplements.

Run scripts from the repository root.

## Paper Application Workflow

Generate response-construction figures for the manuscript:

```sh
Rscript analysis/ecoli_promoter_screen/plot_modeling_step_figures.R
```

Generate compact DStressR diagnostic figures, including the DESeq-style
compound mean-variance diagnostic and the significant-hit heatmap:

```sh
Rscript analysis/ecoli_promoter_screen/plot_binsfeld_compound_mean_variance_diagnostic.R
Rscript analysis/ecoli_promoter_screen/plot_promoter_effect_distributions.R
Rscript analysis/ecoli_promoter_screen/plot_binsfeld_hit_heatmap.R
```

Run the primary three-set comparison used in the manuscript:

```sh
Rscript analysis/ecoli_promoter_screen/run_evc_calibrated_analysis.R
```

This workflow compares:

- the reconstructed reference Wilcoxon/Z-score rule from the original screen;
- DStressR without EVC;
- DStressR with EVC.

Current WT result:

```text
Reference WT hits: 53
DStressR without EVC WT hits: 80
DStressR with EVC WT hits: 92
All three analyses: 35
Reference-only hits: 16
DStressR-without-EVC only hits: 9
DStressR-with-EVC only hits: 19
Both DStressR workflows only: 36
Union significant by at least one primary analysis: 117
```

Run the full promoter-by-compound interaction-model analysis for Appendix E:

```sh
Rscript analysis/ecoli_promoter_screen/run_interaction_sensitivity.R
```

Current WT interaction sensitivity:

```text
Interaction with EVC hits: 56
Interaction without EVC hits: 51
DStressR with EVC hits recovered by interaction with EVC: 56
DStressR without EVC hits recovered by interaction without EVC: 50
All five methods overlap: 28
```

## Supplemental Comparison Scripts

Rebuild the two-method comparison between the reconstructed reference rule and
DStressR without EVC:

```sh
Rscript analysis/ecoli_promoter_screen/compare_reporter_hits.R
Rscript analysis/ecoli_promoter_screen/plot_comparison_figures.R
Rscript analysis/ecoli_promoter_screen/build_shareable_report.R
```

Build the older three-method sensitivity report that also includes fixed
alpha=1 response normalization:

```sh
Rscript analysis/ecoli_promoter_screen/build_three_method_report.R
```

The fixed alpha=1 report is retained as a sensitivity analysis and for
continuity with earlier local reports; it is not the primary manuscript
comparison.

Generated figures, tables, and HTML reports are written under
`analysis/outputs/` and remain ignored by Git.
