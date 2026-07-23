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

Run the primary three-set comparison used in the manuscript:

```sh
Rscript analysis/ecoli_promoter_screen/run_evc_calibrated_analysis.R
```

This workflow compares:

- the reconstructed reference Wilcoxon/Z-score rule from the original screen;
- DStressR with the default modeled response;
- DStressR with the default modeled response plus EVC-Huber calibration.

Current WT result:

```text
Reference WT hits: 53
DStressR default modeled-response WT hits: 80
DStressR EVC-Huber WT hits: 92
All three analyses: 35
Reference-only hits: 16
DStressR modeled-only hits: 9
DStressR EVC-Huber-only hits: 19
Modeled + EVC-Huber only: 36
Union significant by at least one primary analysis: 117
```

## Supplemental Comparison Scripts

Rebuild the two-method comparison between the reconstructed reference rule and
the default DStressR modeled-response analysis:

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
