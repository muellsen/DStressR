# Kochanowski Promoter-Scaling Analysis

This folder contains an exploratory DStressR summary-level analysis of
Kochanowski et al. (2017), *Molecular Systems Biology*. The released
supplementary tables provide processed steady-state promoter activities and
estimated standard deviations, not replicate-level wells. The analysis therefore
uses `fit_destress_summary()` with supplied response variances.

The script expects the Kochanowski supplementary Excel files under
`data-raw/kochanowski2017/source/`. These source files are intentionally kept
out of Git.

```sh
Rscript analysis/kochanowski_promoter_scaling/run_kochanowski_destress_analysis.R
Rscript analysis/kochanowski_promoter_scaling/plot_paper_figures.R
```

Generated tables and figures are written to
`analysis/outputs/kochanowski_promoter_scaling/`, which is also ignored by Git.
The analysis reports a rank-1 low-rank component for comparison with the first
singular-vector analysis in Kochanowski et al.; the automatic rank diagnostic is
saved separately.

The first command saves the fit objects, standardized results, parameters,
diagnostics, and published annotations in `kochanowski_paper.rds`. The second
only renders figures from that artifact. The shared paper rebuild interface
is documented in `analysis/PAPER_FIGURES.md`.
