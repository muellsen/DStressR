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
promoter  compound  effect  pvalue  padj
```

The current comparison scripts look for:

```text
median_polish_pair_results.tsv
destress_standard_pair_results.tsv
destress_moderated_pair_results.tsv
```

These files should be exported by package code, for example from
`fit_workflow(..., workflow = "median_polish")$pair_results` or
`results(fit_workflow(..., workflow = "model", ...))`, after any selected
package-level p-value adjustment.

## Scripts

- `compare_pair_level_pvalues.R`: merges package pair-level outputs and writes
  comparison summaries, hit membership, and p-value scatter plots.
- `plot_bh_hit_venn.R`: plots overlap of package-level hit calls.
- `export_rejected_pair_list.R`: exports the union of differential
  promoter-compound pairs with method-specific values.
- `plot_bh_hit_network.R`: draws a network view of the differential pairs.

If a script needs a statistic that is not already in the package output, add it
to the package output first rather than computing it here.
