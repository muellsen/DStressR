# Binsfeld analysis with direct EVC subtraction

Run from the DStressR repository root:

```sh
Rscript analysis/ecoli_promoter_screen/run_additive_evc_analysis.R
```

This exploratory analysis writes exclusively to
`analysis/outputs/binsfeld_additive_evc/`. The current manuscript, paper figures,
and existing fits are unchanged. The two saved comparator fits must already
exist; their paths are specified in the script.

## Response and inference

We retain the existing WT filtering, reference-based growth-exponent estimation,
and log2 response construction. For each tested observation, we subtract the
mean EVC response matched by compound, dose level, and replicate:
`y_i_adjusted = y_i - mean(y_EVC_matched)`.
The EVC response uses its own estimated growth exponent. Both treatment and
Water-reference observations undergo this operation. This is subtraction on
the constructed log-response scale, not subtraction of raw luminescence or a
post-fit subtraction of total-effect estimates.

The downstream reporter-wise model retains categorical replicate and numeric
dose-level adjustment, empirical-Bayes moderation, and within-reporter BH
adjustment at 0.05. No low-rank correction is used. Growth exponents are checked
against the existing Huber fit. Direct subtraction fixes the background
coefficient at one; Huber calibration estimates a reporter-specific slope.

The package treats the matched EVC response as fixed input. Its sampling
uncertainty and the resulting dependence between observations sharing an EVC
estimate are not explicitly propagated into downstream standard errors. This
comparison evaluates the current subtraction implementation under that
assumption.

## Outputs

- `additive_evc_fit.rds`: reusable fit, separated from result extraction.
- `additive_evc_full_results.tsv`: total and reporter-relative effects, standard errors,
  p-values, adjusted p-values, and hit indicators.
- `additive_evc_significant_total_pairs.tsv`: primary significant total effects.
- `additive_evc_significant_specific_pairs.tsv`: significant reporter-relative effects.
- `three_method_pair_comparison.tsv`: aligned results for all three approaches.
- `hit_counts.tsv`, `total_hit_overlaps.tsv`, `matching_checks.tsv`: numerical
  comparison and matching audit.
- `growth_exponents.tsv`, `background_calibration.tsv`: fitted response parameters.
- `total_effect_comparison.png`: common-scale p-value histograms and volcano plots.
  Blue histogram portions indicate significant pairs; orange volcano points
  indicate significant pairs. The dashed raw-p-value guide is not the BH cutoff.
- `additive_evc_comparison.rds`, `sessionInfo.txt`: plotting inputs and session details.
