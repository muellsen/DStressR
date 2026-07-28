# Dryad E. coli Global-Regulator Reporter Dataset Audit

Dataset: Dash et al. (2024), "A library of reporters of the global regulators of
gene expression of Escherichia coli", Dryad doi:10.5061/dryad.b2rbnzsm8.

Primary article: Dash et al. (2024), mSystems, doi:10.1128/msystems.00065-24.

## Audit Status

The Dryad landing page and API were audited on 2026-07-26. The record exposes a
432 MB `Data.zip` archive and a 5.8 KB `README.md`. Direct file downloads from
the Dryad file-stream/API endpoints returned 401/403 responses in this local
shell session, but the archive was later added locally in unzipped form under
`analysis/dryad_global_regulators/raw/Data`.

## Data Structure From Public Metadata

The dataset contains several modalities:

- flow cytometry single-cell GFP and width measurements
- microscopy images
- plasmid sequences
- RNA-seq TPM summaries
- RT-PCR measurements
- spectrophotometry OD and fluorescence time series

The spectrophotometry data are the relevant part for DStressR-style modeling.
The Dryad README describes separate folders for fluorescence and OD. The
fluorescence folder contains:

- `01_Growth_Transition`
- `02_RNA_Seq_Correlation`
- `03_Weak_Stresses`

The OD folder contains:

- `01_LB_Growth_curves`
- `02_M9_Glucose_Growth_curves`
- `03_Stress_growth_curves`

The primary article states that weak-stress experiments measured time-lapse
reporter activity for Fur, MarA, SoxS, and LexA reporters under iron excess,
tetracycline, hydrogen peroxide, and kanamycin. GFP values were normalized by
the corresponding OD600 values at the same time point, and fold changes were
computed relative to standard growth at the same time point. The Figure 6
description reports SEM from three biological replicates and states that WT
cellular autofluorescence was subtracted.

## DStressR Suitability

This dataset is more structured than the Pseudomonas tailocin Biolog dataset
for formal modeling because it has matched OD and fluorescence and biological
replication for the weak-stress experiment.

It is not a direct replacement for a large compound screen:

- there are only four weak-stress perturbations in the main cross-response
  experiment
- the relevant reporter subset is four reporters, not the full 16, for the
  weak-stress specificity/sensitivity experiment
- the paper's analysis is time-resolved fold-change against standard growth,
  not a high-throughput compound-library hit-calling problem

The best DStressR use case is therefore a small, structured validation example:

- promoters/reporters: Fur, MarA, SoxS, LexA
- perturbations: iron excess, tetracycline, hydrogen peroxide, kanamycin
- reference: standard growth condition
- response signal: fluorescence, with OD as the growth/abundance signal
- candidate response summaries: time-window AUC after 50 min, endpoint/late
  mean response, or full time-point-specific modeling
- expected biological structure: near-diagonal reporter-stress specificity

This could become a valuable manuscript appendix or short application showing
that DStressR is generic beyond Lux/OD compound libraries: it can model
fluorescence/OD reporter responses in a small factorial reporter-stress design.

## Open Local-Audit Items

The weak-stress spectrophotometry workbooks were inspected and can be parsed
into a DStressR input table. The next audit should verify:

- exact workbook names and sheet structure
- whether OD and fluorescence tables share plate/well/time identifiers
- number of biological replicates per reporter-stress combination
- whether standard-growth controls are measured in the same plate/run as stress
  conditions
- whether the full 16-reporter growth-transition data can be used for response
  normalization calibration
- whether time series can be summarized into a clean long table suitable for
  `prepare_assay()`

## DStressR Pipeline

The analysis script is:

```sh
Rscript analysis/dryad_global_regulators/run_default_pipeline.R
```

It expects the unzipped Dryad archive at:

```text
analysis/dryad_global_regulators/raw/Data/
```

The script reads the weak-stress fluorescence workbook and the stress-growth OD
workbook. It subtracts the sheet-specific MG1655/WT autofluorescence trace from
reporter fluorescence, converts matched fluorescence and OD time series into
AUC summaries, and runs the moderated DStressR model. It writes outputs under:

```text
analysis/outputs/dryad_global_regulators/
```

The current run produces:

- `dryad_weak_stress_auc_input.tsv`
- `dryad_default_pair_results.tsv`
- `dryad_default_significant_pairs.tsv`
- `dryad_default_growth_exponents.tsv`
- `dryad_default_pvalue_histogram.{png,pdf}`
- `dryad_default_volcano.{png,pdf}`
- `dryad_default_effect_heatmap.{png,pdf}`
- `dryad_alpha1_pair_results.tsv`
- `dryad_alpha1_significant_pairs.tsv`
- `dryad_alpha1_pvalue_histogram.{png,pdf}`
- `dryad_alpha1_volcano.{png,pdf}`
- `dryad_alpha1_effect_heatmap.{png,pdf}`
- `dryad_estimated_alpha_vs_alpha1_comparison.tsv`

The pipeline enables EVC-style calibration only if the extracted data contain a
background reporter label such as `EVC`, `empty_vector`, or `no_promoter`.
No such reporter is present in the weak-stress spectrophotometry tables, so
the current analysis uses standard/no-stress growth as the reference condition.

Two response-scale choices are exported. `dryad_default_*` uses estimated
growth exponents to mirror the DStressR default. In this small dataset, the
no-stress controls are partly repeated across stress sheets, producing
near-perfect control fits and unstable promoter-specific exponents. Therefore
`dryad_alpha1_*` is also exported as the paper-faithful GFP/OD-style
normalization with fixed `alpha = 1`.

## Growth-Transition Library Analysis

The full reporter-library growth-transition analysis is:

```sh
Rscript analysis/dryad_global_regulators/run_growth_transition_pipeline.R
```

This analysis uses:

- `Spectrophotometry/Fluorescence/Growth_Transition/Raw_GFP_data.xlsx`
- `Spectrophotometry/OD/M9_Glucose_Growth_curves.xlsx`

It subtracts the MG1655/WT autofluorescence trace from each reporter's GFP
time series, summarizes matched GFP and OD curves into four growth-phase
windows, and treats the windows as DStressR conditions:

- `Baseline`: 0--300 min, used as the reference
- `Early transition`: 320--600 min
- `Late transition`: 620--900 min
- `Stationary`: 920--1260 min

The response uses fixed `alpha = 1`, matching the GFP/OD normalization in
Dash et al. rather than estimating promoter-specific growth exponents from
only three baseline-window replicates.

Current output dimensions:

```text
Reporters: 18
Growth-phase windows including reference: 4
Well-window summaries: 216
Tested reporter-window pairs: 54
Significant pairs at within-reporter FDR 0.05: 52
```

The main outputs are:

- `dryad_growth_transition_auc_input.tsv`
- `dryad_growth_transition_alpha1_pair_results.tsv`
- `dryad_growth_transition_alpha1_significant_pairs.tsv`
- `dryad_growth_transition_alpha1_pvalue_histogram.{png,pdf}`
- `dryad_growth_transition_alpha1_volcano.{png,pdf}`
- `dryad_growth_transition_alpha1_effect_heatmap.{png,pdf}`

## Weak-Stress Windowed Reporter Analysis

The weak-stress compound subset can also be analyzed as a time-windowed
DStressR design:

```sh
Rscript analysis/dryad_global_regulators/run_weak_stress_window_pipeline.R
```

This script uses the same weak-stress GFP and OD workbooks as the default
pipeline, but summarizes each time series in three non-overlapping windows
after the initial reporter ramp-up period:

- `Early`: 60--100 min
- `Middle`: 120--180 min
- `Late`: 200--240 min

The four original reporters are crossed with these three windows, giving 12
pseudo-reporters. The model then tests the four weak-stress conditions against
standard/no-stress growth, giving a 12 x 4 reporter-condition result matrix.
WT autofluorescence is subtracted from GFP, Standard controls measured on the
separate stress sheets are averaged within reporter-window-replicate, and OD is
used with fixed `alpha = 1`.

Current output dimensions:

```text
Base reporters: 4
Time windows: 3
Pseudo-reporters: 12
Conditions including reference: 5
Well-window summaries: 180
Tested pseudo-reporter-condition pairs: 48
Significant pairs at within-reporter FDR 0.05: 36
```

The main outputs are:

- `dryad_weak_stress_windows_alpha1_input.tsv`
- `dryad_weak_stress_windows_alpha1_pair_results.tsv`
- `dryad_weak_stress_windows_alpha1_significant_pairs.tsv`
- `dryad_weak_stress_windows_alpha1_summary.tsv`
- `dryad_weak_stress_windows_alpha1_pvalue_histogram.{png,pdf}`
- `dryad_weak_stress_windows_alpha1_volcano.{png,pdf}`
- `dryad_weak_stress_windows_alpha1_effect_heatmap.{png,pdf}`

## Weak-Stress Analysis With No-Stress Response Calibration

The calibrated response-scale analysis is:

```sh
Rscript analysis/dryad_global_regulators/run_weak_stress_window_calibrated_pipeline.R
```

This is a separated two-stage analysis. First, the larger no-stress
growth-transition reporter panel is summarized over the same post-ramp-up
windows used above. These 18 reporters provide 162 calibration rows and are
used to estimate shrunken promoter-level growth exponents. Second, the
calibrated exponents for Fur, MarA, SoxS, and LexA are supplied to the
weak-stress 12 x 4 inference model.

The estimated exponent is a property of the chosen row-level summaries. Here,
it is estimated from window-level GFP AUC and OD AUC summaries, not from a
continuous pointwise time-series model.

Current output dimensions:

```text
Calibration reporters: 18
Calibration rows: 162
Weak-stress pseudo-reporters: 12
Tested pseudo-reporter-condition pairs: 48
Significant pairs at within-reporter FDR 0.05: 42
```

The main outputs are:

- `dryad_growth_transition_matching_windows_calibration_input.tsv`
- `dryad_growth_transition_matching_windows_growth_exponents.tsv`
- `dryad_growth_transition_matching_windows_alpha_shrinkage_estimates.{png,pdf}`
- `dryad_growth_transition_matching_windows_growth_exponents_by_window.tsv`
- `dryad_growth_transition_matching_windows_alpha_heatmap.{png,pdf}`
- `dryad_growth_transition_matching_windows_alpha_curves.{png,pdf}`
- `dryad_weak_stress_windows_calibrated_alpha_supplied.tsv`
- `dryad_weak_stress_windows_calibrated_alpha_pair_results.tsv`
- `dryad_weak_stress_windows_calibrated_alpha_significant_pairs.tsv`
- `dryad_weak_stress_windows_calibrated_alpha_total_effect_heatmap.{png,pdf}`
- `dryad_weak_stress_windows_calibrated_alpha_specific_effect_heatmap.{png,pdf}`
- `dryad_weak_stress_windows_calibrated_alpha_pvalue_histogram_combined.{png,pdf}`
- `dryad_weak_stress_windows_calibrated_alpha_volcano.{png,pdf}`
- `dryad_weak_stress_windows_alpha1_vs_calibrated_alpha_comparison.tsv`

## Weak-Stress Rank-1 Sensitivity Analysis

The calibrated 12 x 4 weak-stress analysis can be rerun with one low-rank
background component removed from the total-effect matrix before testing
rank-adjusted total effects:

```sh
Rscript analysis/dryad_global_regulators/run_weak_stress_window_rank1_sensitivity.R
```

This script uses the calibrated-alpha weak-stress input and supplied
promoter-level exponents from the preceding workflow. It fits the same
moderated DStressR model twice, first with `background_rank = 0` and then with
`background_rank = 1`. The rank-1 term is estimated from the reporter-window
by stress total-effect matrix. The unadjusted total effects remain reported,
and the rank-adjusted total effect is tested as a sensitivity target.

Current output dimensions:

```text
Tested pseudo-reporter-condition pairs: 48
Significant pairs without rank removal: 42
Significant rank-adjusted total pairs after rank-1 removal: 26
First total-effect component variance explained: 83.2%
```

The main outputs are:

- `dryad_weak_stress_windows_calibrated_alpha_rank1_pair_results.tsv`
- `dryad_weak_stress_windows_calibrated_alpha_rank1_significant_pairs.tsv`
- `dryad_weak_stress_windows_calibrated_alpha_rank_sensitivity_summary.tsv`
- `dryad_weak_stress_windows_calibrated_alpha_rank_scree_diagnostics.tsv`
- `dryad_weak_stress_windows_calibrated_alpha_rank_expected_pairs.tsv`
- `dryad_weak_stress_windows_calibrated_alpha_rank1_effect_decomposition.{png,pdf}`
- `dryad_weak_stress_windows_calibrated_alpha_rank1_pvalue_histograms.{png,pdf}`
- `dryad_weak_stress_windows_calibrated_alpha_rank1_volcano.{png,pdf}`
