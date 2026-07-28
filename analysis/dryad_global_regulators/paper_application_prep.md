# Dryad/Dash Fluorescence Example: Paper Integration Prep

This note prepares a second application example for the DStressR manuscript
without editing the manuscript itself.

## Role In The Paper

The Binsfeld/Campylobacter example demonstrates the primary large-scale
chemical-screen use case: many compounds, promoter reporters, matched growth
measurements, DMSO/reference normalization, hit discovery, and comparison with
legacy median-polish-style analyses.

The Dryad/Dash example should play a different role. It is smaller, but more
didactic statistically. It can show that DStressR is a general framework for
paired response/growth measurements, not a workflow restricted to Lux
screens or large antibacterial compound libraries.

The main new points are:

- fluorescence/OD response data rather than Lux/OD data;
- time-window summaries treated as reporter-like analysis units;
- response-scale estimation separated from the testing dataset;
- use of an auxiliary no-stress reporter panel to estimate growth correction;
- total-effect inference as a biologically natural endpoint;
- low-rank background subtraction as a sensitivity/decomposition analysis.

## Proposed Application Subsection

Working title:

```text
Application to a fluorescence/OD regulator-stress reporter design
```

Possible opening logic:

1. Introduce the dataset as a structured reporter-stress experiment in
   Escherichia coli global-regulator reporters.
2. State that the original measurements are GFP fluorescence and OD time
   series.
3. Explain that the weak-stress subset contains four reporters and four
   perturbations, but that a larger no-stress growth-transition panel is
   available for response-scale calibration.
4. Emphasize that this is not a discovery-scale screen; it is a compact
   validation example for the modeling framework.

## Statistical Setup

Let \(a\) denote the biological reporter/promoter and \(w\) denote a selected
time window. The analysis unit is

```text
g = (a, w)
```

with

```text
a in {Fur, MarA, SoxS, LexA}
w in {Early, Middle, Late}
```

The inferential subset therefore has 12 reporter-window units and four
weak-stress perturbations:

```text
Iron, Tetracycline, H2O2, Kanamycin
```

Standard/no-stress growth is the reference condition.

The response scale is defined from window-level summaries:

```text
R_ig = log F_ig - alpha_a log O_i
```

where \(F_{ig}\) is GFP AUC in the selected time window and \(O_i\) is matched
OD AUC in the same window. The exponent is promoter-level, not
window-specific, in the current main calibrated analysis.

Important caveat:

```text
alpha_a is estimated from window-level AUC summaries, not from a continuous
time-series model.
```

## Separation Of Response Modeling And Inference

This is the central didactic point.

Response modeling:

- uses the larger no-stress growth-transition reporter panel;
- includes 18 reporters and the same three post-ramp-up windows;
- estimates shrunken promoter-level growth exponents.

Inference:

- uses only the weak-stress 12 x 4 design;
- tests reporter-window effects under weak stress relative to standard growth;
- receives the calibrated promoter-level exponents as supplied response-scale
  parameters.

This demonstrates that DStressR can use additional control/reference data for
response modeling even when compound/stress measurements are available only
for a subset of reporters.

## Current Numerical Summary

Calibrated promoter-level growth exponents for the tested reporters:

| reporter | calibrated alpha |
|---|---:|
| Fur | -0.006 |
| MarA | 0.994 |
| SoxS | 0.251 |
| LexA | 0.644 |

Three inferential targets:

| target | tested pairs | significant pairs | median raw p-value |
|---|---:|---:|---:|
| Total effect | 48 | 36 | 2.86e-4 |
| Default specific effect | 48 | 42 | 4.58e-8 |
| Rank-adjusted total effect | 48 | 26 | 0.0175 |

Rank diagnostic:

```text
First total-effect component explains 83.2% of the total-effect matrix variance.
Observed first singular value: 4.375
Permutation 99% reference: 3.750
```

## Interpretation

Total effects are closest to the original biological question in the Dash
dataset: does a reporter change under a weak stress relative to standard
growth? These effects show the strongest agreement with the expected
regulator-stress relationships, especially Fur under excess iron.

Default specific effects ask a different question: is a reporter-window effect
large relative to the compound-wide average response? In this dataset the
specific-effect analysis is strongly influenced by a Fur-dominated structure.

Rank-adjusted total effects remove one structured background component from
the total-effect matrix. This absorbs much of the Fur-dominated axis while
retaining a total-effect interpretation for the residual signal. After this
removal, relationships including later-window MarA/tetracycline,
SoxS/hydrogen peroxide, and LexA/kanamycin become easier to inspect as
rank-adjusted total responses.

The manuscript should avoid saying that the Fur signal is an artefact. A
better phrasing is that the Fur response is a strong structured mode. Depending
on the scientific question, this mode may itself be the signal of interest or a
background component that obscures smaller reporter-specific deviations.

## Suggested Figure Plan

Main-text figure candidate:

- panel a: two-stage design schematic;
- panel b: calibrated alpha estimates for the no-stress reporter panel;
- panel c: total-effect heatmap for the 12 x 4 weak-stress analysis;
- panel d: rank-1 total-effect decomposition heatmap;
- panel e: compact p-value histogram or volcano comparison for the three
  inferential targets.

Supplemental/appendix figure candidate:

- full volcano comparison for total, default specific, and rank-adjusted total
  effects;
- rank scree/permutation diagnostic;
- table of expected regulator-stress pairs across windows.

## Existing Scripts

Primary workflows:

```sh
Rscript analysis/dryad_global_regulators/run_weak_stress_window_pipeline.R
Rscript analysis/dryad_global_regulators/run_weak_stress_window_calibrated_pipeline.R
Rscript analysis/dryad_global_regulators/run_weak_stress_window_rank1_sensitivity.R
Rscript analysis/dryad_global_regulators/plot_weak_stress_three_target_inference.R
```

Useful generated outputs:

```text
analysis/outputs/dryad_global_regulators/dryad_growth_transition_matching_windows_alpha_shrinkage_estimates.png
analysis/outputs/dryad_global_regulators/dryad_weak_stress_windows_calibrated_alpha_total_effect_heatmap.png
analysis/outputs/dryad_global_regulators/dryad_weak_stress_windows_calibrated_alpha_specific_effect_heatmap.png
analysis/outputs/dryad_global_regulators/dryad_weak_stress_windows_calibrated_alpha_rank1_effect_decomposition.png
analysis/outputs/dryad_global_regulators/dryad_weak_stress_windows_calibrated_alpha_three_target_pvalue_histograms.png
analysis/outputs/dryad_global_regulators/dryad_weak_stress_windows_calibrated_alpha_three_target_volcano.png
```

Useful tables:

```text
analysis/outputs/dryad_global_regulators/dryad_growth_transition_matching_windows_growth_exponents.tsv
analysis/outputs/dryad_global_regulators/dryad_weak_stress_windows_calibrated_alpha_pair_results.tsv
analysis/outputs/dryad_global_regulators/dryad_weak_stress_windows_calibrated_alpha_rank1_pair_results.tsv
analysis/outputs/dryad_global_regulators/dryad_weak_stress_windows_calibrated_alpha_rank_expected_pairs.tsv
analysis/outputs/dryad_global_regulators/dryad_weak_stress_windows_calibrated_alpha_rank_scree_diagnostics.tsv
analysis/outputs/dryad_global_regulators/dryad_weak_stress_windows_calibrated_alpha_three_target_summary.tsv
```

## Open Decisions For Tomorrow

- Decide whether this belongs in the main Application chapter or Appendix as a
  second application.
- Decide whether the main text should emphasize total effects, with
  specific/rank-1 effects as a sensitivity analysis.
- Decide whether the rank-1 diagnostic should be a main-text panel or
  supplemental figure.
- Decide how much of the original Dash biological expectation to discuss
  explicitly, without treating canonical pairs as ground-truth labels.
- Possibly refine the volcano plot further or replace it by a smaller
  heatmap/table if the figure becomes too dense.
