# French PFIbox integrated-activity analysis

This is the retained DStressR reanalysis of French, Coutts and Brown (2018),
*Cell Systems*, DOI <https://doi.org/10.1016/j.cels.2018.07.004>.
It reconstructs unscaled promoter activity from deposited biomass-normalized
fluorescence and uses its signed integral as the colony-level response.
Discarded scalar-summary, temporal-PCA, alternative-clustering and figure
reconstruction explorations have been removed.

## Data and numerical provenance

Download instructions, source paths and the data licence are in
[`data-raw/french_pfibox/README.md`](../../data-raw/french_pfibox/README.md).
The 1.31-GB source stays local and Git-ignored. The source-stage audit is in
`FRENCH_FIGURE3_AUDIT.md`; response and model details are in
`UNSCALED_ACTIVITY_INTEGRAL.md`. Completed cleanup and validation checks are
recorded in `WORKFLOW_CLEANUP_PLAN.md`.

## Retained pipeline

From the repository root:

```sh
# Source-structure and exact activity-reconstruction checks.
Rscript analysis/french_pfibox/audit_pfibox_data.R
Rscript analysis/french_pfibox/audit_french_figure3_inputs.R
Rscript analysis/french_pfibox/audit_activity_reconstruction.R

# Integrate unscaled activity and fit the reporter-wise DStressR model.
PFIBOX_RESPONSE_SUMMARY=unscaled_integral \
PFIBOX_TIME_WINDOW=all_available \
PFIBOX_REPORTER_LIMIT=2000 PFIBOX_PERTURBATION_LIMIT=100 \
PFIBOX_DESIGN=replicate \
Rscript analysis/french_pfibox/run_promoter_activity_analysis.R

# Re-extract ordinary tests from the saved fit; no refitting.
Rscript analysis/french_pfibox/prepare_integrated_activity_supplement.R
Rscript analysis/french_pfibox/export_public_response_tables.R

# Save the documented 10-positive and 10-negative display selection.
PFIBOX_TIMECOURSE_SELECTION=directional_volcano \
Rscript analysis/french_pfibox/plot_integrated_activity_timecourses.R
```

Set the following environment for clustering and supplementary plots:

```sh
export PFIBOX_ANALYSIS_DIR=analysis/outputs/french_pfibox/temporal_unscaled_integral_analysis
export PFIBOX_RUN_TAG=reporters2000_perturbations100_Tall_available_treatment_replicate_scalable_temporal_unscaled_integral
export PFIBOX_CLUSTER_TOP_REPORTERS=2000
export PFIBOX_CLUSTER_LINKAGE=complete
export PFIBOX_CLUSTER_DISTANCE=correlation
export PFIBOX_CLUSTER_OUTPUT_SUFFIX=total_effect_correlation_complete
export PFIBOX_ADJUSTMENT_SCOPE=global

Rscript analysis/french_pfibox/compute_paper_diagnostics.R
Rscript analysis/french_pfibox/map_french_timecourse_clusters_to_destress.R
Rscript analysis/french_pfibox/build_selected_activity_figure.R --target=total
Rscript analysis/french_pfibox/plot_pfibox_moderation_diagnostic.R --target=total

PFIBOX_DISCOVERY_TARGET=relative PFIBOX_VOLCANO_LABEL_REGION=right \
PFIBOX_FIGURE_OUTPUT_DIR="$PFIBOX_ANALYSIS_DIR/relative_effect_discovery" \
Rscript analysis/french_pfibox/plot_pfibox_model_cluster_and_discovery.R
```

The repository-wide `analysis/build_paper_figures.R` runs these stages with
the other manuscript applications. Figure-only mode reads the saved results
and trajectory data without fitting DStressR or recomputing clusters.

The attributed, compact metadata and long-format scalar tables are under
[`data-raw/french_pfibox/derived/`](../../data-raw/french_pfibox/derived/README.md).
They can be used without loading the 1.31-GB source when only model fitting is
needed. `check_integral_workflow.R` validates their keys and values; add
`--refit` to verify the entire model against the protected paper results.

## Scientific target and output

The response uses all 197 retained activity samples (100 minutes to 18 hours
on the five-minute source grid), with no per-colony min-max normalization.
Each source colony-condition contributes once. The 13 untreated source arrays
form the reference; they are not labelled as independent biological replicates.
The model adjusts for the repeated colony-position index within each reporter.

The current result has 1,797 reporters, 33 antibiotic-dose treatments,
273,792 finite observations and 59,301 total-effect tests. Global BH correction
at 0.05 yields 5,101 calls. Correlation/complete-linkage clustering includes
1,796 nonconstant reporter profiles.

The final Figure 6 is written to:

```text
analysis/outputs/french_pfibox/temporal_unscaled_integral_analysis/selected_four_figure6/
```

The four handpicked illustrations are ompC / Ampicillin-18,
msrB / Sulfamethoxazole-18, recN / Ciprofloxacin-18 and
rmf / NalidixicAcid-18. Their identities are shared across heatmap callouts,
the 2-by-2 time-series panel and volcano annotations. Selection does not alter
inference. No exploratory result is an input to this pipeline.

## Reuse the trajectory style

`activity_plotting.R` defines `plot_activity_trajectories()` for any long-format
trajectory table. Time, value, condition, observation and facet columns are
configurable. The function draws individual observations and optional pointwise
medians or means; it does not smooth trajectories, fit a model or infer effects.
The manuscript generator adds its mathematical annotations afterwards.

```r
source("analysis/french_pfibox/activity_plotting.R")
p <- plot_activity_trajectories(
  curves, time = "hours", response = "activity", condition = "condition",
  trajectory = "colony_id", panel = "pair", ncol = 2, summary = "median"
)
```
