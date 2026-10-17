#!/usr/bin/env Rscript
# Run from the repository root. Figure-only mode never fits a statistical model.
source(file.path("analysis", "_helpers.R"))
mode <- commandArgs(trailingOnly = TRUE)
if (!length(mode)) mode <- "--figures"
stopifnot(length(mode) == 1L)
mode <- match.arg(mode, c("--analysis", "--diagnostics", "--figures", "--all"))
run <- function(script, env = character(), args = character()) {
  message("\n--- ", script, " ---")
  status <- system2(file.path(R.home("bin"), "Rscript"), c(shQuote(script), args), env = env)
  if (status != 0L) stop("Failed: ", script, call. = FALSE)
}
pfibox <- c("PFIBOX_REPORTER_LIMIT=2000", "PFIBOX_PERTURBATION_LIMIT=100",
             "PFIBOX_RESPONSE_SUMMARY=unscaled_integral", "PFIBOX_TIME_WINDOW=all_available", "PFIBOX_DESIGN=replicate",
             "PFIBOX_PERTURBATION_TARGET=treatment", "PFIBOX_MODEL_PATH=scalable")
cluster <- c("PFIBOX_RUN_TAG=reporters2000_perturbations100_Tall_available_treatment_replicate_scalable_temporal_unscaled_integral",
              "PFIBOX_ANALYSIS_DIR=analysis/outputs/french_pfibox/temporal_unscaled_integral_analysis",
              "PFIBOX_ADJUSTMENT_SCOPE=global", "PFIBOX_CLUSTER_TOP_REPORTERS=2000",
              "PFIBOX_CLUSTER_LINKAGE=complete", "PFIBOX_CLUSTER_DISTANCE=correlation",
              "PFIBOX_CLUSTER_OUTPUT_SUFFIX=total_effect_correlation_complete")
if (mode %in% c("--analysis", "--all")) {
  run("analysis/ecoli_promoter_screen/build_workflow_comparison_report.R")
  run("analysis/ecoli_promoter_screen/run_evc_calibrated_analysis.R")
  run("analysis/ecoli_promoter_screen/run_interaction_model.R")
  run("analysis/ecoli_promoter_screen/prepare_response_figure_data.R")
  run("analysis/french_pfibox/run_promoter_activity_analysis.R",
      c(pfibox, "PFIBOX_EMPIRICAL_BAYES=true"))
  run("analysis/french_pfibox/prepare_integrated_activity_supplement.R")
  run("analysis/french_pfibox/export_public_response_tables.R")
  run("analysis/french_pfibox/plot_integrated_activity_timecourses.R",
      "PFIBOX_TIMECOURSE_SELECTION=directional_volcano")
  run("analysis/kochanowski_promoter_scaling/run_kochanowski_destress_analysis.R")
  run("analysis/dryad_global_regulators/run_weak_stress_window_pipeline.R")
  run("analysis/dryad_global_regulators/run_weak_stress_window_calibrated_pipeline.R")
  run("analysis/dryad_global_regulators/run_weak_stress_window_rank1_adjustment.R")
}
if (mode %in% c("--diagnostics", "--all")) {
  run("analysis/ecoli_promoter_screen/compute_paper_diagnostics.R")
  run("analysis/ecoli_promoter_screen/prepare_total_effect_discovery.R")
  run("analysis/ecoli_promoter_screen/prepare_interaction_total_effects.R")
  run("analysis/french_pfibox/compute_paper_diagnostics.R", cluster)
  run("analysis/french_pfibox/map_french_timecourse_clusters_to_destress.R", cluster)
}
if (mode %in% c("--figures", "--all")) {
  run("analysis/check_paper_figure_contract.R")
  scripts <- c(
    "analysis/ecoli_promoter_screen/plot_response_modeling.R",
    "analysis/ecoli_promoter_screen/plot_binsfeld_effect_matrix_views.R",
    "analysis/ecoli_promoter_screen/plot_binsfeld_effect_matrix_and_variance_diagnostic.R",
    "analysis/ecoli_promoter_screen/plot_discovery_figures.R",
    "analysis/ecoli_promoter_screen/plot_interaction_figures.R",
    "analysis/kochanowski_promoter_scaling/plot_paper_figures.R",
    "analysis/dryad_global_regulators/plot_calibration_figures.R",
    "analysis/dryad_global_regulators/plot_rank_adjustment_figures.R",
    "analysis/dryad_global_regulators/plot_weak_stress_three_target_inference.R")
  for (script in scripts) run(script, cluster)
  run("analysis/ecoli_promoter_screen/plot_discovery_figures.R", args = "--target=total")
  run("analysis/french_pfibox/build_selected_activity_figure.R", args = "--target=total")
  run("analysis/french_pfibox/plot_pfibox_moderation_diagnostic.R", cluster, args = "--target=total")
  run("analysis/french_pfibox/plot_pfibox_model_cluster_and_discovery.R",
      c(cluster, "PFIBOX_DISCOVERY_TARGET=relative", "PFIBOX_VOLCANO_LABEL_REGION=right",
        "PFIBOX_FIGURE_OUTPUT_DIR=analysis/outputs/french_pfibox/temporal_unscaled_integral_analysis/relative_effect_discovery"))
}
