#!/usr/bin/env Rscript
# Run from the repository root. Figure-only mode never fits a statistical model.
source(file.path("analysis", "_helpers.R"))
mode <- commandArgs(trailingOnly = TRUE)
if (!length(mode)) mode <- "--figures"
stopifnot(length(mode) == 1L)
mode <- match.arg(mode, c("--analysis", "--diagnostics", "--figures", "--all"))
run <- function(script, env = character()) {
  message("\n--- ", script, " ---")
  status <- system2(file.path(R.home("bin"), "Rscript"), shQuote(script), env = env)
  if (status != 0L) stop("Failed: ", script, call. = FALSE)
}
pfibox <- c("PFIBOX_REPORTER_LIMIT=2000", "PFIBOX_PERTURBATION_LIMIT=100",
             "PFIBOX_T_START_H=2", "PFIBOX_T_END_H=18", "PFIBOX_DESIGN=replicate",
             "PFIBOX_PERTURBATION_TARGET=treatment", "PFIBOX_MODEL_PATH=scalable")
cluster <- c("PFIBOX_RUN_TAG=reporters2000_perturbations100_T2-18h_treatment_replicate_scalable",
              "PFIBOX_CLUSTER_TOP_REPORTERS=2000", "PFIBOX_CLUSTER_LINKAGE=ward.D2",
              "PFIBOX_CLUSTER_DISTANCE=euclidean", "PFIBOX_CLUSTER_OUTPUT_SUFFIX=total_effect_ward")
if (mode %in% c("--analysis", "--all")) {
  run("analysis/ecoli_promoter_screen/build_workflow_comparison_report.R")
  run("analysis/ecoli_promoter_screen/run_evc_calibrated_analysis.R")
  run("analysis/ecoli_promoter_screen/run_interaction_model.R")
  run("analysis/ecoli_promoter_screen/prepare_response_figure_data.R")
  run("analysis/french_pfibox/run_promoter_activity_analysis.R",
      c(pfibox, "PFIBOX_EMPIRICAL_BAYES=true"))
  response <- file.path("analysis", "outputs", "french_pfibox", "promoter_activity_analysis",
    "pfibox_pa_response_table_reporters2000_perturbations100_T2-18h_treatment_replicate_scalable.tsv")
  run("analysis/french_pfibox/run_promoter_activity_analysis.R",
      c(pfibox, "PFIBOX_EMPIRICAL_BAYES=false", paste0("PFIBOX_RESPONSE_TABLE=", response)))
  run("analysis/kochanowski_promoter_scaling/run_kochanowski_destress_analysis.R")
  run("analysis/dryad_global_regulators/run_weak_stress_window_pipeline.R")
  run("analysis/dryad_global_regulators/run_weak_stress_window_calibrated_pipeline.R")
  run("analysis/dryad_global_regulators/run_weak_stress_window_rank1_adjustment.R")
}
if (mode %in% c("--diagnostics", "--all")) {
  run("analysis/ecoli_promoter_screen/compute_paper_diagnostics.R")
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
    "analysis/french_pfibox/plot_clustering_and_inference.R",
    "analysis/french_pfibox/plot_pfibox_moderation_diagnostic.R",
    "analysis/kochanowski_promoter_scaling/plot_paper_figures.R",
    "analysis/dryad_global_regulators/plot_calibration_figures.R",
    "analysis/dryad_global_regulators/plot_rank_adjustment_figures.R",
    "analysis/dryad_global_regulators/plot_weak_stress_three_target_inference.R")
  for (script in scripts) run(script, cluster)
}
