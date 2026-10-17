source(file.path("analysis", "_helpers.R"))
source(file.path("analysis", "_paper_artifacts.R"))
scripts <- c(
  "analysis/ecoli_promoter_screen/plot_response_modeling.R",
  "analysis/ecoli_promoter_screen/plot_binsfeld_effect_matrix_views.R",
  "analysis/ecoli_promoter_screen/plot_binsfeld_effect_matrix_and_variance_diagnostic.R",
  "analysis/ecoli_promoter_screen/plot_discovery_figures.R",
  "analysis/ecoli_promoter_screen/plot_interaction_figures.R",
  "analysis/french_pfibox/plot_clustering_and_inference.R",
  "analysis/french_pfibox/plot_pfibox_model_cluster_and_discovery.R",
  "analysis/french_pfibox/plot_french_cluster_alluvial.R",
  "analysis/french_pfibox/plot_pfibox_moderation_diagnostic.R",
  "analysis/french_pfibox/activity_plotting.R",
  "analysis/french_pfibox/build_selected_activity_figure.R",
  "analysis/kochanowski_promoter_scaling/plot_paper_figures.R",
  "analysis/dryad_global_regulators/plot_calibration_figures.R",
  "analysis/dryad_global_regulators/plot_rank_adjustment_figures.R",
  "analysis/dryad_global_regulators/plot_weak_stress_three_target_inference.R")
forbidden <- c("prepare_assay", "fit_destress", "fit_destress_summary", "results",
               "estimate_growth_exponents", "estimate_eb_prior_variance", "eb_moderate_se",
               "wald_t_test", "p.adjust", "pt", "pnorm", "svd", "hclust",
               "background_rank_diagnostics", "low_rank_effect_decomposition", "effect_tail_scores")
calls <- function(x) {
  if (is.expression(x) || is.pairlist(x)) return(unlist(lapply(as.list(x), calls)))
  if (!is.call(x)) return(character())
  head <- x[[1]]
  name <- if (is.symbol(head)) as.character(head) else if (
    is.call(head) && as.character(head[[1]]) %in% c("::", ":::")) as.character(head[[3]]) else ""
  c(name, unlist(lapply(as.list(x)[-1], calls)))
}
for (path in scripts) {
  bad <- intersect(calls(parse(path)), forbidden)
  if (length(bad)) stop(path, " performs analysis: ", paste(bad, collapse = ", "))
  if (any(grepl("\\$mtime", readLines(path)))) stop(path, " selects inputs by modification time")
}
# Exercise the artifact schema without requiring any assay data.
path <- tempfile(fileext = ".rds")
table <- data.frame(reporter = "a", perturbation = "j", total_effect = 1)
write_paper_artifact(path, results = list(default = table))
stopifnot(identical(read_paper_artifact(path)$results$default, table))
historical <- data.frame(reporter = "a", perturbation = "j", specific_effect = 2)
write_paper_artifact(path, results = list(default = historical),
                     annotations = list(target = "specific"))
renamed <- read_paper_artifact(path)
stopifnot(identical(renamed$results$default$relative_effect, 2),
          !"specific_effect" %in% names(renamed$results$default),
          identical(renamed$annotations$target, "relative"))
bad <- try(write_paper_artifact(path, results = list(data.frame(promoter = "a"))), silent = TRUE)
stopifnot(inherits(bad, "try-error"))
message("Figure contract passed for ", length(scripts), " scripts.")
