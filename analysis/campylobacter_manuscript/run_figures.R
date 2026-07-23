#!/usr/bin/env Rscript

source(file.path("analysis", "_helpers.R"))

project_root <- analysis_project_root()
scripts <- c(
  "analysis/campylobacter_manuscript/plot_growth_parameter_estimates.R",
  "analysis/campylobacter_manuscript/plot_response_heatmap_comparison.R",
  "analysis/campylobacter_manuscript/plot_moderated_hit_network_no_pcmeA_unique.R",
  "analysis/campylobacter_manuscript/plot_moderated_hit_bipartite_heatmap.R",
  "analysis/campylobacter_manuscript/plot_significance_summary.R"
)

for (script in scripts) {
  path <- file.path(project_root, script)
  if (!file.exists(path)) {
    stop("Missing manuscript figure script: ", script, call. = FALSE)
  }
}

message("Generating Campylobacter manuscript figures from package outputs.")
message("Project root: ", project_root)
message("Adjusted p-value family: ", comparison_adjustment())

old_wd <- getwd()
on.exit(setwd(old_wd), add = TRUE)
setwd(project_root)

for (script in scripts) {
  message("\n--- ", script, " ---")
  status <- system2(file.path(R.home("bin"), "Rscript"), script)
  if (!identical(status, 0L)) {
    stop("Figure script failed: ", script, call. = FALSE)
  }
}

message("\nDone. Figures were written under analysis/outputs/comparisons/.")
