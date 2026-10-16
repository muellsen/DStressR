source(file.path("analysis", "_helpers.R"))
source(file.path("analysis", "_paper_artifacts.R"))
load_destress_package()
out_dir <- analysis_output_dir("binsfeld_effect_matrix_views")
without <- read_tsv_base(analysis_path("analysis", "outputs", "binsfeld_workflow_comparison",
                                      "destress_estimated_exponent_full_results.tsv"))
with <- read_tsv_base(analysis_path("analysis", "outputs", "binsfeld_evc_calibrated",
                                   "evc_huber_full_results.tsv"))
comparison <- read_tsv_base(analysis_path("analysis", "outputs", "binsfeld_evc_calibrated",
                                         "evc_huber_comparison_to_binsfeld_and_default.tsv"))
as_display <- function(table, workflow) data.frame(
  workflow = workflow, promoter = table$reporter, compound = table$perturbation,
  mean_response = table$total_effect)
display <- rbind(as_display(without, "DStressR without EVC"),
                 as_display(with, "DStressR with EVC"))
diagnose <- function(table, workflow) {
  d <- perturbation_diagnostics(table, mean_effect = "total_effect",
                                variance_effect = "total_effect", trend_span = 0.45)
  d$workflow <- workflow
  d
}
total_diag <- rbind(diagnose(without, "DStressR without EVC"),
                    diagnose(with, "DStressR with EVC"))
rownames(total_diag) <- NULL
write_paper_artifact(file.path(out_dir, "binsfeld_effect_paper.rds"),
  results = list(without_EVC = without, with_EVC = with),
  diagnostics = list(display = display, total_diag = total_diag),
  annotations = list(comparison = comparison))
