# Re-extract ordinary tests from the saved fit; no response reconstruction or refit.
source(file.path("analysis", "_helpers.R"))
load_destress_package()
out <- analysis_output_dir("french_pfibox", "temporal_unscaled_integral_analysis")
tag <- "reporters2000_perturbations100_Tall_available_treatment_replicate_scalable_temporal_unscaled_integral"
fit <- readRDS(file.path(out, paste0("pfibox_pa_fit_", tag, ".rds")))
fit$empirical_bayes <- FALSE
ordinary <- results(fit)
moderated <- read.delim(file.path(out, paste0("pfibox_pa_destress_pair_results_", tag, ".tsv")))
stopifnot(identical(ordinary$reporter, moderated$reporter),
  identical(ordinary$perturbation, moderated$perturbation),
  isTRUE(all.equal(ordinary$total_effect, moderated$total_effect)))
write.table(ordinary, file.path(out, paste0("pfibox_pa_destress_pair_results_", tag, "_rawse.tsv")),
  sep = "\t", quote = FALSE, row.names = FALSE)
hits <- function(x, target, scope) x[[paste0(target, "_padj_", scope)]] < .05
summary <- data.frame(method = c("Ordinary", "Moderated"),
  total_global = c(sum(hits(ordinary, "total", "global")), sum(hits(moderated, "total", "global"))),
  total_by_reporter = c(sum(hits(ordinary, "total", "by_reporter")), sum(hits(moderated, "total", "by_reporter"))),
  relative_global = c(sum(hits(ordinary, "relative", "global")), sum(hits(moderated, "relative", "global"))),
  relative_by_reporter = c(sum(hits(ordinary, "relative", "by_reporter")), sum(hits(moderated, "relative", "by_reporter"))))
write.table(summary, file.path(out, "integrated_activity_testing_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE)
print(summary)
print(table(total = hits(moderated, "total", "global"), relative = hits(moderated, "relative", "global")))
