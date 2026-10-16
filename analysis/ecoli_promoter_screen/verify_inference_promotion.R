# Verify regenerated manuscript results against the independently audited fits.
source(file.path("analysis", "_helpers.R"))
load_destress_package()
audit <- analysis_output_dir("inference_audit_2026-10-16")
read_output <- function(...) read.delim(analysis_path("analysis", "outputs", ...),
                                        check.names = FALSE)
key <- function(x) paste(x$reporter, x$perturbation, sep = "::")
checks <- list()
check <- function(name, actual, expected, columns) {
  stopifnot(!anyDuplicated(key(actual)), setequal(key(actual), key(expected)))
  actual <- actual[match(key(expected), key(actual)), ]
  for (column in columns) {
    difference <- max(abs(actual[[column]] - expected[[column]]), na.rm = TRUE)
    stopifnot(is.finite(difference), difference < 1e-9)
    checks[[paste(name, column)]] <<- data.frame(
      application = name, quantity = column, pairs = nrow(actual), max_abs_difference = difference)
  }
}
columns <- c("total_effect", "specific_effect", "specific_var", "specific_df_raw",
             "specific_df_for_test", "specific_pvalue", "specific_padj_by_reporter")
bins <- readRDS(file.path(audit, "all_result_tables.rds"))
check("Binsfeld without EVC",
      read_output("binsfeld_workflow_comparison", "destress_estimated_exponent_full_results.tsv"),
      bins[["without_EVC single_WS"]], columns)
check("Binsfeld with EVC",
      read_output("binsfeld_evc_calibrated", "evc_huber_full_results.tsv"),
      bins[["with_EVC single_WS"]], columns)
other <- readRDS(file.path(audit, "other_examples", "candidate_results.rds"))
tag <- "reporters2000_perturbations100_T2-18h_treatment_replicate_scalable"
check("PFIbox", read_output("french_pfibox", "promoter_activity_analysis",
      paste0("pfibox_pa_destress_pair_results_", tag, ".tsv")), other$pfibox, columns)
check("Dash", read_output("dryad_global_regulators",
      "dryad_weak_stress_windows_calibrated_alpha_pair_results.tsv"), other$dash_rank0, columns)
check("Dash rank 1", read_output("dryad_global_regulators",
      "dryad_weak_stress_windows_calibrated_alpha_rank1_pair_results.tsv"), other$dash_rank1,
      c(columns, "rank_adjusted_total_pvalue", "rank_adjusted_total_padj_by_reporter"))
koch <- read_output("kochanowski_promoter_scaling", "kochanowski_destress_summary_effects.tsv")
check("Kochanowski", koch[koch$perturbation != "glucose", ],
      other$kochanowski[other$kochanowski$perturbation != "glucose", ],
      c("total_effect", "specific_effect", "specific_pvalue", "rank_adjusted_total_pvalue"))
stopifnot(nrow(koch) == 1600L, length(unique(koch$reporter)) == 64L,
          length(unique(koch$perturbation)) == 25L)

summary <- do.call(rbind, checks)
write.table(summary, file.path(audit, "promotion_verification.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
sources <- c(list.files("R", "[.]R$", full.names = TRUE),
             "paper/dstressr_package_manuscript/manuscript.qmd", "paper/references.bib")
write.table(data.frame(path = sources, md5 = unname(tools::md5sum(sources))),
            file.path(audit, "promoted_source_manifest.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(audit, "promotion_sessionInfo.txt"))
manuscript <- readLines("paper/dstressr_package_manuscript/manuscript.qmd", warn = FALSE)
images <- regmatches(manuscript, regexpr("[.][.]/[.][.]/analysis/outputs/[^)]+", manuscript))
images <- unique(images[nzchar(images)])
images <- normalizePath(file.path("paper/dstressr_package_manuscript", images), mustWork = TRUE)
write.table(data.frame(path = images, md5 = unname(tools::md5sum(images))),
            file.path(audit, "promoted_figure_manifest.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
message("Verified ", nrow(summary), " regenerated quantities against audited results.")
