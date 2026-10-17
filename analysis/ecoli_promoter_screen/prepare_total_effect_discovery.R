# Reuse saved fits; response construction and model fitting are unchanged.
source(file.path("analysis", "_helpers.R"))
source(file.path("analysis", "_paper_artifacts.R"))
load_destress_package()

source_dir <- analysis_output_dir("binsfeld_evc_calibrated")
out_dir <- analysis_output_dir("binsfeld_total_effect_discovery")
original <- read_paper_artifact(file.path(source_dir, "binsfeld_discovery.rds"))
comparison <- original$annotations$comparison
comparison <- comparison[, !grepl("^fixed_exponent_|^three_method_class$|^evc_huber_class$",
                                  names(comparison)), drop = FALSE]
fit_paths <- c(
  estimated_exponent = analysis_path("analysis", "outputs", "binsfeld_workflow_comparison",
                                     "destress_estimated_exponent_fit.rds"),
  evc_huber = file.path(source_dir, "evc_huber_fit.rds")
)
tables <- lapply(fit_paths, function(path) results(read_analysis_rds(path)))
key <- function(reporter, perturbation) paste(reporter, perturbation, sep = "::")
comparison_key <- key(comparison$promoter, comparison$compound)
stopifnot(!anyDuplicated(comparison_key))
target_counts <- list()
for (method in names(tables)) {
  table <- tables[[method]]
  result_key <- key(table$reporter, table$perturbation)
  stopifnot(!anyDuplicated(result_key), setequal(comparison_key, result_key))
  table <- table[match(comparison_key, result_key), , drop = FALSE]
  expected_padj <- ave(table$total_pvalue, table$reporter,
                       FUN = function(p) stats::p.adjust(p, method = "BH"))
  stopifnot(isTRUE(all.equal(unname(expected_padj),
                            table$total_padj_by_reporter, tolerance = 1e-12)))
  calls <- call_hits(table, fdr = 0.05, effect = "total_effect",
                     padj = "total_padj_by_reporter")$hit
  old_hits <- as.logical(comparison[[paste0(method, "_hit")]])
  comparison[[paste0(method, "_effect")]] <- table$total_effect
  comparison[[paste0(method, "_pvalue")]] <- table$total_pvalue
  comparison[[paste0(method, "_padj_by_reporter")]] <- table$total_padj_by_reporter
  comparison[[paste0(method, "_hit_class")]] <- calls
  comparison[[paste0(method, "_hit")]] <- calls != "Not DE"
  comparison[[paste0("neglog10_", method)]] <-
    -log10(pmax(table$total_pvalue, .Machine$double.xmin))
  new_hits <- comparison[[paste0(method, "_hit")]]
  target_counts[[method]] <- data.frame(
    method = method, pairs = nrow(table), relative_hits = sum(old_hits),
    total_hits = sum(new_hits), shared = sum(old_hits & new_hits),
    total_only = sum(!old_hits & new_hits), relative_only = sum(old_hits & !new_hits)
  )
}
b <- as.logical(comparison$binsfeld_hit)
n <- comparison$estimated_exponent_hit
e <- comparison$evc_huber_hit
comparison$evc_huber_class <- vapply(seq_len(nrow(comparison)), function(i) {
  paste(c(if (b[i]) "Binsfeld reference",
          if (n[i]) "DStressR without EVC",
          if (e[i]) "DStressR with EVC"), collapse = " + ")
}, character(1))
region_counts <- c(
  binsfeld_only = sum(b & !n & !e),
  estimated_exponent_only = sum(!b & n & !e),
  evc_huber_only = sum(!b & !n & e),
  binsfeld_estimated_exponent = sum(b & n & !e),
  binsfeld_evc_huber = sum(b & !n & e),
  estimated_exponent_evc_huber = sum(!b & n & e),
  all_three = sum(b & n & e)
)
set_counts <- c(binsfeld = sum(b), estimated_exponent = sum(n), evc_huber = sum(e))
stopifnot(sum(region_counts) == sum(b | n | e))
summary <- do.call(rbind, target_counts)
rownames(summary) <- NULL
write_paper_artifact(file.path(out_dir, "binsfeld_discovery.rds"),
  results = tables,
  parameters = list(target = "total", empirical_bayes = TRUE,
                    adjustment = "by_reporter", fdr = 0.05),
  annotations = list(comparison = comparison, region_counts = region_counts,
                     set_counts = set_counts, summary = summary),
  provenance = list(fit_paths = fit_paths, original_artifact =
                      file.path(source_dir, "binsfeld_discovery.rds")))
utils::write.table(comparison, file.path(out_dir, "total_effect_comparison.tsv"),
                    sep = "\t", quote = FALSE, row.names = FALSE)
utils::write.table(summary, file.path(out_dir, "total_vs_specific_hit_counts.tsv"),
                    sep = "\t", quote = FALSE, row.names = FALSE)
union <- comparison[b | n | e, , drop = FALSE]
utils::write.table(union, file.path(out_dir, "total_effect_significant_union.tsv"),
                    sep = "\t", quote = FALSE, row.names = FALSE)
# Prior pair/pathway annotations do not establish the direction of a total effect.
annotations <- read_tsv_base(analysis_path("analysis", "outputs",
  "binsfeld_workflow_comparison", "dstressr_only_literature_support.tsv"))
annotation_key <- key(annotations$promoter, annotations$compound)
stopifnot(!anyDuplicated(annotation_key))
candidates <- comparison[!b & (n | e), , drop = FALSE]
index <- match(key(candidates$promoter, candidates$compound), annotation_key)
for (column in c("literature_support", "support_note", "source")) {
  candidates[[column]] <- annotations[[column]][index]
}
unannotated <- is.na(index)
candidates$literature_support[unannotated] <- "Unannotated"
candidates$support_note[unannotated] <-
  "Not covered by the existing pair-level literature annotation table."
candidates$source[unannotated] <- "DStressR analysis"
candidates$annotation_scope <- "Pair/pathway context; total-effect direction requires separate assessment"
utils::write.table(candidates, file.path(out_dir, "total_effect_dstressr_only_annotations.tsv"),
                    sep = "\t", quote = FALSE, row.names = FALSE)
per_reporter <- do.call(rbind, lapply(split(comparison, comparison$promoter), function(d) {
  data.frame(reporter = d$promoter[1], reference = sum(d$binsfeld_hit),
    without_EVC = sum(d$estimated_exponent_hit), with_EVC = sum(d$evc_huber_hit))
}))
utils::write.table(per_reporter, file.path(out_dir, "total_effect_hits_by_reporter.tsv"),
                    sep = "\t", quote = FALSE, row.names = FALSE)
print(summary, row.names = FALSE)
print(region_counts)
message("Prepared total-effect discovery artifact from saved fits: ", out_dir)
