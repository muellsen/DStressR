source(file.path("analysis", "_helpers.R"))
load_destress_package()

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package `ggplot2` is required for interaction-model plots.", call. = FALSE)
}

ggplot2 <- asNamespace("ggplot2")

load_binsfeld_paper_data()

out_dir <- analysis_output_dir("ecoli_interaction_model")
ev_dir <- analysis_output_dir("binsfeld_evc_calibrated")

wt_auc <- binsfeld_reporter_auc[
  binsfeld_reporter_auc$strain == "WT" &
    binsfeld_reporter_auc$removed == "No",
]

safe_neglog10 <- function(x) -log10(pmax(as.numeric(x), .Machine$double.xmin))

fit_interaction <- function(use_ev = FALSE) {
  analysis_data <- if (use_ev) wt_auc else wt_auc[wt_auc$promoter != "EVC", , drop = FALSE]
  assay <- prepare_assay(
    analysis_data,
    reporter = "promoter",
    perturbation = "compound",
    control = "Water",
    lux = "lux_auc",
    growth = "od_auc",
    growth_exponent = "estimate",
    batch = "dose_level",
    replicate = "replicate",
    growth_covariates = "replicate",
    numeric_covariates = "dose_level",
    background_reporter = if (use_ev) "EVC" else NULL,
    background_method = if (use_ev) "huber" else "none",
    background_by = if (use_ev) c("compound", "dose_level", "replicate") else NULL
  )

  fit <- fit_destress(
    assay,
    technical = c("replicate", "dose_level"),
    empirical_bayes = TRUE,
    adjustment = "by_reporter",
    interaction = TRUE
  )

  saveRDS(fit, file.path(out_dir, if (use_ev) {
    "interaction_with_ev_fit.rds"
  } else {
    "interaction_without_ev_fit.rds"
  }))
  res <- results(fit)
  res$promoter <- res$reporter
  res$compound <- res$perturbation
  hit_class <- call_hits(
    res,
    fdr = 0.05,
    effect = "specific_effect",
    padj = "specific_padj_by_reporter"
  )$hit
  res$hit_class <- hit_class
  res$hit <- hit_class != "Not DE"
  res <- res[order(res$promoter, res$compound), , drop = FALSE]

  list(
    fit = fit,
    results = res,
    model_matrix = model_matrix_report(fit)
  )
}

with_ev <- fit_interaction(use_ev = TRUE)
no_ev <- fit_interaction(use_ev = FALSE)

write_pair_results <- function(x, prefix) {
  res <- x$results
  out <- res[, c(
    "promoter", "compound",
    "total_effect", "global_effect", "low_rank_effect",
    "specific_effect", "specific_se", "specific_statistic", "specific_pvalue",
    "specific_padj_by_reporter", "hit", "hit_class"
  )]
  utils::write.table(
    out,
    file.path(out_dir, paste0(prefix, "_pair_results.tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  utils::write.table(
    out[out$hit, , drop = FALSE],
    file.path(out_dir, paste0(prefix, "_significant_pairs.tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  matrix_report <- x$model_matrix
  workflow_suffix <- sub("^interaction_", "", prefix)
  matrix_report$analysis <- paste(matrix_report$model, workflow_suffix, sep = "_")
  utils::write.table(
    matrix_report[, c("analysis", setdiff(names(matrix_report), "analysis"))],
    file.path(out_dir, paste0(prefix, "_model_matrix.tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

write_pair_results(with_ev, "interaction_with_ev")
write_pair_results(no_ev, "interaction_without_ev")

comparison_path <- file.path(ev_dir, "evc_huber_comparison_to_binsfeld_and_default.tsv")
if (!file.exists(comparison_path)) {
  stop("Run analysis/ecoli_promoter_screen/run_evc_calibrated_analysis.R before this script.", call. = FALSE)
}
comparison <- utils::read.delim(comparison_path, check.names = FALSE)
comparison$binsfeld_hit <- as.logical(comparison$binsfeld_hit)
comparison$estimated_exponent_hit <- as.logical(comparison$estimated_exponent_hit)
comparison$evc_huber_hit <- as.logical(comparison$evc_huber_hit)

no_ev_small <- no_ev$results[, c(
  "promoter", "compound", "specific_effect", "specific_pvalue",
  "specific_padj_by_reporter", "hit", "hit_class"
)]
names(no_ev_small) <- c(
  "promoter", "compound", "interaction_without_ev_effect",
  "interaction_without_ev_pvalue", "interaction_without_ev_padj_by_reporter",
  "interaction_without_ev_hit", "interaction_without_ev_hit_class"
)

with_ev_small <- with_ev$results[, c(
  "promoter", "compound", "specific_effect", "specific_pvalue",
  "specific_padj_by_reporter", "hit", "hit_class"
)]
names(with_ev_small) <- c(
  "promoter", "compound", "interaction_with_ev_effect",
  "interaction_with_ev_pvalue", "interaction_with_ev_padj_by_reporter",
  "interaction_with_ev_hit", "interaction_with_ev_hit_class"
)

comparison <- merge(comparison, no_ev_small, by = c("promoter", "compound"), all.x = TRUE, sort = FALSE)
comparison <- merge(comparison, with_ev_small, by = c("promoter", "compound"), all.x = TRUE, sort = FALSE)
comparison$interaction_without_ev_hit[is.na(comparison$interaction_without_ev_hit)] <- FALSE
comparison$interaction_with_ev_hit[is.na(comparison$interaction_with_ev_hit)] <- FALSE
comparison$neglog10_interaction_without_ev <- safe_neglog10(comparison$interaction_without_ev_pvalue)
comparison$neglog10_interaction_with_ev <- safe_neglog10(comparison$interaction_with_ev_pvalue)

utils::write.table(
  comparison,
  file.path(out_dir, "interaction_comparison_to_primary_methods.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

hit_cols <- c(
  binsfeld = "binsfeld_hit",
  dstressr_without_ev = "estimated_exponent_hit",
  dstressr_with_ev = "evc_huber_hit",
  interaction_without_ev = "interaction_without_ev_hit",
  interaction_with_ev = "interaction_with_ev_hit"
)

summary_rows <- data.frame(
  metric = character(),
  count = integer(),
  stringsAsFactors = FALSE
)
for (nm in c(
  "binsfeld",
  "dstressr_with_ev",
  "interaction_with_ev",
  "dstressr_without_ev",
  "interaction_without_ev"
)) {
  summary_rows <- rbind(
    summary_rows,
    data.frame(metric = paste(nm, "hits"), count = sum(comparison[[hit_cols[[nm]]]], na.rm = TRUE))
  )
}
summary_rows <- rbind(
  summary_rows,
  data.frame(
    metric = c(
      "interaction_without_ev and dstressr_without_ev overlap",
      "interaction_with_ev and dstressr_with_ev overlap",
      "all primary methods and interaction variants overlap",
      "binsfeld hits recovered by interaction_without_ev",
      "binsfeld hits recovered by interaction_with_ev",
      "dstressr_without_ev hits recovered by interaction_without_ev",
      "dstressr_with_ev hits recovered by interaction_with_ev"
    ),
    count = c(
      sum(comparison$interaction_without_ev_hit & comparison$estimated_exponent_hit, na.rm = TRUE),
      sum(comparison$interaction_with_ev_hit & comparison$evc_huber_hit, na.rm = TRUE),
      sum(comparison$binsfeld_hit & comparison$estimated_exponent_hit & comparison$evc_huber_hit &
            comparison$interaction_without_ev_hit & comparison$interaction_with_ev_hit, na.rm = TRUE),
      sum(comparison$binsfeld_hit & comparison$interaction_without_ev_hit, na.rm = TRUE),
      sum(comparison$binsfeld_hit & comparison$interaction_with_ev_hit, na.rm = TRUE),
      sum(comparison$estimated_exponent_hit & comparison$interaction_without_ev_hit, na.rm = TRUE),
      sum(comparison$evc_huber_hit & comparison$interaction_with_ev_hit, na.rm = TRUE)
    )
  )
)

utils::write.table(
  summary_rows,
  file.path(out_dir, "interaction_hit_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

source(file.path("analysis", "_paper_artifacts.R"))
write_paper_artifact(file.path(out_dir, "binsfeld_interaction_paper.rds"),
  results = list(without_EVC = no_ev$results, with_EVC = with_ev$results),
  annotations = list(comparison = comparison, hit_cols = hit_cols))
