# Recover total-effect inference from saved contrasts, without refitting models.
source(file.path("analysis", "_helpers.R"))
source(file.path("analysis", "_paper_artifacts.R"))
load_destress_package()

source_dir <- analysis_output_dir("ecoli_interaction_model")
out_dir <- analysis_output_dir("ecoli_interaction_total_effects")
source_path <- file.path(source_dir, "binsfeld_interaction_paper.rds")
original <- read_paper_artifact(source_path)
primary <- read_paper_artifact(analysis_path("analysis", "outputs",
  "binsfeld_total_effect_discovery", "binsfeld_discovery.rds"))
comparison <- primary$annotations$comparison
key <- function(a, j) paste(a, j, sep = "::")
pair_key <- key(comparison$promoter, comparison$compound)
stopifnot(!anyDuplicated(pair_key))
variants <- c(without_EVC = "interaction_without_ev", with_EVC = "interaction_with_ev")
tables <- list()
summaries <- list()
for (variant in names(variants)) {
  prefix <- variants[[variant]]
  table <- original$results[[variant]]
  matrix_path <- file.path(source_dir, paste0(prefix, "_model_matrix.tsv"))
  design <- read_tsv_base(matrix_path)
  df <- design$residual_df[design$model == "interaction"]
  stopifnot(length(df) == 1L, is.finite(df), df > 0,
            all(table$low_rank_effect == 0))
  if (!"total_var_raw" %in% names(table)) {
    # Earlier artifacts stored ordinary total-effect SEs and p-values.
    ordinary <- wald_t_test(table$total_effect, table$total_se, df)
    stopifnot(isTRUE(all.equal(ordinary$pvalue, table$total_pvalue, tolerance = 1e-12)))
    table$total_var_raw <- table$total_se^2
  }
  table$residual_df <- rep(df, nrow(table))
  raw_se <- sqrt(table$total_var_raw)
  moderated <- eb_moderate_se(raw_se, table$residual_df)
  test <- wald_t_test(table$total_effect, moderated$se, moderated$df)
  table$total_se <- moderated$se
  table$total_df_for_test <- moderated$df
  table$total_statistic <- test$statistic
  table$total_pvalue <- test$pvalue
  table$rank_adjusted_total_se <- table$total_se
  table$rank_adjusted_total_statistic <- table$total_statistic
  table$rank_adjusted_total_pvalue <- table$total_pvalue
  for (target in c("total", "rank_adjusted_total")) {
    p <- table[[paste0(target, "_pvalue")]]
    table[[paste0(target, "_padj_global")]] <-
      adjust_destress_pvalues(p, table$reporter, "global")
    table[[paste0(target, "_padj_by_reporter")]] <-
      adjust_destress_pvalues(p, table$reporter, "by_reporter")
    table[[paste0(target, "_padj")]] <- table[[paste0(target, "_padj_by_reporter")]]
  }
  table <- table[, setdiff(names(table), c("promoter", "compound", "hit", "hit_class")),
                 drop = FALSE]
  result_key <- key(table$reporter, table$perturbation)
  stopifnot(!anyDuplicated(result_key), setequal(pair_key, result_key))
  table <- table[match(pair_key, result_key), , drop = FALSE]
  calls <- call_hits(table, fdr = 0.05, effect = "total_effect",
                     padj = "total_padj_by_reporter")$hit
  comparison[[paste0(prefix, "_effect")]] <- table$total_effect
  comparison[[paste0(prefix, "_pvalue")]] <- table$total_pvalue
  comparison[[paste0(prefix, "_padj_by_reporter")]] <- table$total_padj_by_reporter
  comparison[[paste0(prefix, "_hit")]] <- calls != "Not DE"
  comparison[[paste0(prefix, "_hit_class")]] <- calls
  comparison[[paste0("neglog10_", prefix)]] <-
    -log10(pmax(table$total_pvalue, .Machine$double.xmin))
  ordinary <- wald_t_test(table$total_effect, raw_se[match(pair_key, result_key)], df)
  ordinary_q <- adjust_destress_pvalues(ordinary$pvalue, table$reporter, "by_reporter")
  baseline <- if (variant == "without_EVC") "estimated_exponent" else "evc_huber"
  primary_effect <- comparison[[paste0(baseline, "_effect")]]
  primary_hit <- comparison[[paste0(baseline, "_hit")]]
  summaries[[variant]] <- data.frame(
    variant = variant, pairs = nrow(table), residual_df = df,
    prior_df = moderated$prior_df, reporter_wise_hits = sum(primary_hit),
    interaction_hits = sum(calls != "Not DE"),
    ordinary_interaction_hits = sum(ordinary_q < 0.05),
    shared = sum(primary_hit & calls != "Not DE"),
    shared_with_reference = sum(comparison$binsfeld_hit & calls != "Not DE"),
    effect_correlation = stats::cor(primary_effect, table$total_effect),
    mean_absolute_difference = mean(abs(primary_effect - table$total_effect))
  )
  tables[[variant]] <- table
  utils::write.table(table, file.path(out_dir, paste0(prefix, "_total_results.tsv")),
                      sep = "\t", quote = FALSE, row.names = FALSE)
  stopifnot(file.copy(matrix_path, file.path(out_dir, basename(matrix_path)), overwrite = TRUE))
}
hit_cols <- c(binsfeld = "binsfeld_hit", dstressr_without_ev = "estimated_exponent_hit",
  dstressr_with_ev = "evc_huber_hit", interaction_without_ev = "interaction_without_ev_hit",
  interaction_with_ev = "interaction_with_ev_hit")
summary <- do.call(rbind, summaries)
rownames(summary) <- NULL
write_paper_artifact(file.path(out_dir, "binsfeld_interaction_paper.rds"),
  results = tables, parameters = list(target = "total", empirical_bayes = TRUE,
    adjustment = "by_reporter", fdr = 0.05),
  annotations = list(comparison = comparison, hit_cols = hit_cols, summary = summary),
  provenance = list(source_artifact = source_path,
    inference = "Moderated total tests reconstructed from saved raw contrast variances"))
utils::write.table(comparison, file.path(out_dir, "interaction_total_comparison.tsv"),
                    sep = "\t", quote = FALSE, row.names = FALSE)
utils::write.table(summary, file.path(out_dir, "interaction_total_summary.tsv"),
                    sep = "\t", quote = FALSE, row.names = FALSE)
print(summary, row.names = FALSE)
cat("All five methods:", sum(Reduce(`&`, comparison[unname(hit_cols)])), "\n")
