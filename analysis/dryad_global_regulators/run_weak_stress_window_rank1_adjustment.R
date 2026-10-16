source(file.path("analysis", "_helpers.R"))
source(file.path("analysis", "_paper_artifacts.R"))
load_destress_package()

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package `ggplot2` is required for rank-1 adjustment plots.", call. = FALSE)
}
if (!requireNamespace("gridExtra", quietly = TRUE)) {
  stop("Package `gridExtra` is required for rank-1 adjustment plots.", call. = FALSE)
}

ggplot2 <- asNamespace("ggplot2")
gridExtra <- asNamespace("gridExtra")

panel_label <- function(label, size = 15) {
  grid::textGrob(
    label,
    x = grid::unit(0, "npc"),
    y = grid::unit(0.08, "npc"),
    hjust = 0,
    vjust = 0,
    gp = grid::gpar(fontsize = size, fontface = "bold", col = "#111827")
  )
}
out_dir <- analysis_output_dir("dryad_global_regulators")

input_path <- file.path(out_dir, "dryad_weak_stress_windows_calibrated_alpha_input.tsv")
alpha_path <- file.path(out_dir, "dryad_weak_stress_windows_calibrated_alpha_supplied.tsv")
if (!file.exists(input_path) || !file.exists(alpha_path)) {
  stop(
    "Missing calibrated-alpha weak-stress inputs. Run ",
    "analysis/dryad_global_regulators/run_weak_stress_window_calibrated_pipeline.R first.",
    call. = FALSE
  )
}

screen <- utils::read.delim(input_path, check.names = FALSE)
alpha <- utils::read.delim(alpha_path, check.names = FALSE)
stress_alpha <- stats::setNames(alpha$alpha_calibrated, alpha$pseudo_reporter)

promoter_levels <- c("Fur", "MarA", "SoxS", "LexA")
window_levels <- c("Early", "Middle", "Late")
compound_levels <- c("Iron", "Tetracycline", "H2O2", "Kanamycin")
pseudo_reporter_levels <- unlist(
  lapply(promoter_levels, function(promoter) paste(promoter, window_levels, sep = " | ")),
  use.names = FALSE
)

screen$promoter <- factor(screen$promoter, levels = promoter_levels)
screen$window <- factor(screen$window, levels = window_levels)
screen$pseudo_reporter <- factor(screen$pseudo_reporter, levels = pseudo_reporter_levels)
screen$compound <- factor(screen$compound, levels = c("Standard", compound_levels))

assay <- prepare_assay(
  screen,
  reporter = "pseudo_reporter",
  perturbation = "compound",
  control = "Standard",
  lux = "gfp_auc",
  growth = "od_auc",
  growth_exponent = stress_alpha,
  replicate = "replicate"
)

fit_rank0 <- fit_destress(
  assay,
  technical = "replicate",
  empirical_bayes = TRUE,
  adjustment = "by_reporter",
  interaction = FALSE,
  background_rank = 0
)
fit_rank1 <- fit_destress(
  assay,
  technical = "replicate",
  empirical_bayes = TRUE,
  adjustment = "by_reporter",
  interaction = FALSE,
  background_rank = 1
)

annotate_results <- function(res, rank) {
  res$promoter <- res$reporter
  res$compound <- res$perturbation
  parts <- do.call(rbind, strsplit(as.character(res$promoter), " \\| "))
  res$base_promoter <- parts[, 1]
  res$window <- parts[, 2]
  hit_effect <- if (rank > 0) "rank_adjusted_total_effect" else "specific_effect"
  hit_padj <- if (rank > 0) "rank_adjusted_total_padj_by_reporter" else "specific_padj_by_reporter"
  hit_pvalue <- if (rank > 0) "rank_adjusted_total_pvalue" else "specific_pvalue"
  hit_table <- call_hits(
    res,
    fdr = 0.05,
    effect = hit_effect,
    padj = hit_padj
  )
  res$hit_class <- hit_table$hit
  res$hit <- res$hit_class != "Not DE"
  res$background_rank <- rank
  res$rank_target_effect <- res[[hit_effect]]
  res$rank_target_pvalue <- res[[hit_pvalue]]
  res$rank_target_padj <- res[[hit_padj]]
  res
}

saveRDS(list(rank0 = fit_rank0, rank1 = fit_rank1),
        file.path(out_dir, "dryad_weak_stress_windows_calibrated_alpha_rank_fits.rds"))
res_rank0 <- annotate_results(results(fit_rank0), 0)
res_rank1 <- annotate_results(results(fit_rank1), 1)
rank_results <- rbind(res_rank0, res_rank1)

rank_diagnostics <- background_rank_diagnostics(
  res_rank0,
  effect = "total_effect",
  reporter = "promoter",
  perturbation = "compound",
  rank_max = 4,
  permutations = 999,
  seed = 1
)
utils::write.table(
  rank_diagnostics,
  file.path(out_dir, "dryad_weak_stress_windows_calibrated_alpha_rank_scree_diagnostics.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

utils::write.table(
  res_rank1,
  file.path(out_dir, "dryad_weak_stress_windows_calibrated_alpha_rank1_pair_results.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
utils::write.table(
  res_rank1[res_rank1$hit, , drop = FALSE],
  file.path(out_dir, "dryad_weak_stress_windows_calibrated_alpha_rank1_significant_pairs.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

rank_summary <- do.call(rbind, lapply(split(rank_results, rank_results$background_rank), function(d) {
  data.frame(
    background_rank = unique(d$background_rank),
    tested_pairs = nrow(d),
    significant_pairs = sum(d$hit, na.rm = TRUE),
    median_pvalue = stats::median(d$rank_target_pvalue, na.rm = TRUE),
    min_pvalue = min(d$rank_target_pvalue, na.rm = TRUE),
    max_pvalue = max(d$rank_target_pvalue, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))
rank_summary <- rank_summary[order(rank_summary$background_rank), , drop = FALSE]
utils::write.table(
  rank_summary,
  file.path(out_dir, "dryad_weak_stress_windows_calibrated_alpha_rank_adjustment_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

expected <- data.frame(
  base_promoter = c("Fur", "MarA", "SoxS", "LexA"),
  compound = c("Iron", "Tetracycline", "H2O2", "Kanamycin"),
  expected_pair = TRUE,
  stringsAsFactors = FALSE
)
expected_results <- merge(rank_results, expected, by = c("base_promoter", "compound"))
expected_results <- expected_results[
  order(
    expected_results$background_rank,
    match(expected_results$base_promoter, promoter_levels),
    match(expected_results$window, window_levels)
  ),
  ,
  drop = FALSE
]
utils::write.table(
  expected_results,
  file.path(out_dir, "dryad_weak_stress_windows_calibrated_alpha_rank_expected_pairs.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write_paper_artifact(file.path(out_dir, "dash_rank_paper.rds"),
  results = list(rank0 = res_rank0, rank1 = res_rank1),
  diagnostics = list(rank_diagnostics = rank_diagnostics),
  annotations = list(expected = expected, promoter_levels = promoter_levels,
    window_levels = window_levels, compound_levels = compound_levels,
    pseudo_reporter_levels = pseudo_reporter_levels))
