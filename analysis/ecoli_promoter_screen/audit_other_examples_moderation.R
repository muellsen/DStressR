# Extend the frozen-code inference comparison to the three other manuscript examples.
source(file.path("analysis", "_helpers.R"))
root <- analysis_project_root()
audit <- analysis_output_dir("inference_audit_2026-10-16")
out_dir <- file.path(audit, "other_examples")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
load_sources <- function(path) {
  e <- new.env(parent = .GlobalEnv)
  for (f in list.files(path, pattern = "[.]R$", full.names = TRUE)) sys.source(f, e)
  e
}
old <- load_sources(file.path(audit, "baseline_R"))
new <- load_sources(file.path(root, "R"))
read_tab <- function(path) read_analysis_delim(path, check.names = FALSE)
write_tab <- function(x, name) write.table(x, file.path(out_dir, name),
                                          sep = "\t", quote = FALSE, row.names = FALSE)
key <- function(x) paste(x$reporter, x$perturbation, sep = "::")
paths <- list(
  pfibox_response = analysis_path("analysis", "outputs", "french_pfibox", "promoter_activity_analysis",
    "pfibox_pa_response_table_reporters2000_perturbations100_T2-18h_treatment_replicate_scalable.tsv"),
  pfibox_result = analysis_path("analysis", "outputs", "french_pfibox", "promoter_activity_analysis",
    "pfibox_pa_destress_pair_results_reporters2000_perturbations100_T2-18h_treatment_replicate_scalable.tsv"),
  dash_input = analysis_path("analysis", "outputs", "dryad_global_regulators",
    "dryad_weak_stress_windows_calibrated_alpha_input.tsv"),
  dash_alpha = analysis_path("analysis", "outputs", "dryad_global_regulators",
    "dryad_weak_stress_windows_calibrated_alpha_supplied.tsv"),
  dash_result = analysis_path("analysis", "outputs", "dryad_global_regulators",
    "dryad_weak_stress_windows_calibrated_alpha_pair_results.tsv"),
  dash_rank1 = analysis_path("analysis", "outputs", "dryad_global_regulators",
    "dryad_weak_stress_windows_calibrated_alpha_rank1_pair_results.tsv"),
  kochanowski_input = analysis_path("analysis", "outputs", "kochanowski_promoter_scaling",
    "kochanowski_supplied_promoter_activity_long.tsv"),
  kochanowski_result = analysis_path("analysis", "outputs", "kochanowski_promoter_scaling",
    "kochanowski_destress_summary_effects.tsv")
)
stopifnot(all(file.exists(unlist(paths))))
archive_paths <- sub(paste0(root, "/"), "", unlist(paths), fixed = TRUE)
archive_paths <- sub("promoter_activity_analysis", "promoter_activity_subset",
                     archive_paths, fixed = TRUE)
original_dir <- tempfile("dstressr-original-examples-")
dir.create(original_dir)
utils::untar(file.path(audit, "other_examples_prechange_outputs.tar.gz"),
             files = archive_paths, exdir = original_dir)
paths <- as.list(stats::setNames(file.path(original_dir, archive_paths), names(paths)))
stopifnot(all(file.exists(unlist(paths))))
original_hashes <- tools::md5sum(unlist(paths))

get_fit <- function(name, build) {
  path <- file.path(out_dir, paste0(name, "_baseline_fit.rds"))
  if (file.exists(path)) return(read_analysis_rds(path))
  fit <- build()
  saveRDS(fit, path)
  fit
}

message("PFIbox: reconstructing the effect fit from the preserved scalar responses.")
pfibox_fit <- get_fit("pfibox", function() {
  d <- read_tab(paths$pfibox_response)
  assay <- old$prepare_assay(d, reporter = "reporter", perturbation = "perturbation",
    control = "Control", response = "y_pa", batch = "dose_code_centered",
    plate = "drug_factor", replicate = "replicate", numeric_covariates = "dose_code_centered")
  old$fit_destress(assay, technical = "replicate", empirical_bayes = TRUE,
                   adjustment = "by_reporter", interaction = FALSE, verbose = TRUE)
})
message("PFIbox: evaluating old and candidate inference on identical coefficients.")
pfibox_old <- old$results(pfibox_fit)
pfibox_new <- new$results(pfibox_fit)
focal <- load_sources(file.path(root, "R"))
focal$centered_contrast_df <- function(variance, df) df
pfibox_focal <- focal$results(pfibox_fit)
ablation <- data.frame(method = c("Prechange", "Single moderation, focal df", "Single moderation, effective df"),
  reporter_wise_hits = c(sum(pfibox_old$relative_padj_by_reporter < 0.05),
                         sum(pfibox_focal$relative_padj_by_reporter < 0.05),
                         sum(pfibox_new$relative_padj_by_reporter < 0.05)),
  global_hits = c(sum(pfibox_old$relative_padj_global < 0.05),
                   sum(pfibox_focal$relative_padj_global < 0.05),
                   sum(pfibox_new$relative_padj_global < 0.05)))
write_tab(ablation, "pfibox_moderation_ablation.tsv")

message("Dash: reconstructing the effect fit using the saved calibrated exponents.")
dash_fit <- get_fit("dash", function() {
  d <- read_tab(paths$dash_input)
  alpha <- read_tab(paths$dash_alpha)
  promoters <- c("Fur", "MarA", "SoxS", "LexA")
  windows <- c("Early", "Middle", "Late")
  reporter_levels <- unlist(lapply(promoters, function(a) paste(a, windows, sep = " | ")))
  d$pseudo_reporter <- factor(d$pseudo_reporter, levels = reporter_levels)
  d$compound <- factor(d$compound, levels = c("Standard", "Iron", "Tetracycline", "H2O2", "Kanamycin"))
  assay <- old$prepare_assay(d, reporter = "pseudo_reporter", perturbation = "compound",
    control = "Standard", lux = "gfp_auc", growth = "od_auc", replicate = "replicate",
    growth_exponent = setNames(alpha$alpha_calibrated, alpha$pseudo_reporter))
  old$fit_destress(assay, technical = "replicate", empirical_bayes = TRUE,
                   adjustment = "by_reporter", interaction = FALSE)
})
dash_rank1_fit <- dash_fit
dash_rank1_fit$background_rank <- 1L
dash_old <- old$results(dash_fit)
dash_new <- new$results(dash_fit)
dash_rank_old <- old$results(dash_rank1_fit)
dash_rank_new <- new$results(dash_rank1_fit)

message("Kochanowski: independently rebuilding the unchanged supplied-summary path.")
summary_data <- read_tab(paths$kochanowski_input)
summary_fit <- function(e) e$fit_destress_summary(summary_data, response = "response",
  response_variance = "response_variance", reference_perturbation = "glucose",
  background_rank = 1, global_method = "weighted_mean", adjustment = "global")
kochanowski_old <- old$results(summary_fit(old))
kochanowski_new <- new$results(summary_fit(new))

comparisons <- list()
pair_tables <- list()
compare <- function(example, target, before, after, archived, adjustment) {
  if (example == "Kochanowski") before <- before[before$perturbation != "glucose", ]
  after <- after[match(key(before), key(after)), ]
  archived <- archived[match(key(before), key(archived)), ]
  effect <- paste0(target, "_effect")
  p <- paste0(target, "_pvalue")
  q <- paste0(target, "_padj_", adjustment)
  stopifnot(nrow(after) == nrow(before), !anyNA(after$reporter), !anyNA(archived$reporter))
  effect_change <- max(abs(after[[effect]] - before[[effect]]))
  stopifnot(effect_change < 1e-12)
  old_hit <- before[[q]] < 0.05
  new_hit <- after[[q]] < 0.05
  archived_hit <- archived[[q]] < 0.05
  label <- paste(example, target, adjustment, sep = "_")
  pairs <- data.frame(reporter = before$reporter, perturbation = before$perturbation,
    example = example, target = target, adjustment = adjustment, effect = before[[effect]],
    archived_pvalue = archived[[p]], archived_padj = archived[[q]],
    baseline_pvalue = before[[p]], baseline_padj = before[[q]],
    candidate_pvalue = after[[p]], candidate_padj = after[[q]],
    archived_hit = archived_hit, baseline_hit = old_hit, candidate_hit = new_hit)
  pairs$change_vs_baseline <- ifelse(new_hit & !old_hit, "Gained",
                                    ifelse(!new_hit & old_hit, "Lost", "Unchanged"))
  pairs$change_vs_archive <- ifelse(new_hit & !archived_hit, "Gained",
                                   ifelse(!new_hit & archived_hit, "Lost", "Unchanged"))
  pair_tables[[label]] <<- pairs
  write_tab(pairs, paste0(label, "_pair_comparison.tsv"))
  write_tab(pairs[pairs$change_vs_archive != "Unchanged", ], paste0(label, "_changed_pairs.tsv"))
  comparisons[[label]] <<- data.frame(example = example, target = target,
    adjustment = adjustment, n_pairs = nrow(before),
    archived_hits = sum(archived_hit), baseline_hits = sum(old_hit), candidate_hits = sum(new_hit),
    gained_vs_baseline = sum(new_hit & !old_hit), lost_vs_baseline = sum(!new_hit & old_hit),
    gained_vs_archive = sum(new_hit & !archived_hit), lost_vs_archive = sum(!new_hit & archived_hit),
    max_abs_effect_difference = effect_change,
    max_abs_p_difference = max(abs(after[[p]] - before[[p]])),
    median_abs_log10_p_difference = median(abs(log10(pmax(after[[p]], 1e-300)) - log10(pmax(before[[p]], 1e-300)))),
    archive_max_abs_effect_difference = max(abs(before[[effect]] - archived[[effect]])),
    archive_max_abs_p_difference = max(abs(before[[p]] - archived[[p]])))
}

pfibox_archive <- read_tab(paths$pfibox_result)
dash_archive <- read_tab(paths$dash_result)
dash_rank_archive <- read_tab(paths$dash_rank1)
kochanowski_archive <- read_tab(paths$kochanowski_result)
rank_ix <- match(key(dash_rank_archive), key(dash_rank_old))
fit_ix <- match(key(dash_rank_archive), key(dash_fit$reporter_effects))
historical_rank_p <- 2 * pt(abs(dash_rank_old$rank_adjusted_total_statistic[rank_ix]),
                           dash_fit$reporter_effects$residual_df[fit_ix], lower.tail = FALSE)
rank_check <- data.frame(convention = "Unaugmented residual df for rank-adjusted total test",
  max_abs_p_difference = max(abs(historical_rank_p - dash_rank_archive$rank_adjusted_total_pvalue)))
write_tab(rank_check, "dash_historical_rank_reconstruction.tsv")
stopifnot(rank_check$max_abs_p_difference < 1e-10)
for (target in c("total", "relative")) {
  compare("PFIbox", target, pfibox_old, pfibox_new, pfibox_archive, "by_reporter")
  compare("Dash", target, dash_old, dash_new, dash_archive, "by_reporter")
}
compare("PFIbox", "relative", pfibox_old, pfibox_new, pfibox_archive, "global")
compare("Dash", "rank_adjusted_total", dash_rank_old, dash_rank_new, dash_rank_archive, "by_reporter")
for (target in c("total", "relative", "rank_adjusted_total")) {
  compare("Kochanowski", target, kochanowski_old, kochanowski_new, kochanowski_archive, "global")
}
summary <- do.call(rbind, comparisons)
write_tab(summary, "other_examples_inference_summary.tsv")
saveRDS(pair_tables, file.path(out_dir, "all_pair_comparisons.rds"))
saveRDS(list(pfibox = pfibox_new, dash_rank0 = dash_new, dash_rank1 = dash_rank_new,
             kochanowski = kochanowski_new), file.path(out_dir, "candidate_results.rds"))
preservation <- data.frame(path = names(original_hashes),
  unchanged = unname(original_hashes) == unname(tools::md5sum(unlist(paths))))
write_tab(preservation, "source_preservation_checks.tsv")
stopifnot(all(preservation$unchanged))

if (requireNamespace("ggplot2", quietly = TRUE)) {
  d <- do.call(rbind, pair_tables[c("PFIbox_relative_by_reporter", "Dash_relative_by_reporter",
                                    "Kochanowski_relative_global")])
  plot_data <- rbind(data.frame(example = d$example, method = "Prechange code", pvalue = d$baseline_pvalue),
                     data.frame(example = d$example, method = "Candidate", pvalue = d$candidate_pvalue))
  plot_data$example <- factor(plot_data$example, levels = c("PFIbox", "Dash", "Kochanowski"))
  plot_data$method <- factor(plot_data$method, levels = c("Prechange code", "Candidate"))
  hist <- ggplot2::ggplot(plot_data, ggplot2::aes(pvalue)) +
    ggplot2::geom_histogram(binwidth = 0.05, boundary = 0, fill = "#AEB5BC", colour = "white") +
    ggplot2::facet_grid(example ~ method, scales = "free_y") +
    ggplot2::labs(x = "Raw p-value", y = "Number of reporter-perturbation pairs") +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), strip.background = ggplot2::element_blank(),
                    strip.text = ggplot2::element_text(face = "bold"))
  ggplot2::ggsave(file.path(out_dir, "specific_pvalue_histograms.png"), hist,
                  width = 8.5, height = 8, dpi = 250)
}
print(summary, row.names = FALSE)
message("Original example inputs and results remain unchanged.")
