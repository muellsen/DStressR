#!/usr/bin/env Rscript

# DStressR-style supplied-response analysis of Kochanowski et al. 2017.
#
# The public supplementary tables contain processed steady-state promoter
# activities and estimated standard deviations, not individual replicate wells.
# We therefore analyze total/reporter-relative effect matrices and low-rank diagnostics,
# while keeping the p-values below as summary-level approximations.

source(file.path("analysis", "_helpers.R"))
source(file.path("analysis", "_paper_artifacts.R"))

stopifnot(requireNamespace("readxl", quietly = TRUE))
stopifnot(requireNamespace("ggplot2", quietly = TRUE))
stopifnot(requireNamespace("gridExtra", quietly = TRUE))

if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(".", quiet = TRUE)
} else {
  library(DStressR)
}

source_dir <- file.path("data-raw", "kochanowski2017", "source")
ev2_file <- file.path(source_dir, "kochanowski2017_TableEV2.xlsx")
ev3_file <- file.path(source_dir, "kochanowski2017_TableEV3.xlsx")
ev4_file <- file.path(source_dir, "kochanowski2017_TableEV4.xlsx")
out_dir <- file.path("analysis", "outputs", "kochanowski_promoter_scaling")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

theme_koch <- destress_diagnostic_theme(base_size = 9)

read_condition_table <- function() {
  raw <- readxl::read_excel(ev2_file, sheet = 1, col_names = FALSE)
  d <- raw[7:32, 1:5]
  names(d) <- c(
    "condition",
    "perturbation",
    "growth_rate",
    "growth_rate_sd",
    "metabolomics_available"
  )
  d$condition_number <- seq_len(nrow(d))
  d$growth_rate <- as.numeric(d$growth_rate)
  d$growth_rate_sd <- as.numeric(d$growth_rate_sd)
  d$perturbation <- gsub("mueM", "uM", d$perturbation, fixed = TRUE)
  d
}

read_ev3_matrix_pair <- function(sheet, condition_numbers, condition_labels) {
  raw <- readxl::read_excel(ev3_file, sheet = sheet, col_names = FALSE)
  header_row <- which(raw[[1]] == "Promoter name")[1]
  if (!is.finite(header_row)) {
    stop("Could not find `Promoter name` header in ", sheet, call. = FALSE)
  }
  data_rows <- seq.int(header_row + 1L, nrow(raw))
  reporters <- as.character(raw[[1]][data_rows])
  keep <- !is.na(reporters) & nzchar(reporters)
  reporters <- reporters[keep]
  data_rows <- data_rows[keep]

  activity <- as.data.frame(raw[data_rows, seq_along(condition_numbers) + 1L])
  sd <- as.data.frame(raw[data_rows, seq_along(condition_numbers) + 32L])
  names(activity) <- condition_labels
  names(sd) <- condition_labels
  activity[] <- lapply(activity, as.numeric)
  sd[] <- lapply(sd, as.numeric)

  activity_long <- reshape(
    data.frame(reporter = reporters, activity, check.names = FALSE),
    direction = "long",
    varying = condition_labels,
    v.names = "promoter_activity",
    timevar = "perturbation",
    times = condition_labels
  )
  sd_long <- reshape(
    data.frame(reporter = reporters, sd, check.names = FALSE),
    direction = "long",
    varying = condition_labels,
    v.names = "promoter_activity_sd",
    timevar = "perturbation",
    times = condition_labels
  )
  row.names(activity_long) <- NULL
  row.names(sd_long) <- NULL
  merge(activity_long, sd_long, by = c("reporter", "perturbation"), sort = FALSE)
}

read_original_svd_axis <- function() {
  raw <- readxl::read_excel(ev4_file, sheet = "EV table 4 - 1", col_names = FALSE)
  d <- raw[6:31, c(1, 2, 3, 4, 5)]
  names(d) <- c(
    "condition_number",
    "perturbation",
    "growth_rate_ev4",
    "growth_rate_sd_ev4",
    "paper_first_singular_vector"
  )
  d$condition_number <- as.integer(d$condition_number)
  d$growth_rate_ev4 <- as.numeric(d$growth_rate_ev4)
  d$growth_rate_sd_ev4 <- as.numeric(d$growth_rate_sd_ev4)
  d$paper_first_singular_vector <- as.numeric(d$paper_first_singular_vector)
  d$perturbation <- gsub("mueM", "uM", d$perturbation, fixed = TRUE)
  d
}

read_paper_global_explanation <- function() {
  raw <- readxl::read_excel(ev4_file, sheet = "EV table 4 - 3", col_names = FALSE)
  header_row <- which(raw[[1]] == "Promoter")[1]
  if (!is.finite(header_row)) {
    stop("Could not find `Promoter` header in EV table 4 - 3.", call. = FALSE)
  }
  d <- raw[(header_row + 1L):nrow(raw), c(1, 3, 4)]
  names(d) <- c("reporter", "paper_global_correlation", "paper_global_pvalue")
  d$reporter <- as.character(d$reporter)
  d <- d[!is.na(d$reporter) & nzchar(d$reporter), , drop = FALSE]
  d$paper_global_correlation <- as.numeric(d$paper_global_correlation)
  d$paper_global_pvalue <- as.numeric(d$paper_global_pvalue)
  d$paper_global_r2 <- d$paper_global_correlation^2
  d
}

conditions <- read_condition_table()
condition_labels <- conditions$perturbation
activity <- read_ev3_matrix_pair(
  "EV table 3 - 1",
  condition_numbers = conditions$condition_number,
  condition_labels = condition_labels
)
activity <- merge(activity, conditions, by = "perturbation", all.x = TRUE, sort = FALSE)
activity$response <- log2(activity$promoter_activity)
activity$response_variance <-
  (activity$promoter_activity_sd / (activity$promoter_activity * log(2)))^2
activity$response_variance[!is.finite(activity$response_variance)] <- NA_real_

reference_perturbation <- "glucose"
analysis_rank <- 1L

fit_rank0 <- fit_destress_summary(
  activity,
  response = "response",
  response_variance = "response_variance",
  reference_perturbation = reference_perturbation,
  background_rank = 0,
  global_method = "weighted_mean",
  adjustment = "global"
)
effects_rank0 <- results(fit_rank0)
rank_diag <- background_rank_diagnostics(
  effects_rank0,
  effect = "total_effect",
  reporter = "reporter",
  perturbation = "perturbation",
  rank_max = 6,
  permutations = 500,
  threshold = 0.99,
  impute = "column_mean",
  seed = 1
)
rank_pass <- rank_diag$observed > rank_diag$null_threshold
first_fail <- match(FALSE, rank_pass)
selected_rank <- if (is.na(first_fail)) sum(rank_pass) else first_fail - 1L
selected_rank <- max(1L, selected_rank)

fit_summary <- fit_destress_summary(
  activity,
  response = "response",
  response_variance = "response_variance",
  reference_perturbation = reference_perturbation,
  background_rank = analysis_rank,
  global_method = "weighted_mean",
  adjustment = "global"
)
saveRDS(list(rank0 = fit_rank0, rank1 = fit_summary),
        file.path(out_dir, "kochanowski_summary_fits.rds"))
effects <- results(fit_summary)
effect_metadata <- activity[, c(
  "reporter", "perturbation", "promoter_activity", "promoter_activity_sd",
  "condition", "growth_rate", "growth_rate_sd", "metabolomics_available",
  "condition_number"
)]
effects <- merge(
  effects,
  effect_metadata,
  by = c("reporter", "perturbation"),
  all.x = TRUE,
  sort = FALSE
)
effects$reporter_label <- effects$reporter
effects$perturbation_label <- effects$perturbation

decomp <- low_rank_effect_decomposition(
  effects,
  effect = "total_effect",
  reporter = "reporter",
  perturbation = "perturbation",
  rank = analysis_rank,
  impute = "column_mean",
  reporter_label = "reporter_label",
  perturbation_label = "perturbation_label"
)
decomp_for_plot <- decomp
names(decomp)[names(decomp) == "effect"] <- "decomp_total_effect"
effects <- merge(
  effects[, setdiff(names(effects), c("low_rank_effect", "rank_adjusted_effect",
                                      "reporter_score_rank1", "perturbation_score_rank1"))],
  decomp[, c(
    "reporter", "perturbation", "low_rank_effect", "rank_adjusted_effect",
    "reporter_score_rank1", "perturbation_score_rank1"
  )],
  by = c("reporter", "perturbation"),
  all.x = TRUE,
  sort = FALSE
)
effects$rank_adjusted_total_effect <- effects$rank_adjusted_effect
effects$total_effect_se <- effects$total_se
effects$approx_total_pvalue <- effects$total_pvalue
effects$approx_total_padj <- effects$total_padj
effects$approx_relative_pvalue <- effects$relative_pvalue
effects$approx_relative_padj <- effects$relative_padj
effects$approx_rank_adjusted_total_pvalue <- effects$rank_adjusted_total_pvalue
effects$approx_rank_adjusted_total_padj <- effects$rank_adjusted_total_padj

paper_axis <- read_original_svd_axis()
axis_compare <- merge(
  unique(effects[, c(
    "perturbation", "growth_rate", "perturbation_score_rank1", "global_effect"
  )]),
  paper_axis[, c("perturbation", "paper_first_singular_vector")],
  by = "perturbation",
  all.x = TRUE,
  sort = FALSE
)
axis_compare$perturbation_score_rank1_flipped <-
  -axis_compare$perturbation_score_rank1
axis_compare$axis_label <- ""
axis_label_idx <- unique(c(
  head(order(abs(axis_compare$paper_first_singular_vector), decreasing = TRUE), 7),
  head(order(abs(axis_compare$perturbation_score_rank1_flipped), decreasing = TRUE), 7),
  head(order(axis_compare$growth_rate), 3),
  head(order(axis_compare$growth_rate, decreasing = TRUE), 3)
))
axis_compare$axis_label[axis_label_idx] <- axis_compare$perturbation[axis_label_idx]
axis_cor <- data.frame(
  comparison = c(
    "Sign-flipped DStressR rank-1 perturbation score vs paper first singular vector",
    "DStressR global effect vs paper first singular vector",
    "Sign-flipped DStressR rank-1 perturbation score vs growth rate"
  ),
  spearman = c(
    stats::cor(axis_compare$perturbation_score_rank1_flipped, axis_compare$paper_first_singular_vector,
               method = "spearman", use = "complete.obs"),
    stats::cor(axis_compare$global_effect, axis_compare$paper_first_singular_vector,
               method = "spearman", use = "complete.obs"),
    stats::cor(axis_compare$perturbation_score_rank1_flipped, axis_compare$growth_rate,
               method = "spearman", use = "complete.obs")
  ),
  pearson = c(
    stats::cor(axis_compare$perturbation_score_rank1_flipped, axis_compare$paper_first_singular_vector,
               method = "pearson", use = "complete.obs"),
    stats::cor(axis_compare$global_effect, axis_compare$paper_first_singular_vector,
               method = "pearson", use = "complete.obs"),
    stats::cor(axis_compare$perturbation_score_rank1_flipped, axis_compare$growth_rate,
               method = "pearson", use = "complete.obs")
  )
)

paper_global <- read_paper_global_explanation()
promoter_global <- do.call(rbind, lapply(split(effects, effects$reporter), function(dg) {
  finite <- is.finite(dg$total_effect) & is.finite(dg$low_rank_effect)
  if (sum(finite) < 3) {
    r <- NA_real_
  } else {
    r <- stats::cor(dg$total_effect[finite], dg$low_rank_effect[finite])
  }
  data.frame(
    reporter = dg$reporter[1],
    destress_global_correlation = r,
    destress_global_r2 = r^2,
    mean_abs_rank_adjusted_total_effect =
      mean(abs(dg$rank_adjusted_total_effect), na.rm = TRUE),
    max_abs_rank_adjusted_total_effect =
      max(abs(dg$rank_adjusted_total_effect), na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))
promoter_global <- merge(
  promoter_global,
  paper_global,
  by = "reporter",
  all.x = TRUE,
  sort = FALSE
)
promoter_global <- promoter_global[
  order(promoter_global$destress_global_r2, promoter_global$reporter),
  , drop = FALSE
]
promoter_global$rank <- seq_len(nrow(promoter_global))
promoter_global$label <- ""
label_idx <- unique(c(
  head(order(promoter_global$destress_global_r2), 8),
  tail(order(promoter_global$destress_global_r2), 8),
  head(order(-promoter_global$max_abs_rank_adjusted_total_effect), 8)
))
promoter_global$label[label_idx] <- promoter_global$reporter[label_idx]
global_explanation_cor <- data.frame(
  comparison = "DStressR promoter R2 vs Kochanowski promoter R2",
  spearman = stats::cor(
    promoter_global$destress_global_r2,
    promoter_global$paper_global_r2,
    method = "spearman",
    use = "complete.obs"
  ),
  pearson = stats::cor(
    promoter_global$destress_global_r2,
    promoter_global$paper_global_r2,
    method = "pearson",
    use = "complete.obs"
  )
)
condition_diagnostics_total <- perturbation_diagnostics(
  effects,
  mean_effect = "total_effect",
  variance_effect = "total_effect",
  perturbation_label = "perturbation_label",
  trend_span = 0.55
)
condition_diagnostics_total$effect_matrix <- "total"
condition_diagnostics_residual <- perturbation_diagnostics(
  effects,
  mean_effect = "rank_adjusted_total_effect",
  variance_effect = "rank_adjusted_total_effect",
  perturbation_label = "perturbation_label",
  trend_span = 0.55
)
condition_diagnostics_residual$effect_matrix <- "rank_adjusted_total"
condition_diagnostics <- rbind(
  condition_diagnostics_total,
  condition_diagnostics_residual
)
condition_meta <- unique(effects[, c(
  "perturbation", "growth_rate", "perturbation_score_rank1", "global_effect"
)])
condition_diagnostics <- merge(
  condition_diagnostics,
  condition_meta,
  by = "perturbation",
  all.x = TRUE,
  sort = FALSE
)
condition_diagnostics$effect_matrix_label <- ifelse(
  condition_diagnostics$effect_matrix == "total",
  "Total effects",
  "Rank-adjusted total effects"
)
condition_diagnostics$effect_matrix_label <- factor(
  condition_diagnostics$effect_matrix_label,
  levels = c("Total effects", "Rank-adjusted total effects")
)

tail_total <- effect_tail_scores(
  effects,
  effect = "total_effect",
  group = NULL,
  alternative = "two_sided"
)
tail_rank_adjusted <- effect_tail_scores(
  effects,
  effect = "rank_adjusted_total_effect",
  group = NULL,
  alternative = "two_sided"
)
tail_total$tail_analysis <- "Total effect"
tail_rank_adjusted$tail_analysis <- "Rank-adjusted total effect"
tail_scores <- rbind(
  tail_total[, c(names(effects), "diagnostic_center", "diagnostic_scale",
                "diagnostic_z", "tail_probability", "tail_score", "tail_analysis")],
  tail_rank_adjusted[, c(names(effects), "diagnostic_center", "diagnostic_scale",
                         "diagnostic_z", "tail_probability", "tail_score", "tail_analysis")]
)

utils::write.table(
  activity,
  file.path(out_dir, "kochanowski_supplied_promoter_activity_long.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
utils::write.table(
  effects,
  file.path(out_dir, "kochanowski_destress_summary_effects.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
utils::write.table(
  rank_diag,
  file.path(out_dir, "kochanowski_total_effect_rank_diagnostics.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
utils::write.table(
  axis_compare,
  file.path(out_dir, "kochanowski_rank1_axis_comparison.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
utils::write.table(
  axis_cor,
  file.path(out_dir, "kochanowski_rank1_axis_correlations.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
utils::write.table(
  promoter_global,
  file.path(out_dir, "kochanowski_promoter_global_explanation.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
utils::write.table(
  global_explanation_cor,
  file.path(out_dir, "kochanowski_promoter_global_explanation_correlation.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
utils::write.table(
  condition_diagnostics,
  file.path(out_dir, "kochanowski_condition_effect_diagnostics.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
utils::write.table(
  tail_scores,
  file.path(out_dir, "kochanowski_effect_tail_scores.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write_paper_artifact(file.path(out_dir, "kochanowski_paper.rds"),
  results = list(rank0 = effects_rank0, rank1 = effects),
  assays = list(activity = activity), parameters = list(analysis_rank = analysis_rank, selected_rank = selected_rank,
                                                       reference_perturbation = reference_perturbation),
  diagnostics = list(rank_diag = rank_diag, decomp_for_plot = decomp_for_plot,
    axis_compare = axis_compare, axis_cor = axis_cor, promoter_global = promoter_global,
    global_explanation_cor = global_explanation_cor, condition_diagnostics = condition_diagnostics,
    condition_diagnostics_total = condition_diagnostics_total,
    condition_diagnostics_residual = condition_diagnostics_residual,
    tail_scores = tail_scores), annotations = list(condition_labels = condition_labels, conditions = conditions))
