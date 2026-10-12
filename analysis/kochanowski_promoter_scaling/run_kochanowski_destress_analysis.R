#!/usr/bin/env Rscript

# DStressR-style supplied-response analysis of Kochanowski et al. 2017.
#
# The public supplementary tables contain processed steady-state promoter
# activities and estimated standard deviations, not individual replicate wells.
# We therefore analyze total/specific effect matrices and low-rank diagnostics,
# while keeping the p-values below as summary-level approximations.

source(file.path("analysis", "_helpers.R"))

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
effects$approx_specific_pvalue <- effects$specific_pvalue
effects$approx_specific_padj <- effects$specific_padj
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

p_rank <- plot_background_rank_diagnostics(
  diagnostics = rank_diag,
  rank = selected_rank
) +
  theme_koch
ggplot2::ggsave(
  file.path(out_dir, "kochanowski_total_effect_rank_diagnostic.png"),
  p_rank,
  width = 4.8,
  height = 3.2,
  dpi = 300
)
ggplot2::ggsave(
  file.path(out_dir, "kochanowski_total_effect_rank_diagnostic.pdf"),
  p_rank,
  width = 4.8,
  height = 3.2
)

perturbation_order <- conditions$perturbation[conditions$perturbation != reference_perturbation]
p_low_rank <- plot_low_rank_effect_heatmap(
  decomp_for_plot,
  matrices = c("effect", "low_rank_effect", "rank_adjusted_effect"),
  perturbation_order = perturbation_order,
  legend_title = expression(hat(delta)[aj]),
  ylab = "Promoters ordered by rank-1 score"
) +
  theme_koch +
  ggplot2::theme(
    legend.position = "right",
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1, size = 6.5),
    strip.text = ggplot2::element_text(size = 8, face = "bold")
  )
ggplot2::ggsave(
  file.path(out_dir, "kochanowski_total_low_rank_heatmap.png"),
  p_low_rank,
  width = 8.0,
  height = 4.2,
  dpi = 300
)
ggplot2::ggsave(
  file.path(out_dir, "kochanowski_total_low_rank_heatmap.pdf"),
  p_low_rank,
  width = 8.0,
  height = 4.2
)

p_axis <- ggplot2::ggplot(
  axis_compare,
  ggplot2::aes(x = paper_first_singular_vector, y = perturbation_score_rank1_flipped)
) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.25, color = "#808080") +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.25, color = "#808080") +
  ggplot2::geom_point(
    ggplot2::aes(fill = growth_rate),
    shape = 21,
    size = 2.8,
    color = "white",
    stroke = 0.25
  ) +
  ggrepel::geom_text_repel(
    ggplot2::aes(label = axis_label),
    size = 2.5,
    min.segment.length = 0,
    box.padding = 0.2,
    point.padding = 0.1,
    max.overlaps = Inf,
    seed = 1
  ) +
  ggplot2::scale_fill_gradient(
    low = "#9DCDB8",
    high = "#D95F02",
    name = expression("Growth rate"~(h^{-1}))
  ) +
  ggplot2::labs(
    x = "First singular vector reported by Kochanowski et al.",
    y = "Sign-flipped DStressR rank-1 perturbation score"
  ) +
  theme_koch +
  ggplot2::theme(legend.position = "right")
ggplot2::ggsave(
  file.path(out_dir, "kochanowski_rank1_axis_comparison.png"),
  p_axis,
  width = 5.2,
  height = 4.0,
  dpi = 300
)
ggplot2::ggsave(
  file.path(out_dir, "kochanowski_rank1_axis_comparison.pdf"),
  p_axis,
  width = 5.2,
  height = 4.0
)

p_promoter_global <- ggplot2::ggplot(
  promoter_global,
  ggplot2::aes(
    x = destress_global_r2,
    y = paper_global_r2
  )
) +
  ggplot2::geom_abline(
    slope = 1,
    intercept = 0,
    linewidth = 0.3,
    linetype = "longdash",
    color = "#808080"
  ) +
  ggplot2::geom_point(
    ggplot2::aes(fill = max_abs_rank_adjusted_total_effect),
    shape = 21,
    size = 2.2,
    color = "white",
    stroke = 0.25,
    alpha = 0.9
  ) +
  ggrepel::geom_text_repel(
    ggplot2::aes(label = label),
    size = 2.4,
    min.segment.length = 0,
    box.padding = 0.2,
    point.padding = 0.1,
    max.overlaps = Inf,
    seed = 1
  ) +
  ggplot2::scale_fill_gradient(
    low = "#9DCDB8",
    high = "#D95F02",
    name = expression(max~"|"*hat(delta)[aj]^res*"|")
  ) +
  ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = TRUE) +
  ggplot2::labs(
    x = expression("DStressR global explanation " * R[a]^2),
    y = expression("Kochanowski global explanation " * R[a]^2)
  ) +
  theme_koch +
  ggplot2::theme(legend.position = "right")
ggplot2::ggsave(
  file.path(out_dir, "kochanowski_promoter_global_explanation_comparison.png"),
  p_promoter_global,
  width = 5.6,
  height = 4.2,
  dpi = 300
)
ggplot2::ggsave(
  file.path(out_dir, "kochanowski_promoter_global_explanation_comparison.pdf"),
  p_promoter_global,
  width = 5.6,
  height = 4.2
)

p_promoter_global_rank <- ggplot2::ggplot(
  promoter_global,
  ggplot2::aes(x = rank, y = destress_global_r2)
) +
  ggplot2::geom_point(
    ggplot2::aes(fill = max_abs_rank_adjusted_total_effect),
    shape = 21,
    size = 2.1,
    color = "white",
    stroke = 0.25,
    alpha = 0.9
  ) +
  ggrepel::geom_text_repel(
    ggplot2::aes(label = label),
    size = 2.4,
    min.segment.length = 0,
    box.padding = 0.2,
    point.padding = 0.1,
    max.overlaps = Inf,
    seed = 1
  ) +
  ggplot2::scale_fill_gradient(
    low = "#9DCDB8",
    high = "#D95F02",
    name = expression(max~"|"*hat(delta)[aj]^res*"|")
  ) +
  ggplot2::labs(
    x = "Promoters ranked by global explanation",
    y = expression("DStressR global explanation " * R[a]^2)
  ) +
  theme_koch +
  ggplot2::theme(legend.position = "right")
ggplot2::ggsave(
  file.path(out_dir, "kochanowski_promoter_global_explanation_rank.png"),
  p_promoter_global_rank,
  width = 6.2,
  height = 3.6,
  dpi = 300
)
ggplot2::ggsave(
  file.path(out_dir, "kochanowski_promoter_global_explanation_rank.pdf"),
  p_promoter_global_rank,
  width = 6.2,
  height = 3.6
)

condition_label_idx <- unique(c(
  head(order(-condition_diagnostics_total$abs_mean_effect), 5),
  head(order(-condition_diagnostics_total$effect_variance), 5),
  head(order(-condition_diagnostics_residual$effect_variance), 5)
))
condition_diagnostics$label <- ""
label_names <- unique(condition_diagnostics_total$perturbation[condition_label_idx])
condition_diagnostics$label[
  condition_diagnostics$perturbation %in% label_names
] <- condition_diagnostics$perturbation_label[
  condition_diagnostics$perturbation %in% label_names
]
p_condition_rank <- ggplot2::ggplot(
  condition_diagnostics,
  ggplot2::aes(x = rank_abs_mean_effect, y = effect_variance)
) +
  ggplot2::geom_point(
    ggplot2::aes(fill = mean_effect),
    shape = 21,
    size = 2.4,
    color = "white",
    stroke = 0.25,
    alpha = 0.9
  ) +
  ggplot2::geom_line(
    ggplot2::aes(y = variance_trend),
    color = "#008A5B",
    linewidth = 0.55,
    na.rm = TRUE
  ) +
  ggrepel::geom_text_repel(
    ggplot2::aes(label = label),
    size = 2.4,
    min.segment.length = 0,
    box.padding = 0.2,
    point.padding = 0.1,
    max.overlaps = Inf,
    seed = 1
  ) +
  ggplot2::facet_wrap(~effect_matrix_label, nrow = 1, scales = "free_y") +
  ggplot2::scale_fill_gradient2(
    low = "#2166AC",
    mid = "#D1D5DB",
    high = "#D95F02",
    midpoint = 0,
    name = expression(bar(hat(delta))[j])
  ) +
  ggplot2::labs(
    x = expression("Rank of " * "|" * bar(hat(delta))[j] * "|"),
    y = expression("Variance across promoters")
  ) +
  theme_koch +
  ggplot2::theme(
    legend.position = "right",
    strip.text = ggplot2::element_text(face = "bold")
  )
ggplot2::ggsave(
  file.path(out_dir, "kochanowski_condition_rank_variance_diagnostic.png"),
  p_condition_rank,
  width = 7.4,
  height = 3.6,
  dpi = 300
)
ggplot2::ggsave(
  file.path(out_dir, "kochanowski_condition_rank_variance_diagnostic.pdf"),
  p_condition_rank,
  width = 7.4,
  height = 3.6
)

heatmap_reporter_order <- unique(decomp_for_plot[, c("reporter", "reporter_score_rank1")])
heatmap_reporter_order <- heatmap_reporter_order[
  order(heatmap_reporter_order$reporter_score_rank1, heatmap_reporter_order$reporter),
  , drop = FALSE
]
heatmap_reporter_order$reporter_order <- seq_len(nrow(heatmap_reporter_order))
heatmap_limit <- stats::quantile(
  abs(c(
    decomp_for_plot$effect,
    decomp_for_plot$low_rank_effect,
    decomp_for_plot$rank_adjusted_effect
  )),
  probs = 0.985,
  names = FALSE,
  na.rm = TRUE
)
if (!is.finite(heatmap_limit) || heatmap_limit <= 0) {
  heatmap_limit <- max(
    abs(c(
      decomp_for_plot$effect,
      decomp_for_plot$low_rank_effect,
      decomp_for_plot$rank_adjusted_effect
    )),
    na.rm = TRUE
  )
}
if (!is.finite(heatmap_limit) || heatmap_limit <= 0) {
  heatmap_limit <- 1
}

plot_koch_heatmap_panel <- function(effect_col,
                                    panel_tag,
                                    panel_label,
                                    show_y = FALSE,
                                    show_legend = FALSE) {
  d <- decomp_for_plot
  d$reporter_order <- heatmap_reporter_order$reporter_order[
    match(d$reporter, heatmap_reporter_order$reporter)
  ]
  d$perturbation_label <- factor(d$perturbation_label, levels = perturbation_order)
  d$.effect <- pmax(pmin(d[[effect_col]], heatmap_limit), -heatmap_limit)

  ggplot2::ggplot(
    d,
    ggplot2::aes(x = perturbation_label, y = reporter_order, fill = .effect)
  ) +
    ggplot2::geom_tile(width = 0.98, height = 0.95, color = "white", linewidth = 0.05) +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-heatmap_limit, heatmap_limit),
      breaks = c(-heatmap_limit, 0, heatmap_limit),
      labels = signif(c(-heatmap_limit, 0, heatmap_limit), 2),
      name = expression(hat(Delta*y)[aj])
    ) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::scale_x_discrete(expand = c(0, 0)) +
    ggplot2::labs(
      tag = panel_tag,
      title = panel_label,
      x = NULL,
      y = if (show_y) "Promoters ordered by rank-1 score" else NULL
    ) +
    theme_koch +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 8.5, face = "bold", hjust = 0.5),
      plot.tag = ggplot2::element_text(size = 10, face = "bold", color = "#111827"),
      plot.tag.position = c(0.01, 1.01),
      legend.position = if (show_legend) "right" else "none",
      legend.justification = c(0.5, 0.5),
      axis.text.x = ggplot2::element_text(angle = 55, hjust = 1, vjust = 1, size = 5.4),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(6, 4, 3, 4)
    )
}

p_composite_heatmap_total <- plot_koch_heatmap_panel(
  "effect",
  "a",
  "Total effect",
  show_y = TRUE,
  show_legend = FALSE
)
p_composite_heatmap_rank1 <- plot_koch_heatmap_panel(
  "low_rank_effect",
  "b",
  "Rank-1 component",
  show_y = FALSE,
  show_legend = FALSE
)
p_composite_heatmap_residual <- plot_koch_heatmap_panel(
  "rank_adjusted_effect",
  "c",
  "Rank-adjusted total effect",
  show_y = FALSE,
  show_legend = TRUE
)

p_composite_rank <- p_condition_rank +
  ggplot2::labs(tag = "d") +
  ggplot2::theme(
    plot.tag = ggplot2::element_text(size = 10, face = "bold", color = "#111827"),
    plot.tag.position = c(0.01, 1.01),
    legend.position = "right",
    plot.margin = ggplot2::margin(6, 6, 4, 4)
  )
p_composite_axis <- p_axis +
  ggplot2::labs(tag = "e") +
  ggplot2::theme(
    plot.tag = ggplot2::element_text(size = 10, face = "bold", color = "#111827"),
    plot.tag.position = c(0.01, 1.01),
    legend.position = "right",
    plot.margin = ggplot2::margin(6, 6, 4, 4)
  )

koch_composite <- gridExtra::arrangeGrob(
  gridExtra::arrangeGrob(
    p_composite_heatmap_total,
    p_composite_heatmap_rank1,
    p_composite_heatmap_residual,
    ncol = 3,
    widths = c(1, 1, 1.18)
  ),
  gridExtra::arrangeGrob(
    p_composite_rank,
    p_composite_axis,
    ncol = 2,
    widths = c(1.18, 1)
  ),
  ncol = 1,
  heights = c(1.05, 1)
)
ggplot2::ggsave(
  file.path(out_dir, "kochanowski_dstressr_composite_figure.png"),
  koch_composite,
  width = 10.8,
  height = 7.6,
  dpi = 300,
  bg = "white"
)
ggplot2::ggsave(
  file.path(out_dir, "kochanowski_dstressr_composite_figure.pdf"),
  koch_composite,
  width = 10.8,
  height = 7.6,
  bg = "white"
)

condition_order_residual <- condition_diagnostics_residual$perturbation[
  order(condition_diagnostics_residual$mean_effect)
]
residual_distribution <- effects
residual_distribution$perturbation_label <- factor(
  residual_distribution$perturbation_label,
  levels = condition_diagnostics_residual$perturbation_label[
    match(condition_order_residual, condition_diagnostics_residual$perturbation)
  ]
)
p_condition_residual <- ggplot2::ggplot(
  residual_distribution,
  ggplot2::aes(x = perturbation_label, y = rank_adjusted_total_effect)
) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.25, color = "#808080") +
  ggplot2::geom_boxplot(
    width = 0.58,
    outlier.shape = NA,
    fill = "#E5E7EB",
    color = "#4B5563",
    linewidth = 0.28
  ) +
  ggplot2::geom_jitter(
    ggplot2::aes(color = rank_adjusted_total_effect),
    width = 0.14,
    height = 0,
    size = 0.75,
    alpha = 0.65
  ) +
  ggplot2::scale_color_gradient2(
    low = "#2166AC",
    mid = "#9CA3AF",
    high = "#D95F02",
    midpoint = 0,
    guide = "none"
  ) +
  ggplot2::coord_flip() +
  ggplot2::labs(
    x = NULL,
    y = expression("Rank-adjusted total effect " * hat(delta)[aj]^res)
  ) +
  theme_koch +
  ggplot2::theme(
    axis.text.y = ggplot2::element_text(size = 7),
    legend.position = "none"
  )
ggplot2::ggsave(
  file.path(out_dir, "kochanowski_condition_residual_distributions.png"),
  p_condition_residual,
  width = 5.8,
  height = 6.0,
  dpi = 300
)
ggplot2::ggsave(
  file.path(out_dir, "kochanowski_condition_residual_distributions.pdf"),
  p_condition_residual,
  width = 5.8,
  height = 6.0
)

p_hist <- ggplot2::ggplot(
  effects,
  ggplot2::aes(x = approx_rank_adjusted_total_pvalue)
) +
  ggplot2::geom_histogram(
    bins = 40,
    boundary = 0,
    fill = "#4F8F7B",
    color = "white",
    linewidth = 0.15
  ) +
  ggplot2::labs(
    x = "Approximate p-value",
    y = "Reporter-condition pairs"
  ) +
  theme_koch
ggplot2::ggsave(
  file.path(out_dir, "kochanowski_approx_rank_adjusted_pvalue_histogram.png"),
  p_hist,
  width = 4.8,
  height = 3.2,
  dpi = 300
)
ggplot2::ggsave(
  file.path(out_dir, "kochanowski_approx_rank_adjusted_pvalue_histogram.pdf"),
  p_hist,
  width = 4.8,
  height = 3.2
)

volcano_data <- effects
volcano_data$label <- ""
label_idx <- order(volcano_data$approx_rank_adjusted_total_padj, -abs(volcano_data$rank_adjusted_total_effect))
label_idx <- label_idx[is.finite(volcano_data$approx_rank_adjusted_total_padj[label_idx])]
label_idx <- head(label_idx, 15)
volcano_data$label[label_idx] <- paste(volcano_data$reporter[label_idx], volcano_data$perturbation[label_idx], sep = " / ")
volcano_data$minus_log10_p <- -log10(pmax(volcano_data$approx_rank_adjusted_total_pvalue, .Machine$double.xmin))
p_volcano <- ggplot2::ggplot(
  volcano_data,
  ggplot2::aes(x = rank_adjusted_total_effect, y = minus_log10_p)
) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.25, color = "#808080") +
  ggplot2::geom_point(
    ggplot2::aes(color = approx_rank_adjusted_total_padj < 0.05),
    alpha = 0.75,
    size = 1.4
  ) +
  ggplot2::scale_color_manual(
    values = c("FALSE" = "#A8A8A8", "TRUE" = "#D95F02"),
    guide = "none"
  ) +
  ggrepel::geom_text_repel(
    ggplot2::aes(label = label),
    size = 2.4,
    min.segment.length = 0,
    box.padding = 0.2,
    point.padding = 0.08,
    max.overlaps = Inf,
    seed = 1
  ) +
  ggplot2::labs(
    x = expression("Rank-adjusted total effect " * hat(delta)[aj]^{"res"}),
    y = expression(-log[10]("approximate p-value"))
  ) +
  theme_koch
ggplot2::ggsave(
  file.path(out_dir, "kochanowski_approx_rank_adjusted_volcano.png"),
  p_volcano,
  width = 5.4,
  height = 3.8,
  dpi = 300
)
ggplot2::ggsave(
  file.path(out_dir, "kochanowski_approx_rank_adjusted_volcano.pdf"),
  p_volcano,
  width = 5.4,
  height = 3.8
)

plot_koch_volcano_panel <- function(effect_col,
                                    pvalue_col,
                                    padj_col,
                                    panel_tag,
                                    panel_label,
                                    xlab,
                                    show_y = FALSE,
                                    show_legend = FALSE) {
  plot_volcano(
    effects,
    effect = effect_col,
    pvalue = pvalue_col,
    padj = padj_col,
    y = "pvalue",
    reporter = "reporter",
    perturbation = "perturbation",
    perturbation_label = "perturbation_label",
    fdr = 0.05,
    top_n = 10,
    color_by = "significance",
    title = NULL,
    subtitle = "",
    xlab = xlab,
    ylab = if (show_y) expression(-log[10]("approximate p-value")) else NULL,
    point_size = 0.95,
    point_alpha = 0.72,
    label_by = "pair",
    max_label_chars = 34,
    call_labels = c("Not called", "FDR < 0.05")
  ) +
    ggplot2::labs(
      tag = panel_tag,
      title = panel_label,
      color = "Approximate call"
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 8.5, face = "bold", hjust = 0.5),
      plot.tag = ggplot2::element_text(size = 10, face = "bold", color = "#111827"),
      plot.tag.position = c(0.01, 1.01),
      legend.position = if (show_legend) "right" else "none",
      legend.justification = c(0.5, 0.5),
      axis.title.y = if (show_y) ggplot2::element_text() else ggplot2::element_blank(),
      axis.text.y = if (show_y) ggplot2::element_text() else ggplot2::element_blank(),
      axis.ticks.y = if (show_y) ggplot2::element_line() else ggplot2::element_blank(),
      plot.margin = ggplot2::margin(6, 4, 4, 4)
    )
}

p_volcano_total <- plot_koch_volcano_panel(
  effect_col = "total_effect",
  pvalue_col = "approx_total_pvalue",
  padj_col = "approx_total_padj",
  panel_tag = "a",
  panel_label = "Total effect",
  xlab = expression(hat(Delta*y)[aj]^"tot"),
  show_y = TRUE,
  show_legend = FALSE
)
p_volcano_specific <- plot_koch_volcano_panel(
  effect_col = "specific_effect",
  pvalue_col = "approx_specific_pvalue",
  padj_col = "approx_specific_padj",
  panel_tag = "b",
  panel_label = "Specific effect",
  xlab = expression(hat(Delta*y)[aj]^"spec"),
  show_y = FALSE,
  show_legend = FALSE
)
p_volcano_rank_adjusted <- plot_koch_volcano_panel(
  effect_col = "rank_adjusted_total_effect",
  pvalue_col = "approx_rank_adjusted_total_pvalue",
  padj_col = "approx_rank_adjusted_total_padj",
  panel_tag = "c",
  panel_label = "Rank-adjusted total effect",
  xlab = expression(hat(Delta*y)[aj]^"tot,k"),
  show_y = FALSE,
  show_legend = FALSE
)
koch_volcano_composite <- gridExtra::arrangeGrob(
  p_volcano_total,
  p_volcano_specific,
  p_volcano_rank_adjusted,
  ncol = 3,
  widths = c(1, 1, 1)
)
ggplot2::ggsave(
  file.path(out_dir, "kochanowski_approx_volcano_three_targets.png"),
  koch_volcano_composite,
  width = 10.2,
  height = 3.8,
  dpi = 300,
  bg = "white"
)
ggplot2::ggsave(
  file.path(out_dir, "kochanowski_approx_volcano_three_targets.pdf"),
  koch_volcano_composite,
  width = 10.2,
  height = 3.8,
  bg = "white"
)

p_tail <- plot_effect_tail_histogram(
  tail_scores,
  effect = "rank_adjusted_total_effect",
  tail_probability = "tail_probability",
  facet = "tail_analysis",
  bins = 35,
  xlab = expression("Effect " * hat(delta)[aj]),
  ylab = "Reporter-condition pairs"
) +
  theme_koch
ggplot2::ggsave(
  file.path(out_dir, "kochanowski_effect_tail_histograms.png"),
  p_tail,
  width = 6.2,
  height = 3.2,
  dpi = 300
)
ggplot2::ggsave(
  file.path(out_dir, "kochanowski_effect_tail_histograms.pdf"),
  p_tail,
  width = 6.2,
  height = 3.2
)

summary <- data.frame(
  n_reporters = length(unique(effects$reporter)),
  n_perturbations = length(unique(effects$perturbation)),
  reference_perturbation = reference_perturbation,
  selected_rank = selected_rank,
  analysis_rank = analysis_rank,
  finite_pairs = sum(is.finite(effects$total_effect)),
  approx_total_fdr_0_05 = sum(effects$approx_total_padj < 0.05, na.rm = TRUE),
  approx_specific_fdr_0_05 = sum(effects$approx_specific_padj < 0.05, na.rm = TRUE),
  approx_rank_adjusted_total_fdr_0_05 =
    sum(effects$approx_rank_adjusted_total_padj < 0.05, na.rm = TRUE)
)
utils::write.table(
  summary,
  file.path(out_dir, "kochanowski_analysis_summary.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

message("Wrote Kochanowski supplied-response outputs to: ", out_dir)
message("Reporters: ", summary$n_reporters)
message("Perturbations vs reference: ", summary$n_perturbations)
message("Selected background rank by diagnostic: ", selected_rank)
message("Reported low-rank analysis rank: ", analysis_rank)
message(
  "Approximate FDR < 0.05 after rank adjustment: ",
  summary$approx_rank_adjusted_total_fdr_0_05
)
