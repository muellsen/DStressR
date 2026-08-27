# Exploratory DStressR-style analysis of Zaslaver et al. promoter activities.
#
# This script deliberately stays outside the package vignettes/manuscript. The
# fixed-growth-rate table from Zaslaver et al. is an author-derived summary with
# one value per reporter, condition, and growth rate. We therefore estimate
# effect matrices and diagnostics, but we do not treat these summaries as
# replicated well-level observations for formal DStressR hit calling.

stopifnot(requireNamespace("ggplot2", quietly = TRUE))
stopifnot(requireNamespace("pkgload", quietly = TRUE))
pkgload::load_all(".", quiet = TRUE)

out_dir <- file.path("analysis", "outputs", "zaslaver_promoter_activity")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

load(file.path("data", "zaslaver_promoter_activity.rda"))
load(file.path("data", "zaslaver_promoter_annotations.rda"))
load(file.path("data", "zaslaver_promoter_timecourse.rda"))

eps <- 1
reference_condition <- "Glucose"
top_n_reporters <- 60

make_response_effects <- function(activity, reference = "Glucose", eps = 1) {
  d <- activity
  d$response <- log2(d$promoter_activity + eps)
  ref <- d[d$condition_label == reference, c(
    "reporter_id", "growth_rate", "response", "promoter_activity", "reporter"
  )]
  names(ref) <- c(
    "reporter_id", "growth_rate", "reference_response",
    "reference_promoter_activity", "reference_reporter"
  )
  d <- merge(d, ref, by = c("reporter_id", "growth_rate"), all.x = TRUE, sort = FALSE)
  d <- d[d$condition_label != reference, , drop = FALSE]
  d$total_effect <- d$response - d$reference_response
  global <- stats::aggregate(
    total_effect ~ growth_rate + condition_label,
    d,
    mean,
    na.rm = TRUE
  )
  names(global)[3] <- "global_effect"
  d <- merge(d, global, by = c("growth_rate", "condition_label"), all.x = TRUE, sort = FALSE)
  d$specific_effect <- d$total_effect - d$global_effect
  d
}

effect_table <- make_response_effects(
  zaslaver_promoter_activity,
  reference = reference_condition,
  eps = eps
)
annotation_focus <- c(
  "reporter_id", "ribosomes", "ribosomal_proteins", "metabolism",
  "carbon_utilization", "nitrogen_metabolism", "phosphorous_metabolism",
  "other_.mechanical._nutritional._oxidative_stress."
)
annotation_focus <- intersect(annotation_focus, names(zaslaver_promoter_annotations))

effect_table <- merge(
  effect_table,
  zaslaver_promoter_annotations[, annotation_focus, drop = FALSE],
  by = "reporter_id",
  all.x = TRUE,
  sort = FALSE
)

utils::write.table(
  effect_table,
  file.path(out_dir, "zaslaver_fixed_growth_effects.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

global_effects <- unique(effect_table[, c("growth_rate", "condition_label", "global_effect")])
global_effects <- global_effects[order(global_effects$growth_rate, global_effects$condition_label), ]
utils::write.table(
  global_effects,
  file.path(out_dir, "zaslaver_condition_global_effects.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

top_specific <- effect_table[
  order(-abs(effect_table$specific_effect)),
  c("growth_rate", "condition_label", "reporter_id", "reporter", "reference_reporter",
    "promoter_activity", "reference_promoter_activity", "total_effect",
    "global_effect", "specific_effect", "ribosomes", "ribosomal_proteins",
    "metabolism", "carbon_utilization", "nitrogen_metabolism", "phosphorous_metabolism")
]
top_specific <- head(top_specific, 300)
utils::write.table(
  top_specific,
  file.path(out_dir, "zaslaver_top_specific_effects.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

matrix_from_effects <- function(d, growth_rate, effect) {
  dg <- d[d$growth_rate == growth_rate, , drop = FALSE]
  reporters <- unique(dg$reporter_id)
  conditions <- unique(dg$condition_label)
  mat <- matrix(
    NA_real_,
    nrow = length(reporters),
    ncol = length(conditions),
    dimnames = list(reporters, conditions)
  )
  idx <- cbind(match(dg$reporter_id, reporters), match(dg$condition_label, conditions))
  mat[idx] <- dg[[effect]]
  mat
}

diag_list <- lapply(sort(unique(effect_table$growth_rate)), function(gr) {
  dg <- effect_table[effect_table$growth_rate == gr, , drop = FALSE]
  d <- background_rank_diagnostics(
    dg,
    effect = "total_effect",
    reporter = "reporter_id",
    perturbation = "condition_label",
    rank_max = 5,
    permutations = 200,
    seed = 100 + as.integer(gr * 100)
  )
  d$growth_rate <- gr
  d$group_label <- paste0("g = ", gr)
  d
})
rank_diag <- do.call(rbind, diag_list)
utils::write.table(
  rank_diag,
  file.path(out_dir, "zaslaver_total_effect_rank_diagnostics.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

annotation_cols <- setdiff(names(zaslaver_promoter_annotations), c("reporter_index", "reporter_id", "reporter"))
top_by_reporter <- stats::aggregate(
  abs(specific_effect) ~ growth_rate + reporter_id,
  effect_table,
  max,
  na.rm = TRUE
)
names(top_by_reporter)[3] <- "max_abs_specific_effect"

enrichment_rows <- lapply(split(top_by_reporter, top_by_reporter$growth_rate), function(dg) {
  cutoff <- stats::quantile(dg$max_abs_specific_effect, probs = 0.95, na.rm = TRUE)
  selected <- dg$reporter_id[dg$max_abs_specific_effect >= cutoff]
  ann <- zaslaver_promoter_annotations
  ann$selected <- ann$reporter_id %in% selected
  rows <- lapply(annotation_cols, function(term) {
    x <- ann[[term]]
    if (!all(x %in% c(0L, 1L, NA_integer_), na.rm = TRUE)) {
      return(NULL)
    }
    tab <- table(factor(ann$selected, c(FALSE, TRUE)), factor(x == 1L, c(FALSE, TRUE)))
    if (any(dim(tab) != c(2, 2)) || sum(tab[, 2]) < 5) {
      return(NULL)
    }
    ft <- stats::fisher.test(tab, alternative = "greater")
    data.frame(
      growth_rate = unique(dg$growth_rate),
      term = term,
      selected_with_term = tab[2, 2],
      selected_without_term = tab[2, 1],
      background_with_term = tab[1, 2],
      background_without_term = tab[1, 1],
      odds_ratio = unname(ft$estimate),
      pvalue = ft$p.value,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
})
enrichment <- do.call(rbind, enrichment_rows)
enrichment$padj <- stats::p.adjust(enrichment$pvalue, method = "BH")
enrichment <- enrichment[order(enrichment$growth_rate, enrichment$padj, -enrichment$odds_ratio), ]
utils::write.table(
  enrichment,
  file.path(out_dir, "zaslaver_top_specific_annotation_enrichment.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

plot_theme <- function(base_size = 9) {
  ggplot2::theme_light(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_blank(),
      legend.position = "right",
      strip.background = ggplot2::element_rect(fill = "#F5F5F5", color = "#D0D0D0"),
      strip.text = ggplot2::element_text(face = "bold", color = "#333333")
    )
}

p_global <- ggplot2::ggplot(
  global_effects,
  ggplot2::aes(x = condition_label, y = global_effect, fill = factor(growth_rate))
) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.25, color = "#606060") +
  ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.72), width = 0.64) +
  ggplot2::scale_fill_manual(values = c("0.25" = "#76B7B2", "0.8" = "#E15759"), name = "Growth rate") +
  ggplot2::labs(x = NULL, y = expression("Global condition effect " * hat(delta)[j])) +
  ggplot2::coord_flip() +
  plot_theme(9)

ggplot2::ggsave(
  file.path(out_dir, "zaslaver_global_condition_effects.png"),
  p_global,
  width = 6.0,
  height = 3.5,
  dpi = 300
)
ggplot2::ggsave(
  file.path(out_dir, "zaslaver_global_condition_effects.pdf"),
  p_global,
  width = 6.0,
  height = 3.5
)

top_ids <- unique(top_specific$reporter_id)[seq_len(min(top_n_reporters, length(unique(top_specific$reporter_id))))]
heat_data <- effect_table[effect_table$reporter_id %in% top_ids, , drop = FALSE]
label_map <- unique(heat_data[, c("reporter_id", "reporter")])
label_map$label <- paste0(label_map$reporter, " (", sub("^reporter_", "", label_map$reporter_id), ")")
label_map <- label_map[!duplicated(label_map$reporter_id), ]
heat_data$reporter_label <- label_map$label[match(heat_data$reporter_id, label_map$reporter_id)]
heat_data$reporter_label <- factor(heat_data$reporter_label, levels = rev(label_map$label))

p_heat <- ggplot2::ggplot(
  heat_data,
  ggplot2::aes(x = condition_label, y = reporter_label, fill = specific_effect)
) +
  ggplot2::geom_tile(color = "white", linewidth = 0.12) +
  ggplot2::facet_wrap(~growth_rate, nrow = 1, labeller = ggplot2::label_both) +
  ggplot2::scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    name = expression(hat(delta)[aj]^spec)
  ) +
  ggplot2::labs(x = NULL, y = NULL) +
  plot_theme(7) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
    axis.text.y = ggplot2::element_text(size = 5.4),
    legend.key.height = grid::unit(12, "mm")
  )

ggplot2::ggsave(
  file.path(out_dir, "zaslaver_top_specific_effect_heatmap.png"),
  p_heat,
  width = 8.0,
  height = 8.5,
  dpi = 300
)
ggplot2::ggsave(
  file.path(out_dir, "zaslaver_top_specific_effect_heatmap.pdf"),
  p_heat,
  width = 8.0,
  height = 8.5
)

p_rank <- plot_background_rank_diagnostics(
  diagnostics = rank_diag,
  xlab = "Component",
  ylab = "Singular value"
)

ggplot2::ggsave(
  file.path(out_dir, "zaslaver_total_effect_rank_diagnostics.png"),
  p_rank,
  width = 5.4,
  height = 3.8,
  dpi = 300
)
ggplot2::ggsave(
  file.path(out_dir, "zaslaver_total_effect_rank_diagnostics.pdf"),
  p_rank,
  width = 5.4,
  height = 3.8
)

rank1_table <- low_rank_effect_decomposition(
  effect_table,
  effect = "total_effect",
  reporter = "reporter_id",
  perturbation = "condition_label",
  group = "growth_rate",
  rank = 1,
  reporter_label = "reporter",
  perturbation_label = "condition_label"
)
names(rank1_table)[names(rank1_table) == "reporter"] <- "reporter_id"
names(rank1_table)[names(rank1_table) == "perturbation"] <- "condition_label"
names(rank1_table)[names(rank1_table) == "reporter_label"] <- "reporter"
rank1_table$total_effect <- rank1_table$effect
rank1_table$rank1_effect <- rank1_table$low_rank_effect
rank1_table$residual_effect <- rank1_table$rank_adjusted_effect
rank1_table$rank1_row_score <- rank1_table$reporter_score_rank1
rank1_order <- unique(rank1_table[, c("growth_rate", "reporter_id", "rank1_row_score")])
rank1_order <- rank1_order[order(rank1_order$growth_rate, rank1_order$rank1_row_score, rank1_order$reporter_id), ]
rank1_order$row_order <- ave(seq_len(nrow(rank1_order)), rank1_order$growth_rate, FUN = seq_along)
rank1_table$row_order <- rank1_order$row_order[
  match(paste(rank1_table$growth_rate, rank1_table$reporter_id), paste(rank1_order$growth_rate, rank1_order$reporter_id))
]

rank1_table <- merge(
  rank1_table,
  unique(effect_table[, c("growth_rate", "condition_label", "reporter_id", "reporter", "specific_effect")]),
  by = c("growth_rate", "condition_label", "reporter_id"),
  all.x = TRUE,
  sort = FALSE
)
if ("reporter.x" %in% names(rank1_table)) {
  rank1_table$reporter <- rank1_table$reporter.x
  missing_reporter <- is.na(rank1_table$reporter) | !nzchar(rank1_table$reporter)
  if ("reporter.y" %in% names(rank1_table)) {
    rank1_table$reporter[missing_reporter] <- rank1_table$reporter.y[missing_reporter]
  }
  rank1_table$reporter.x <- NULL
  rank1_table$reporter.y <- NULL
}
names(rank1_table)[names(rank1_table) == "specific_effect"] <- "specific_effect_rank0"
rank1_table$rank_adjusted_total_effect <- rank1_table$residual_effect
rank1_global <- stats::aggregate(
  rank_adjusted_total_effect ~ growth_rate + condition_label,
  rank1_table,
  mean,
  na.rm = TRUE
)
names(rank1_global)[3] <- "rank_adjusted_global_effect"
rank1_table <- merge(
  rank1_table,
  rank1_global,
  by = c("growth_rate", "condition_label"),
  all.x = TRUE,
  sort = FALSE
)
rank1_table$specific_effect_rank1 <- rank1_table$rank_adjusted_total_effect -
  rank1_table$rank_adjusted_global_effect

rank0_tail_scores <- effect_tail_scores(
  rank1_table,
  effect = "specific_effect_rank0",
  group = c("growth_rate", "condition_label")
)
rank1_tail_scores <- effect_tail_scores(
  rank1_table,
  effect = "specific_effect_rank1",
  group = c("growth_rate", "condition_label")
)
rank1_table$diagnostic_pvalue_rank0 <- rank0_tail_scores$tail_probability
rank1_table$diagnostic_pvalue_rank1 <- rank1_tail_scores$tail_probability
rank1_table$diagnostic_padj_global_rank0 <- stats::p.adjust(rank1_table$diagnostic_pvalue_rank0, method = "BH")
rank1_table$diagnostic_padj_global_rank1 <- stats::p.adjust(rank1_table$diagnostic_pvalue_rank1, method = "BH")

utils::write.table(
  rank1_table,
  file.path(out_dir, "zaslaver_total_effect_rank1_matrix.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

rank1_plot_table <- rank1_table
rank1_plot_table$reporter <- rank1_plot_table$reporter_id
rank1_plot_table$perturbation <- rank1_plot_table$condition_label
rank1_plot_table$reporter_label <- rank1_plot_table$reporter
rank1_plot_table$perturbation_label <- rank1_plot_table$condition_label
p_destress_rank1_heatmap <- plot_low_rank_effect_heatmap(
  rank1_plot_table,
  matrices = c("effect", "low_rank_effect", "rank_adjusted_effect"),
  perturbation_order = c("Ethanol", "Nitrogen limited", "no AA", "no Glucose", "Phosphate limited"),
  xlab = NULL,
  ylab = "Reporters ordered by rank-1 score",
  legend_title = "Effect"
)

ggplot2::ggsave(
  file.path(out_dir, "zaslaver_destress_low_rank_effect_heatmap.png"),
  p_destress_rank1_heatmap,
  width = 8.6,
  height = 5.6,
  dpi = 300
)
ggplot2::ggsave(
  file.path(out_dir, "zaslaver_destress_low_rank_effect_heatmap.pdf"),
  p_destress_rank1_heatmap,
  width = 8.6,
  height = 5.6
)

rank1_pvalue_long <- rbind(
  data.frame(
    rank1_table[, c("growth_rate", "condition_label", "reporter_id", "reporter")],
    model = "before rank-1 adjustment",
    effect = rank1_table$specific_effect_rank0,
    diagnostic_pvalue = rank1_table$diagnostic_pvalue_rank0,
    diagnostic_padj_global = rank1_table$diagnostic_padj_global_rank0
  ),
  data.frame(
    rank1_table[, c("growth_rate", "condition_label", "reporter_id", "reporter")],
    model = "after rank-1 adjustment",
    effect = rank1_table$specific_effect_rank1,
    diagnostic_pvalue = rank1_table$diagnostic_pvalue_rank1,
    diagnostic_padj_global = rank1_table$diagnostic_padj_global_rank1
  )
)
rank1_pvalue_long$model <- factor(
  rank1_pvalue_long$model,
  levels = c("before rank-1 adjustment", "after rank-1 adjustment")
)
rank1_pvalue_long$growth_label <- paste0("g = ", rank1_pvalue_long$growth_rate)
utils::write.table(
  rank1_pvalue_long,
  file.path(out_dir, "zaslaver_fixed_growth_rank1_diagnostic_pvalues.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

p_destress_tail_hist <- plot_effect_tail_histogram(
  rank1_pvalue_long,
  effect = "effect",
  tail_probability = "diagnostic_pvalue",
  facet = c("model", "growth_label", "condition_label"),
  tail_threshold = 0.05,
  bins = 40,
  xlab = "Effect",
  ylab = "Promoter-condition effects"
)

ggplot2::ggsave(
  file.path(out_dir, "zaslaver_destress_effect_tail_histograms.png"),
  p_destress_tail_hist,
  width = 8.2,
  height = 6.8,
  dpi = 300
)
ggplot2::ggsave(
  file.path(out_dir, "zaslaver_destress_effect_tail_histograms.pdf"),
  p_destress_tail_hist,
  width = 8.2,
  height = 6.8
)

p_fixed_rank_hist <- ggplot2::ggplot(
  rank1_pvalue_long,
  ggplot2::aes(x = diagnostic_pvalue)
) +
  ggplot2::geom_histogram(bins = 50, fill = "#6BAED6", color = "white", linewidth = 0.15) +
  ggplot2::facet_grid(model + growth_label ~ condition_label) +
  ggplot2::labs(x = "Diagnostic tail probability", y = "Promoter-condition effects") +
  plot_theme(7) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

ggplot2::ggsave(
  file.path(out_dir, "zaslaver_fixed_growth_rank1_pvalue_histograms.png"),
  p_fixed_rank_hist,
  width = 8.2,
  height = 6.8,
  dpi = 300
)
ggplot2::ggsave(
  file.path(out_dir, "zaslaver_fixed_growth_rank1_pvalue_histograms.pdf"),
  p_fixed_rank_hist,
  width = 8.2,
  height = 6.8
)

volcano_fixed <- rank1_pvalue_long[
  rank1_pvalue_long$model == "after rank-1 adjustment" &
    is.finite(rank1_pvalue_long$diagnostic_pvalue) &
    rank1_pvalue_long$diagnostic_pvalue > 0,
]
volcano_fixed$neg_log10_p <- -log10(volcano_fixed$diagnostic_pvalue)
volcano_fixed$hit <- volcano_fixed$diagnostic_padj_global < 0.05
volcano_fixed$growth_label <- paste0("g = ", volcano_fixed$growth_rate)
label_fixed <- do.call(rbind, lapply(split(volcano_fixed, paste(volcano_fixed$growth_rate, volcano_fixed$condition_label)), function(d) {
  head(d[order(d$diagnostic_pvalue, -abs(d$effect)), ], 2)
}))
label_fixed$label <- paste(label_fixed$reporter, label_fixed$condition_label, sep = " / ")

p_fixed_volcano <- ggplot2::ggplot(
  volcano_fixed,
  ggplot2::aes(x = effect, y = neg_log10_p)
) +
  ggplot2::geom_point(
    ggplot2::aes(color = hit),
    size = 0.65,
    alpha = 0.6
  ) +
  ggplot2::scale_color_manual(values = c("FALSE" = "#BDBDBD", "TRUE" = "#D55E00"), guide = "none") +
  ggplot2::facet_grid(growth_label ~ condition_label) +
  ggplot2::labs(
    x = expression("Rank-1 adjusted specific effect " * hat(delta)[aj]^spec),
    y = expression(-log[10] * " diagnostic tail probability")
  ) +
  plot_theme(7)

if (requireNamespace("ggrepel", quietly = TRUE)) {
  p_fixed_volcano <- p_fixed_volcano +
    ggrepel::geom_text_repel(
      data = label_fixed,
      ggplot2::aes(label = label),
      size = 1.8,
      min.segment.length = 0,
      box.padding = 0.18,
      max.overlaps = Inf
    )
}

ggplot2::ggsave(
  file.path(out_dir, "zaslaver_fixed_growth_rank1_volcano.png"),
  p_fixed_volcano,
  width = 8.8,
  height = 5.4,
  dpi = 300
)
ggplot2::ggsave(
  file.path(out_dir, "zaslaver_fixed_growth_rank1_volcano.pdf"),
  p_fixed_volcano,
  width = 8.8,
  height = 5.4
)

rank1_score_wide <- stats::reshape(
  unique(rank1_table[, c("reporter_id", "reporter", "growth_rate", "rank1_row_score")]),
  idvar = c("reporter_id", "reporter"),
  timevar = "growth_rate",
  direction = "wide"
)
names(rank1_score_wide) <- sub("rank1_row_score\\.0\\.25", "rank1_score_g0_25", names(rank1_score_wide))
names(rank1_score_wide) <- sub("rank1_row_score\\.0\\.8", "rank1_score_g0_8", names(rank1_score_wide))
embedding_annotation_terms <- c(
  "reporter_id", "ribosomes", "translation", "ribosomal_proteins",
  "amino_acids", "biosynthesis_of_building_blocks", "metabolism",
  "membrane", "transport", "SOS_response",
  "other_.mechanical._nutritional._oxidative_stress."
)
embedding_annotation_terms <- intersect(embedding_annotation_terms, names(zaslaver_promoter_annotations))
rank1_embedding <- merge(
  rank1_score_wide,
  zaslaver_promoter_annotations[, embedding_annotation_terms, drop = FALSE],
  by = "reporter_id",
  all.x = TRUE,
  sort = FALSE
)
rank1_embedding$physiology <- "other"
if ("metabolism" %in% names(rank1_embedding)) {
  rank1_embedding$physiology[rank1_embedding$metabolism == 1L] <- "metabolism"
}
if ("membrane" %in% names(rank1_embedding)) {
  rank1_embedding$physiology[rank1_embedding$membrane == 1L] <- "membrane"
}
if ("transport" %in% names(rank1_embedding)) {
  rank1_embedding$physiology[rank1_embedding$transport == 1L] <- "transport"
}
if ("amino_acids" %in% names(rank1_embedding)) {
  rank1_embedding$physiology[rank1_embedding$amino_acids == 1L] <- "amino acids"
}
if ("biosynthesis_of_building_blocks" %in% names(rank1_embedding)) {
  rank1_embedding$physiology[rank1_embedding$biosynthesis_of_building_blocks == 1L] <- "biosynthesis"
}
if ("SOS_response" %in% names(rank1_embedding)) {
  rank1_embedding$physiology[rank1_embedding$SOS_response == 1L] <- "SOS response"
}
if ("translation" %in% names(rank1_embedding)) {
  rank1_embedding$physiology[rank1_embedding$translation == 1L] <- "translation"
}
if ("ribosomes" %in% names(rank1_embedding)) {
  rank1_embedding$physiology[rank1_embedding$ribosomes == 1L] <- "ribosomes"
}
rank1_embedding$physiology <- factor(
  rank1_embedding$physiology,
  levels = c("ribosomes", "translation", "biosynthesis", "amino acids",
             "metabolism", "transport", "membrane", "SOS response", "other")
)
utils::write.table(
  rank1_embedding,
  file.path(out_dir, "zaslaver_rank1_growth_rate_embedding.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

embedding_colors <- c(
  "ribosomes" = "#D55E00",
  "translation" = "#E69F00",
  "biosynthesis" = "#009E73",
  "amino acids" = "#56B4E9",
  "metabolism" = "#0072B2",
  "transport" = "#CC79A7",
  "membrane" = "#7F7F7F",
  "SOS response" = "#B2182B",
  "other" = "#D0D0D0"
)

p_rank_embedding <- ggplot2::ggplot(
  rank1_embedding,
  ggplot2::aes(x = rank1_score_g0_25, y = rank1_score_g0_8)
) +
  ggplot2::geom_hline(yintercept = 0, color = "#707070", linewidth = 0.25) +
  ggplot2::geom_vline(xintercept = 0, color = "#707070", linewidth = 0.25) +
  ggplot2::geom_point(ggplot2::aes(color = physiology), size = 1.0, alpha = 0.75) +
  ggplot2::scale_color_manual(values = embedding_colors, drop = FALSE, name = NULL) +
  ggplot2::labs(
    x = expression("Rank-1 reporter score at " * g == 0.25),
    y = expression("Rank-1 reporter score at " * g == 0.8)
  ) +
  plot_theme(8) +
  ggplot2::theme(
    legend.position = c(0.98, 0.02),
    legend.justification = c(1, 0),
    legend.background = ggplot2::element_rect(fill = grDevices::adjustcolor("white", 0.86), color = "#D0D0D0"),
    legend.key.size = grid::unit(3.6, "mm"),
    legend.text = ggplot2::element_text(size = 6.4)
  )

ggplot2::ggsave(
  file.path(out_dir, "zaslaver_rank1_growth_rate_embedding.png"),
  p_rank_embedding,
  width = 5.2,
  height = 4.8,
  dpi = 300
)
ggplot2::ggsave(
  file.path(out_dir, "zaslaver_rank1_growth_rate_embedding.pdf"),
  p_rank_embedding,
  width = 5.2,
  height = 4.8
)

highlight_terms <- c(
  "ribosomes" = "ribosomes",
  "translation" = "translation",
  "biosynthesis" = "biosynthesis_of_building_blocks",
  "amino acids" = "amino_acids",
  "metabolism" = "metabolism",
  "transport" = "transport",
  "membrane" = "membrane",
  "SOS response" = "SOS_response"
)
highlight_terms <- highlight_terms[highlight_terms %in% names(rank1_embedding)]
highlight_labels <- names(highlight_terms)
embedding_limits <- range(
  c(rank1_embedding$rank1_score_g0_25, rank1_embedding$rank1_score_g0_8),
  finite = TRUE
)
embedding_limits <- embedding_limits + c(-1, 1) * diff(embedding_limits) * 0.04

make_highlight_embedding <- function(label) {
  d <- rank1_embedding
  d$highlight <- d[[highlight_terms[[label]]]] == 1L
  d$highlight[is.na(d$highlight)] <- FALSE
  d <- d[order(d$highlight), , drop = FALSE]
  ggplot2::ggplot(
    d,
    ggplot2::aes(x = rank1_score_g0_25, y = rank1_score_g0_8)
  ) +
    ggplot2::geom_hline(yintercept = 0, color = "#808080", linewidth = 0.25) +
    ggplot2::geom_vline(xintercept = 0, color = "#808080", linewidth = 0.25) +
    ggplot2::geom_point(
      data = d[!d$highlight, , drop = FALSE],
      color = "#D3D3D3",
      size = 0.85,
      alpha = 0.42
    ) +
    ggplot2::geom_point(
      data = d[d$highlight, , drop = FALSE],
      color = unname(embedding_colors[label]),
      size = 1.25,
      alpha = 0.92
    ) +
    ggplot2::coord_cartesian(xlim = embedding_limits, ylim = embedding_limits) +
    ggplot2::labs(
      x = expression("Rank-1 reporter score at " * g == 0.25),
      y = expression("Rank-1 reporter score at " * g == 0.8)
    ) +
    plot_theme(8) +
    ggplot2::theme(
      legend.position = "none",
      plot.title = ggplot2::element_text(face = "bold", size = 9, hjust = 0),
      aspect.ratio = 1
    ) +
    ggplot2::ggtitle(label)
}

highlight_plots <- stats::setNames(lapply(highlight_labels, make_highlight_embedding), highlight_labels)
for (label in highlight_labels) {
  file_label <- gsub("[^A-Za-z0-9]+", "_", tolower(label))
  ggplot2::ggsave(
    file.path(out_dir, paste0("zaslaver_rank1_embedding_highlight_", file_label, ".png")),
    highlight_plots[[label]],
    width = 4.4,
    height = 4.2,
    dpi = 300
  )
  ggplot2::ggsave(
    file.path(out_dir, paste0("zaslaver_rank1_embedding_highlight_", file_label, ".pdf")),
    highlight_plots[[label]],
    width = 4.4,
    height = 4.2
  )
}

highlight_long <- do.call(rbind, lapply(highlight_labels, function(label) {
  d <- rank1_embedding[, c("reporter_id", "reporter", "rank1_score_g0_25", "rank1_score_g0_8", "physiology")]
  d$highlight_label <- label
  d$highlight <- rank1_embedding[[highlight_terms[[label]]]] == 1L
  d$highlight[is.na(d$highlight)] <- FALSE
  d
}))
highlight_long$highlight_label <- factor(highlight_long$highlight_label, levels = highlight_labels)
highlight_long <- highlight_long[order(highlight_long$highlight_label, highlight_long$highlight), , drop = FALSE]

p_highlight_grid <- ggplot2::ggplot(
  highlight_long,
  ggplot2::aes(x = rank1_score_g0_25, y = rank1_score_g0_8)
) +
  ggplot2::geom_hline(yintercept = 0, color = "#808080", linewidth = 0.22) +
  ggplot2::geom_vline(xintercept = 0, color = "#808080", linewidth = 0.22) +
  ggplot2::geom_point(
    data = highlight_long[!highlight_long$highlight, , drop = FALSE],
    color = "#D3D3D3",
    size = 0.45,
    alpha = 0.35
  ) +
  ggplot2::geom_point(
    data = highlight_long[highlight_long$highlight, , drop = FALSE],
    ggplot2::aes(color = highlight_label),
    size = 0.85,
    alpha = 0.9
  ) +
  ggplot2::scale_color_manual(values = embedding_colors[highlight_labels], guide = "none") +
  ggplot2::coord_cartesian(xlim = embedding_limits, ylim = embedding_limits) +
  ggplot2::facet_wrap(~highlight_label, ncol = 4) +
  ggplot2::labs(
    x = expression("Rank-1 reporter score at " * g == 0.25),
    y = expression("Rank-1 reporter score at " * g == 0.8)
  ) +
  plot_theme(7) +
  ggplot2::theme(
    strip.text = ggplot2::element_text(face = "bold", color = "#333333", size = 7.4),
    aspect.ratio = 1,
    panel.grid.minor = ggplot2::element_blank()
  )

ggplot2::ggsave(
  file.path(out_dir, "zaslaver_rank1_embedding_highlight_grid.png"),
  p_highlight_grid,
  width = 8.4,
  height = 5.9,
  dpi = 300
)
ggplot2::ggsave(
  file.path(out_dir, "zaslaver_rank1_embedding_highlight_grid.pdf"),
  p_highlight_grid,
  width = 8.4,
  height = 5.9
)

offdiag_fit <- stats::lm(rank1_score_g0_8 ~ rank1_score_g0_25, data = rank1_embedding)
rank1_embedding$expected_rank1_score_g0_8 <- stats::predict(offdiag_fit, newdata = rank1_embedding)
rank1_embedding$growth_shift_residual <- stats::residuals(offdiag_fit)
rank1_embedding$abs_growth_shift_residual <- abs(rank1_embedding$growth_shift_residual)
rank1_embedding <- rank1_embedding[order(-rank1_embedding$abs_growth_shift_residual), ]
utils::write.table(
  rank1_embedding,
  file.path(out_dir, "zaslaver_rank1_growth_rate_embedding.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
top_shift <- head(rank1_embedding[is.finite(rank1_embedding$growth_shift_residual), ], 18)
top_shift$label <- top_shift$reporter
shift_limit <- stats::quantile(abs(rank1_embedding$growth_shift_residual), probs = 0.985, na.rm = TRUE)

p_offdiag <- ggplot2::ggplot(
  rank1_embedding,
  ggplot2::aes(x = rank1_score_g0_25, y = rank1_score_g0_8)
) +
  ggplot2::geom_hline(yintercept = 0, color = "#808080", linewidth = 0.24) +
  ggplot2::geom_vline(xintercept = 0, color = "#808080", linewidth = 0.24) +
  ggplot2::geom_abline(
    intercept = stats::coef(offdiag_fit)[1],
    slope = stats::coef(offdiag_fit)[2],
    color = "#303030",
    linewidth = 0.45
  ) +
  ggplot2::geom_point(
    ggplot2::aes(color = growth_shift_residual),
    size = 0.95,
    alpha = 0.78
  ) +
  ggplot2::scale_color_gradient2(
    low = "#2166AC",
    mid = "#D9D9D9",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-shift_limit, shift_limit),
    oob = function(x, range) pmin(pmax(x, range[1]), range[2]),
    name = expression(e[a])
  ) +
  ggplot2::coord_cartesian(xlim = embedding_limits, ylim = embedding_limits) +
  ggplot2::labs(
    x = expression("Rank-1 reporter score at " * g == 0.25),
    y = expression("Rank-1 reporter score at " * g == 0.8)
  ) +
  plot_theme(8) +
  ggplot2::theme(
    legend.position = c(0.98, 0.03),
    legend.justification = c(1, 0),
    legend.background = ggplot2::element_rect(fill = grDevices::adjustcolor("white", 0.86), color = "#D0D0D0"),
    aspect.ratio = 1
  )

if (requireNamespace("ggrepel", quietly = TRUE)) {
  p_offdiag <- p_offdiag +
    ggrepel::geom_text_repel(
      data = top_shift,
      ggplot2::aes(label = label),
      size = 2.1,
      min.segment.length = 0,
      box.padding = 0.18,
      max.overlaps = Inf
    )
}

ggplot2::ggsave(
  file.path(out_dir, "zaslaver_rank1_growth_rate_offdiagonal_shift.png"),
  p_offdiag,
  width = 5.4,
  height = 5.0,
  dpi = 300
)
ggplot2::ggsave(
  file.path(out_dir, "zaslaver_rank1_growth_rate_offdiagonal_shift.pdf"),
  p_offdiag,
  width = 5.4,
  height = 5.0
)

residual_quantiles <- stats::quantile(
  rank1_embedding$growth_shift_residual,
  probs = c(0.05, 0.95),
  na.rm = TRUE
)
shift_density <- stats::density(
  rank1_embedding$growth_shift_residual[is.finite(rank1_embedding$growth_shift_residual)],
  adjust = 1.0
)
shift_density_df <- data.frame(
  growth_shift_residual = shift_density$x,
  density = shift_density$y
)
shift_rug <- rbind(
  data.frame(
    rank1_embedding[rank1_embedding$growth_shift_residual <= residual_quantiles[1], c("reporter", "growth_shift_residual")],
    side = "negative"
  ),
  data.frame(
    rank1_embedding[rank1_embedding$growth_shift_residual >= residual_quantiles[2], c("reporter", "growth_shift_residual")],
    side = "positive"
  )
)
shift_label_density <- rbind(
  head(rank1_embedding[order(rank1_embedding$growth_shift_residual), ], 6),
  head(rank1_embedding[order(-rank1_embedding$growth_shift_residual), ], 6)
)
shift_label_density$density <- stats::approx(
  shift_density_df$growth_shift_residual,
  shift_density_df$density,
  xout = shift_label_density$growth_shift_residual,
  rule = 2
)$y
shift_label_density$label <- shift_label_density$reporter

p_shift_distribution <- ggplot2::ggplot(
  shift_density_df,
  ggplot2::aes(x = growth_shift_residual, y = density)
) +
  ggplot2::geom_area(fill = "#D9D9D9", color = "#555555", linewidth = 0.35, alpha = 0.9) +
  ggplot2::geom_vline(xintercept = 0, color = "#303030", linewidth = 0.35) +
  ggplot2::geom_vline(xintercept = residual_quantiles, color = "#B2182B", linetype = "dashed", linewidth = 0.3) +
  ggplot2::geom_rug(
    data = shift_rug,
    ggplot2::aes(x = growth_shift_residual, color = side),
    sides = "b",
    alpha = 0.65,
    linewidth = 0.28,
    inherit.aes = FALSE
  ) +
  ggplot2::scale_color_manual(values = c("negative" = "#2166AC", "positive" = "#B2182B"), guide = "none") +
  ggplot2::labs(
    x = expression("Off-diagonal growth-rate shift " * e[a]),
    y = "Density"
  ) +
  plot_theme(8)

if (requireNamespace("ggrepel", quietly = TRUE)) {
  p_shift_distribution <- p_shift_distribution +
    ggrepel::geom_text_repel(
      data = shift_label_density,
      ggplot2::aes(label = label),
      size = 2.1,
      min.segment.length = 0,
      box.padding = 0.18,
      max.overlaps = Inf
    )
}

ggplot2::ggsave(
  file.path(out_dir, "zaslaver_rank1_growth_rate_shift_residual_distribution.png"),
  p_shift_distribution,
  width = 5.4,
  height = 3.8,
  dpi = 300
)
ggplot2::ggsave(
  file.path(out_dir, "zaslaver_rank1_growth_rate_shift_residual_distribution.pdf"),
  p_shift_distribution,
  width = 5.4,
  height = 3.8
)

shift_enrichment <- do.call(rbind, lapply(c("negative_shift", "positive_shift"), function(side) {
  selected <- if (side == "negative_shift") {
    rank1_embedding$growth_shift_residual <= residual_quantiles[1]
  } else {
    rank1_embedding$growth_shift_residual >= residual_quantiles[2]
  }
  rows <- lapply(highlight_labels, function(label) {
    x <- rank1_embedding[[highlight_terms[[label]]]] == 1L
    x[is.na(x)] <- FALSE
    tab <- table(factor(selected, c(FALSE, TRUE)), factor(x, c(FALSE, TRUE)))
    if (any(dim(tab) != c(2, 2)) || sum(tab[, 2]) < 5) {
      return(NULL)
    }
    ft <- stats::fisher.test(tab, alternative = "greater")
    data.frame(
      side = side,
      term = label,
      selected_with_term = tab[2, 2],
      selected_without_term = tab[2, 1],
      background_with_term = tab[1, 2],
      background_without_term = tab[1, 1],
      odds_ratio = unname(ft$estimate),
      pvalue = ft$p.value,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}))
shift_enrichment$padj <- stats::p.adjust(shift_enrichment$pvalue, method = "BH")
shift_enrichment <- shift_enrichment[order(shift_enrichment$side, shift_enrichment$padj), ]
utils::write.table(
  shift_enrichment,
  file.path(out_dir, "zaslaver_rank1_growth_rate_shift_enrichment.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

density_grid <- seq(embedding_limits[1], embedding_limits[2], length.out = 256)
ridge_rows <- do.call(rbind, lapply(highlight_labels, function(label) {
  selected <- rank1_embedding[[highlight_terms[[label]]]] == 1L
  selected[is.na(selected)] <- FALSE
  do.call(rbind, lapply(c("g = 0.25", "g = 0.8"), function(g_label) {
    score_col <- if (g_label == "g = 0.25") "rank1_score_g0_25" else "rank1_score_g0_8"
    vals <- rank1_embedding[[score_col]][selected]
    vals <- vals[is.finite(vals)]
    if (length(vals) < 2) {
      dens_y <- rep(0, length(density_grid))
    } else {
      dens <- stats::density(vals, from = min(density_grid), to = max(density_grid), n = length(density_grid), adjust = 1.0)
      dens_y <- dens$y / max(dens$y, na.rm = TRUE)
    }
    data.frame(
      highlight_label = label,
      growth_label = g_label,
      x = density_grid,
      density_scaled = dens_y,
      n = length(vals),
      stringsAsFactors = FALSE
    )
  }))
}))
ridge_rows$highlight_label <- factor(ridge_rows$highlight_label, levels = highlight_labels)
ridge_rows$growth_label <- factor(ridge_rows$growth_label, levels = c("g = 0.25", "g = 0.8"))
ridge_rows$y_base <- ifelse(ridge_rows$growth_label == "g = 0.25", 1, 0)
ridge_rows$y <- ridge_rows$y_base + 0.72 * ridge_rows$density_scaled

p_ridges <- ggplot2::ggplot(ridge_rows, ggplot2::aes(x = x)) +
  ggplot2::geom_ribbon(
    ggplot2::aes(ymin = y_base, ymax = y, fill = growth_label),
    alpha = 0.70,
    color = "white",
    linewidth = 0.15
  ) +
  ggplot2::geom_hline(yintercept = c(0, 1), color = "#808080", linewidth = 0.22) +
  ggplot2::scale_fill_manual(values = c("g = 0.25" = "#76B7B2", "g = 0.8" = "#E15759"), name = NULL) +
  ggplot2::facet_wrap(~highlight_label, ncol = 4) +
  ggplot2::scale_y_continuous(
    breaks = c(0, 1),
    labels = c("g = 0.8", "g = 0.25"),
    limits = c(-0.05, 1.9)
  ) +
  ggplot2::labs(
    x = "Rank-1 reporter score",
    y = NULL
  ) +
  plot_theme(7) +
  ggplot2::theme(
    legend.position = c(0.98, 0.02),
    legend.justification = c(1, 0),
    legend.background = ggplot2::element_rect(fill = grDevices::adjustcolor("white", 0.86), color = "#D0D0D0"),
    strip.text = ggplot2::element_text(face = "bold", color = "#333333", size = 7.4),
    panel.grid.minor = ggplot2::element_blank()
  )

ggplot2::ggsave(
  file.path(out_dir, "zaslaver_rank1_physiology_paired_ridgelines.png"),
  p_ridges,
  width = 8.4,
  height = 5.5,
  dpi = 300
)
ggplot2::ggsave(
  file.path(out_dir, "zaslaver_rank1_physiology_paired_ridgelines.pdf"),
  p_ridges,
  width = 8.4,
  height = 5.5
)

rank1_long <- rbind(
  data.frame(
    rank1_table[, c("growth_rate", "reporter_id", "condition_label", "row_order")],
    matrix = "Total effect",
    effect = rank1_table$total_effect
  ),
  data.frame(
    rank1_table[, c("growth_rate", "reporter_id", "condition_label", "row_order")],
    matrix = "Rank-1 component",
    effect = rank1_table$rank1_effect
  )
)
rank1_long$matrix <- factor(rank1_long$matrix, levels = c("Total effect", "Rank-1 component"))
rank1_long$condition_label <- factor(
  rank1_long$condition_label,
  levels = c("Ethanol", "Nitrogen limited", "no AA", "no Glucose", "Phosphate limited")
)
effect_limit <- stats::quantile(abs(rank1_long$effect), probs = 0.985, na.rm = TRUE)

p_total_rank1 <- ggplot2::ggplot(
  rank1_long,
  ggplot2::aes(x = condition_label, y = row_order, fill = effect)
) +
  ggplot2::geom_raster() +
  ggplot2::facet_grid(growth_rate ~ matrix, scales = "free_y", labeller = ggplot2::label_both) +
  ggplot2::scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-effect_limit, effect_limit),
    oob = function(x, range) pmin(pmax(x, range[1]), range[2]),
    name = expression(hat(Delta * y)[aj]^tot)
  ) +
  ggplot2::labs(x = NULL, y = "Promoters ordered by rank-1 loading") +
  plot_theme(8) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
    axis.text.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_blank(),
    legend.key.height = grid::unit(14, "mm")
  )

ggplot2::ggsave(
  file.path(out_dir, "zaslaver_total_effect_and_rank1_heatmap.png"),
  p_total_rank1,
  width = 7.4,
  height = 6.4,
  dpi = 300
)
ggplot2::ggsave(
  file.path(out_dir, "zaslaver_total_effect_and_rank1_heatmap.pdf"),
  p_total_rank1,
  width = 7.4,
  height = 6.4
)

rank1_loading <- unique(rank1_table[, c("growth_rate", "reporter_id", "rank1_row_score")])
annotation_profile_terms <- c(
  "ribosomes",
  "ribosomal_proteins",
  "metabolism",
  "carbon_utilization",
  "nitrogen_metabolism",
  "phosphorous_metabolism",
  "other_.mechanical._nutritional._oxidative_stress.",
  "SOS_response",
  "drug_resistance.sensitivity",
  "translation",
  "Transcription_related"
)
annotation_profile_terms <- intersect(annotation_profile_terms, names(zaslaver_promoter_annotations))
ann_profile <- merge(
  rank1_loading,
  zaslaver_promoter_annotations[, c("reporter_id", annotation_profile_terms), drop = FALSE],
  by = "reporter_id",
  all.x = TRUE,
  sort = FALSE
)
ann_profile$loading_bin <- ave(
  ann_profile$rank1_row_score,
  ann_profile$growth_rate,
  FUN = function(x) {
    cut(
      rank(x, ties.method = "first"),
      breaks = seq(0, length(x), length.out = 21),
      labels = FALSE,
      include.lowest = TRUE
    )
  }
)

annotation_fraction <- do.call(rbind, lapply(annotation_profile_terms, function(term) {
  stats::aggregate(
    ann_profile[[term]] == 1L,
    by = list(growth_rate = ann_profile$growth_rate, loading_bin = ann_profile$loading_bin),
    FUN = mean,
    na.rm = TRUE
  ) |>
    transform(term = term)
}))
names(annotation_fraction)[3] <- "fraction"
annotation_fraction$term_label <- gsub("_", " ", annotation_fraction$term)
annotation_fraction$term_label <- gsub("\\.+", " ", annotation_fraction$term_label)
annotation_fraction$term_label <- trimws(annotation_fraction$term_label)

utils::write.table(
  annotation_fraction,
  file.path(out_dir, "zaslaver_rank1_loading_annotation_profiles.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

extreme_enrichment <- do.call(rbind, lapply(split(ann_profile, ann_profile$growth_rate), function(dg) {
  q <- stats::quantile(dg$rank1_row_score, probs = c(0.05, 0.95), na.rm = TRUE)
  rows <- lapply(c("negative_loading", "positive_loading"), function(side) {
    selected <- if (side == "negative_loading") dg$rank1_row_score <= q[1] else dg$rank1_row_score >= q[2]
    do.call(rbind, lapply(annotation_profile_terms, function(term) {
      x <- dg[[term]] == 1L
      tab <- table(factor(selected, c(FALSE, TRUE)), factor(x, c(FALSE, TRUE)))
      if (any(dim(tab) != c(2, 2)) || sum(tab[, 2]) < 5) {
        return(NULL)
      }
      ft <- stats::fisher.test(tab, alternative = "greater")
      data.frame(
        growth_rate = unique(dg$growth_rate),
        side = side,
        term = term,
        selected_with_term = tab[2, 2],
        selected_without_term = tab[2, 1],
        background_with_term = tab[1, 2],
        background_without_term = tab[1, 1],
        odds_ratio = unname(ft$estimate),
        pvalue = ft$p.value,
        stringsAsFactors = FALSE
      )
    }))
  })
  do.call(rbind, rows)
}))
extreme_enrichment$padj <- stats::p.adjust(extreme_enrichment$pvalue, method = "BH")
extreme_enrichment <- extreme_enrichment[
  order(extreme_enrichment$growth_rate, extreme_enrichment$side, extreme_enrichment$padj),
]
utils::write.table(
  extreme_enrichment,
  file.path(out_dir, "zaslaver_rank1_extreme_annotation_enrichment.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

p_ann_profile <- ggplot2::ggplot(
  annotation_fraction,
  ggplot2::aes(x = loading_bin, y = fraction, color = term_label)
) +
  ggplot2::geom_line(linewidth = 0.45, alpha = 0.9) +
  ggplot2::facet_wrap(~growth_rate, nrow = 2, labeller = ggplot2::label_both) +
  ggplot2::labs(
    x = "Rank-1 loading bin",
    y = "Fraction of promoters",
    color = NULL
  ) +
  plot_theme(8) +
  ggplot2::theme(
    legend.position = "right",
    legend.text = ggplot2::element_text(size = 6.2)
  )

ggplot2::ggsave(
  file.path(out_dir, "zaslaver_rank1_loading_annotation_profiles.png"),
  p_ann_profile,
  width = 7.2,
  height = 5.2,
  dpi = 300
)
ggplot2::ggsave(
  file.path(out_dir, "zaslaver_rank1_loading_annotation_profiles.pdf"),
  p_ann_profile,
  width = 7.2,
  height = 5.2
)

p_dist <- ggplot2::ggplot(
  effect_table,
  ggplot2::aes(x = specific_effect, color = condition_label)
) +
  ggplot2::geom_density(linewidth = 0.45, adjust = 1.1) +
  ggplot2::facet_wrap(~growth_rate, nrow = 2, labeller = ggplot2::label_both) +
  ggplot2::labs(x = expression("Specific condition effect " * hat(delta)[aj]^spec), y = "Density", color = NULL) +
  plot_theme(9) +
  ggplot2::theme(legend.position = c(0.98, 0.98), legend.justification = c(1, 1))

ggplot2::ggsave(
  file.path(out_dir, "zaslaver_specific_effect_density.png"),
  p_dist,
  width = 6.4,
  height = 4.2,
  dpi = 300
)
ggplot2::ggsave(
  file.path(out_dir, "zaslaver_specific_effect_density.pdf"),
  p_dist,
  width = 6.4,
  height = 4.2
)

cat("Wrote Zaslaver exploratory outputs to:", out_dir, "\n")
cat("Effect rows:", nrow(effect_table), "\n")
cat("Top specific-effect table:", file.path(out_dir, "zaslaver_top_specific_effects.tsv"), "\n")

# Exploratory time-course model ------------------------------------------------
#
# This section uses the S1 condition time courses as repeated measurements. It
# is useful for triage and visualization, but the nominal p-values should be
# interpreted with care because adjacent time points are autocorrelated.

if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(export_all = TRUE, quiet = TRUE)
}

tc <- zaslaver_promoter_timecourse
tc$response <- log2(tc$promoter_activity + eps)
tc$log_od <- log2(tc$od + 1e-8)
tc$time_min_numeric <- tc$time_min
tc <- tc[is.finite(tc$response) & is.finite(tc$log_od), , drop = FALSE]

assay_tc <- prepare_assay(
  tc,
  reporter = "reporter_id",
  perturbation = "condition_label",
  control = reference_condition,
  response = "response",
  numeric_covariates = c("time_min_numeric", "log_od")
)

label_lookup <- tc[order(tc$condition_label != reference_condition, tc$reporter_id), c("reporter_id", "reporter")]
label_lookup <- label_lookup[!duplicated(label_lookup$reporter_id), , drop = FALSE]
names(label_lookup)[2] <- "reporter_label"

fit_timecourse_rank <- function(background_rank) {
  fit <- fit_destress(
    assay_tc,
    preset = "model",
    technical = c("time_min_numeric", "log_od"),
    empirical_bayes = TRUE,
    adjustment = "global",
    background_rank = background_rank
  )
  res <- results(fit)
  res <- merge(
    res,
    label_lookup,
    by.x = "reporter",
    by.y = "reporter_id",
    all.x = TRUE,
    sort = FALSE
  )
  res$background_rank <- background_rank
  res$model <- if (background_rank == 0) "rank 0" else paste0("rank ", background_rank)
  res[order(res$specific_padj_global, -abs(res$specific_effect)), ]
}

res_tc_rank0 <- fit_timecourse_rank(0)
res_tc_rank1 <- fit_timecourse_rank(1)
res_tc <- rbind(res_tc_rank0, res_tc_rank1)

utils::write.table(
  res_tc_rank0,
  file.path(out_dir, "zaslaver_timecourse_destress_results_rank0.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
utils::write.table(
  res_tc_rank1,
  file.path(out_dir, "zaslaver_timecourse_destress_results_rank1.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
utils::write.table(
  res_tc,
  file.path(out_dir, "zaslaver_timecourse_destress_results.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

tc_rank_table <- data.frame(
  growth_rate = 1,
  reporter_id = res_tc_rank0$reporter,
  condition_label = res_tc_rank0$perturbation,
  total_effect = res_tc_rank0$total_effect,
  stringsAsFactors = FALSE
)
tc_rank <- background_rank_diagnostics(
  tc_rank_table,
  effect = "total_effect",
  reporter = "reporter_id",
  perturbation = "condition_label",
  rank_max = 5,
  permutations = 200,
  seed = 811
)
utils::write.table(
  tc_rank,
  file.path(out_dir, "zaslaver_timecourse_total_effect_rank_diagnostics.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

rank_compare <- merge(
  res_tc_rank0[, c("reporter", "perturbation", "reporter_label", "total_effect", "specific_effect",
                   "specific_pvalue", "specific_padj_global")],
  res_tc_rank1[, c("reporter", "perturbation", "low_rank_effect", "rank_adjusted_total_effect",
                   "specific_effect", "specific_pvalue", "specific_padj_global")],
  by = c("reporter", "perturbation"),
  suffixes = c("_rank0", "_rank1"),
  all = FALSE,
  sort = FALSE
)
rank_compare$abs_specific_change <- abs(rank_compare$specific_effect_rank1 - rank_compare$specific_effect_rank0)
rank_compare$log10_pvalue_change <- log10(rank_compare$specific_pvalue_rank1) -
  log10(rank_compare$specific_pvalue_rank0)
rank_compare <- rank_compare[order(-abs(rank_compare$low_rank_effect), -rank_compare$abs_specific_change), ]
utils::write.table(
  rank_compare,
  file.path(out_dir, "zaslaver_timecourse_rank1_comparison.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

low_rank_matrix <- res_tc_rank1[, c(
  "reporter", "reporter_label", "perturbation", "total_effect",
  "low_rank_effect", "rank_adjusted_total_effect", "specific_effect",
  "specific_pvalue", "specific_padj_global"
)]
low_rank_matrix <- low_rank_matrix[order(low_rank_matrix$reporter, low_rank_matrix$perturbation), ]
utils::write.table(
  low_rank_matrix,
  file.path(out_dir, "zaslaver_timecourse_low_rank_effect_matrix.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

timecourse_low_rank_score <- stats::aggregate(
  low_rank_effect ~ reporter,
  low_rank_matrix,
  mean,
  na.rm = TRUE
)
names(timecourse_low_rank_score) <- c("reporter_id", "low_rank_score")
timecourse_low_rank_enrichment_data <- merge(
  timecourse_low_rank_score,
  zaslaver_promoter_annotations,
  by = "reporter_id",
  all.x = TRUE,
  sort = FALSE
)
timecourse_low_rank_quantiles <- stats::quantile(
  timecourse_low_rank_enrichment_data$low_rank_score,
  probs = c(0.05, 0.95),
  na.rm = TRUE
)
timecourse_low_rank_enrichment <- do.call(rbind, lapply(c("negative_loading", "positive_loading"), function(side) {
  selected <- if (side == "negative_loading") {
    timecourse_low_rank_enrichment_data$low_rank_score <= timecourse_low_rank_quantiles[1]
  } else {
    timecourse_low_rank_enrichment_data$low_rank_score >= timecourse_low_rank_quantiles[2]
  }
  rows <- lapply(annotation_cols, function(term) {
    x <- timecourse_low_rank_enrichment_data[[term]]
    if (!all(x %in% c(0L, 1L, NA_integer_), na.rm = TRUE)) {
      return(NULL)
    }
    tab <- table(factor(selected, c(FALSE, TRUE)), factor(x == 1L, c(FALSE, TRUE)))
    if (any(dim(tab) != c(2, 2)) || sum(tab[, 2]) < 5) {
      return(NULL)
    }
    ft <- stats::fisher.test(tab, alternative = "greater")
    data.frame(
      side = side,
      term = term,
      selected_with_term = tab[2, 2],
      selected_without_term = tab[2, 1],
      background_with_term = tab[1, 2],
      background_without_term = tab[1, 1],
      odds_ratio = unname(ft$estimate),
      pvalue = ft$p.value,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}))
timecourse_low_rank_enrichment$padj <- stats::p.adjust(timecourse_low_rank_enrichment$pvalue, method = "BH")
timecourse_low_rank_enrichment <- timecourse_low_rank_enrichment[
  order(timecourse_low_rank_enrichment$side, timecourse_low_rank_enrichment$padj),
]
utils::write.table(
  timecourse_low_rank_enrichment,
  file.path(out_dir, "zaslaver_timecourse_low_rank_annotation_enrichment.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

pvalue_summary <- do.call(rbind, lapply(split(res_tc, res_tc$model), function(d) {
  data.frame(
    model = unique(d$model),
    n_tests = nrow(d),
    median_pvalue = stats::median(d$specific_pvalue, na.rm = TRUE),
    prop_pvalue_lt_0.01 = mean(d$specific_pvalue < 0.01, na.rm = TRUE),
    n_fdr_0.05 = sum(d$specific_padj_global < 0.05, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))
utils::write.table(
  pvalue_summary,
  file.path(out_dir, "zaslaver_timecourse_rank_pvalue_summary.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

top_tc <- head(res_tc_rank1[is.finite(res_tc_rank1$specific_pvalue), ], 100)
utils::write.table(
  top_tc,
  file.path(out_dir, "zaslaver_timecourse_top_specific_hits_rank1.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

p_tc_hist <- ggplot2::ggplot(res_tc, ggplot2::aes(x = specific_pvalue)) +
  ggplot2::geom_histogram(bins = 50, fill = "#6BAED6", color = "white", linewidth = 0.15) +
  ggplot2::facet_grid(model ~ perturbation) +
  ggplot2::labs(x = "Raw p-value", y = "Reporter-condition tests") +
  plot_theme(8) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

ggplot2::ggsave(
  file.path(out_dir, "zaslaver_timecourse_specific_pvalue_histograms.png"),
  p_tc_hist,
  width = 8.0,
  height = 4.6,
  dpi = 300
)
ggplot2::ggsave(
  file.path(out_dir, "zaslaver_timecourse_specific_pvalue_histograms.pdf"),
  p_tc_hist,
  width = 8.0,
  height = 4.6
)

volcano_data <- res_tc[is.finite(res_tc$specific_pvalue) & res_tc$specific_pvalue > 0, ]
volcano_data$neg_log10_p <- -log10(volcano_data$specific_pvalue)
volcano_data$hit <- volcano_data$specific_padj_global < 0.05
label_volcano <- head(volcano_data[volcano_data$model == "rank 1", ][
  order(volcano_data$specific_pvalue[volcano_data$model == "rank 1"],
        -abs(volcano_data$specific_effect[volcano_data$model == "rank 1"])),
], 18)
label_name <- if ("reporter_label" %in% names(label_volcano)) label_volcano$reporter_label else label_volcano$reporter
label_volcano$label <- paste(label_name, label_volcano$perturbation, sep = " / ")

p_tc_volcano <- ggplot2::ggplot(
  volcano_data,
  ggplot2::aes(specific_effect, neg_log10_p)
) +
  ggplot2::geom_point(
    ggplot2::aes(color = hit),
    size = 0.75,
    alpha = 0.65
  ) +
  ggplot2::scale_color_manual(values = c("FALSE" = "#BDBDBD", "TRUE" = "#D55E00"), guide = "none") +
  ggplot2::facet_grid(model ~ perturbation) +
  ggplot2::labs(
    x = expression("Specific condition effect " * hat(delta)[aj]^spec),
    y = expression(-log[10] * " raw p-value")
  ) +
  plot_theme(8)

if (requireNamespace("ggrepel", quietly = TRUE)) {
  p_tc_volcano <- p_tc_volcano +
    ggrepel::geom_text_repel(
      data = label_volcano,
      ggplot2::aes(label = label),
      size = 2.1,
      min.segment.length = 0,
      box.padding = 0.25,
      max.overlaps = Inf
    )
}

ggplot2::ggsave(
  file.path(out_dir, "zaslaver_timecourse_specific_volcano.png"),
  p_tc_volcano,
  width = 9.0,
  height = 5.8,
  dpi = 300
)
ggplot2::ggsave(
  file.path(out_dir, "zaslaver_timecourse_specific_volcano.pdf"),
  p_tc_volcano,
  width = 9.0,
  height = 5.8
)

low_rank_plot_data <- res_tc_rank1
low_rank_plot_data$reporter_order <- ave(
  low_rank_plot_data$low_rank_effect,
  low_rank_plot_data$reporter,
  FUN = mean
)
low_rank_plot_data <- low_rank_plot_data[order(low_rank_plot_data$reporter_order), ]
low_rank_plot_data$reporter_index <- match(low_rank_plot_data$reporter, unique(low_rank_plot_data$reporter))
low_rank_limit <- stats::quantile(abs(low_rank_plot_data$low_rank_effect), probs = 0.995, na.rm = TRUE)

p_low_rank <- ggplot2::ggplot(
  low_rank_plot_data,
  ggplot2::aes(x = perturbation, y = reporter_index, fill = low_rank_effect)
) +
  ggplot2::geom_tile() +
  ggplot2::scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-low_rank_limit, low_rank_limit),
    oob = function(x, range) pmin(pmax(x, range[1]), range[2]),
    name = expression(hat(delta)[aj]^rank)
  ) +
  ggplot2::labs(x = NULL, y = "Reporters ordered by fitted low-rank effect") +
  plot_theme(8) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
    axis.text.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    legend.key.height = grid::unit(13, "mm")
  )

ggplot2::ggsave(
  file.path(out_dir, "zaslaver_timecourse_low_rank_effect_heatmap.png"),
  p_low_rank,
  width = 5.6,
  height = 5.4,
  dpi = 300
)
ggplot2::ggsave(
  file.path(out_dir, "zaslaver_timecourse_low_rank_effect_heatmap.pdf"),
  p_low_rank,
  width = 5.6,
  height = 5.4
)

cat("Exploratory time-course DStressR rows:", nrow(res_tc), "\n")
cat("Rank-0 FDR < 0.05:", pvalue_summary$n_fdr_0.05[pvalue_summary$model == "rank 0"], "\n")
cat("Rank-1 FDR < 0.05:", pvalue_summary$n_fdr_0.05[pvalue_summary$model == "rank 1"], "\n")
