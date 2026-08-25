source(file.path("analysis", "_helpers.R"))
load_destress_package()

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package `ggplot2` is required for interaction-sensitivity plots.", call. = FALSE)
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
comparison$modeled_hit <- as.logical(comparison$modeled_hit)
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
  dstressr_without_ev = "modeled_hit",
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
      sum(comparison$interaction_without_ev_hit & comparison$modeled_hit, na.rm = TRUE),
      sum(comparison$interaction_with_ev_hit & comparison$evc_huber_hit, na.rm = TRUE),
      sum(comparison$binsfeld_hit & comparison$modeled_hit & comparison$evc_huber_hit &
            comparison$interaction_without_ev_hit & comparison$interaction_with_ev_hit, na.rm = TRUE),
      sum(comparison$binsfeld_hit & comparison$interaction_without_ev_hit, na.rm = TRUE),
      sum(comparison$binsfeld_hit & comparison$interaction_with_ev_hit, na.rm = TRUE),
      sum(comparison$modeled_hit & comparison$interaction_without_ev_hit, na.rm = TRUE),
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

set_labels <- c(
  binsfeld = "Binsfeld reference",
  dstressr_without_ev = "DStressR without EV",
  dstressr_with_ev = "DStressR with EV",
  interaction_without_ev = "Interaction without EV",
  interaction_with_ev = "Interaction with EV"
)

hit_matrix <- comparison[, hit_cols, drop = FALSE]
hit_matrix[] <- lapply(hit_matrix, function(x) as.logical(x) & !is.na(x))
names(hit_matrix) <- names(hit_cols)
hit_matrix$pattern <- apply(hit_matrix[, names(hit_cols), drop = FALSE], 1, function(x) {
  paste(as.integer(x), collapse = "")
})
intersection_counts <- stats::aggregate(
  count ~ pattern,
  data.frame(pattern = hit_matrix$pattern, count = 1L),
  sum
)
intersection_counts <- intersection_counts[intersection_counts$pattern != "00000", , drop = FALSE]
intersection_counts <- intersection_counts[order(intersection_counts$count, decreasing = TRUE), , drop = FALSE]
intersection_counts$intersection <- seq_len(nrow(intersection_counts))
intersection_counts$pattern <- factor(
  intersection_counts$pattern,
  levels = intersection_counts$pattern
)

membership <- do.call(rbind, lapply(seq_len(nrow(intersection_counts)), function(i) {
  active <- strsplit(as.character(intersection_counts$pattern[i]), "", fixed = TRUE)[[1]] == "1"
  data.frame(
    intersection = intersection_counts$intersection[i],
    method = factor(names(hit_cols), levels = rev(names(hit_cols))),
    active = active,
    stringsAsFactors = FALSE
  )
}))
membership$method_label <- factor(
  set_labels[as.character(membership$method)],
  levels = set_labels[rev(names(hit_cols))]
)

set_sizes <- data.frame(
  method = factor(names(hit_cols), levels = rev(names(hit_cols))),
  method_label = factor(set_labels[names(hit_cols)], levels = set_labels[rev(names(hit_cols))]),
  size = vapply(names(hit_cols), function(nm) sum(hit_matrix[[nm]], na.rm = TRUE), numeric(1)),
  stringsAsFactors = FALSE
)

if (requireNamespace("patchwork", quietly = TRUE)) {
  patchwork <- asNamespace("patchwork")
  bar_plot <- ggplot2$ggplot(intersection_counts, ggplot2$aes(intersection, count)) +
    ggplot2$geom_col(fill = "#374151", width = 0.72) +
    ggplot2$geom_text(ggplot2$aes(label = count), vjust = -0.25, size = 3.2) +
    ggplot2$scale_x_continuous(expand = ggplot2$expansion(mult = c(0.01, 0.01))) +
    ggplot2$scale_y_continuous(expand = ggplot2$expansion(mult = c(0, 0.12))) +
    ggplot2$labs(x = NULL, y = "Intersection size") +
    ggplot2$theme_bw(base_size = 11) +
    ggplot2$theme(
      axis.text.x = ggplot2$element_blank(),
      axis.ticks.x = ggplot2$element_blank(),
      panel.grid.major.x = ggplot2$element_blank(),
      panel.grid.minor = ggplot2$element_blank()
    )

  matrix_plot <- ggplot2$ggplot(membership, ggplot2$aes(intersection, method_label)) +
    ggplot2$geom_line(
      data = membership[membership$active, , drop = FALSE],
      ggplot2$aes(group = intersection),
      color = "#374151",
      linewidth = 0.45
    ) +
    ggplot2$geom_point(
      data = membership[!membership$active, , drop = FALSE],
      color = "#d1d5db",
      size = 2.5
    ) +
    ggplot2$geom_point(
      data = membership[membership$active, , drop = FALSE],
      color = "#374151",
      size = 3.1
    ) +
    ggplot2$scale_x_continuous(expand = ggplot2$expansion(mult = c(0.01, 0.01))) +
    ggplot2$labs(x = "Hit-set intersection", y = NULL) +
    ggplot2$theme_bw(base_size = 11) +
    ggplot2$theme(
      panel.grid.major.x = ggplot2$element_blank(),
      panel.grid.minor = ggplot2$element_blank(),
      axis.text.x = ggplot2$element_blank(),
      axis.ticks.x = ggplot2$element_blank()
    )

  set_plot <- ggplot2$ggplot(set_sizes, ggplot2$aes(size, method_label)) +
    ggplot2$geom_col(fill = "#9ca3af", width = 0.6) +
    ggplot2$geom_text(ggplot2$aes(label = size), hjust = -0.15, size = 3.1) +
    ggplot2$scale_x_continuous(expand = ggplot2$expansion(mult = c(0, 0.2))) +
    ggplot2$labs(x = "Set size", y = NULL) +
    ggplot2$theme_bw(base_size = 11) +
    ggplot2$theme(
      panel.grid.major.y = ggplot2$element_blank(),
      panel.grid.minor = ggplot2$element_blank()
    )

  joint_intersection_plot <- (bar_plot / matrix_plot) | set_plot
  joint_intersection_plot <- joint_intersection_plot +
    patchwork$plot_layout(widths = c(4.6, 1.35), heights = c(2.2, 1.45)) +
    patchwork$plot_annotation(
      title = "Joint hit-set intersections for primary and interaction analyses"
    )

  ggplot2$ggsave(
    file.path(out_dir, "interaction_joint_intersection_plot.png"),
    joint_intersection_plot,
    width = 9.2,
    height = 5.4,
    dpi = 300
  )
  ggplot2$ggsave(
    file.path(out_dir, "interaction_joint_intersection_plot.pdf"),
    joint_intersection_plot,
    width = 9.2,
    height = 5.4
  )
}

intersection_ledger <- intersection_counts
venn_short <- c("B", "D-", "D+", "I-", "I+")
names(venn_short) <- names(hit_cols)
decode_pattern <- function(pattern) {
  active <- strsplit(as.character(pattern), "", fixed = TRUE)[[1]] == "1"
  paste(venn_short[names(hit_cols)[active]], collapse = " + ")
}
intersection_ledger$Combination <- vapply(intersection_ledger$pattern, decode_pattern, character(1))
intersection_ledger$rank <- seq_len(nrow(intersection_ledger))
utils::write.table(
  intersection_ledger[, c("Combination", "count", "pattern")],
  file.path(out_dir, "interaction_joint_intersections.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

if (!requireNamespace("VennDiagram", quietly = TRUE)) {
  stop("Package `VennDiagram` is required to draw the five-set Venn diagram.", call. = FALSE)
}

pair_id <- paste(comparison$promoter, comparison$compound, sep = "\r")
venn_order <- c(
  "binsfeld",
  "dstressr_without_ev",
  "interaction_without_ev",
  "interaction_with_ev",
  "dstressr_with_ev"
)
venn_hit_sets <- lapply(venn_order, function(nm) {
  unique(pair_id[comparison[[hit_cols[[nm]]]]])
})
names(venn_hit_sets) <- c(
  "Binsfeld\nref.",
  "DStressR\nno EV",
  "Interaction\nno EV",
  "Interaction\nwith EV",
  "DStressR\nwith EV"
)
venn_fill <- c("#5b8fd9", "#d8b56d", "#f1d99b", "#bd8068", "#9b6a55")
venn_outline <- c("#254f91", "#9f7625", "#b58c3d", "#7c4b38", "#633b2d")
venn_object <- VennDiagram::venn.diagram(
  x = venn_hit_sets,
  filename = NULL,
  category.names = rep("", 5),
  fill = venn_fill,
  col = venn_outline,
  alpha = 0.42,
  lwd = 2,
  lty = "solid",
  cex = 1.1,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 0.9,
  cat.fontface = "bold",
  cat.fontfamily = "sans",
  cat.col = venn_outline,
  cat.default.pos = "outer",
  cat.dist = c(0.24, 0.13, 0.17, 0.17, 0.13),
  margin = 0.17
)
venn_object <- lapply(venn_object, function(grob) {
  if (inherits(grob, "text") && identical(grob$label, "0")) {
    grob$label <- ""
  }
  grob
})
class(venn_object) <- c("VennDiagram", "gList")

draw_joint_venn <- function() {
  grid::grid.newpage()
  grid::grid.draw(venn_object)
  grid::grid.text("Binsfeld\nreference", x = 0.50, y = 0.965, gp = grid::gpar(col = venn_outline[1], fontsize = 15, fontface = "bold", fontfamily = "sans", lineheight = 0.9))
  grid::grid.text("DStressR\nwithout EV", x = 0.075, y = 0.565, gp = grid::gpar(col = venn_outline[2], fontsize = 15, fontface = "bold", fontfamily = "sans", lineheight = 0.9))
  grid::grid.text("Interaction\nwithout EV", x = 0.235, y = 0.065, gp = grid::gpar(col = venn_outline[3], fontsize = 15, fontface = "bold", fontfamily = "sans", lineheight = 0.9))
  grid::grid.text("Interaction\nwith EV", x = 0.805, y = 0.065, gp = grid::gpar(col = venn_outline[4], fontsize = 15, fontface = "bold", fontfamily = "sans", lineheight = 0.9))
  grid::grid.text("DStressR\nwith EV", x = 0.925, y = 0.605, gp = grid::gpar(col = venn_outline[5], fontsize = 15, fontface = "bold", fontfamily = "sans", lineheight = 0.9))
}

grDevices::png(
  file.path(out_dir, "interaction_joint_venn_plot.png"),
  width = 2400,
  height = 1900,
  res = 300,
  bg = "white"
)
draw_joint_venn()
grDevices::dev.off()

grDevices::pdf(
  file.path(out_dir, "interaction_joint_venn_plot.pdf"),
  width = 8,
  height = 6.3,
  bg = "white"
)
draw_joint_venn()
grDevices::dev.off()

key_pairs <- data.frame(
  promoter = c("soxSp", "micFp", "micFp", "ompFp", "marRABp", "marRABp", "robp", "soxSp"),
  compound = c("Paraquat", "Caffeine", "Procaine", "Procaine", "Chloramphenicol", "Vanillin", "Paraquat", "Caffeine"),
  stringsAsFactors = FALSE
)
key_pairs <- merge(key_pairs, comparison, by = c("promoter", "compound"), all.x = TRUE, sort = FALSE)
key_cols <- c(
  "promoter", "compound",
  "binsfeld_hit", "modeled_hit", "evc_huber_hit",
  "interaction_without_ev_hit", "interaction_with_ev_hit",
  "mean_z", "modeled_effect", "evc_huber_effect",
  "interaction_without_ev_effect", "interaction_with_ev_effect",
  "binsfeld_padj", "modeled_padj_by_reporter", "evc_huber_padj_by_reporter",
  "interaction_without_ev_padj_by_reporter", "interaction_with_ev_padj_by_reporter"
)
utils::write.table(
  key_pairs[, key_cols, drop = FALSE],
  file.path(out_dir, "interaction_key_pair_comparison.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

plot_data <- rbind(
  data.frame(
    variant = "without EV",
    primary_effect = comparison$modeled_effect,
    interaction_effect = comparison$interaction_without_ev_effect,
    primary_pvalue = comparison$modeled_pvalue,
    interaction_pvalue = comparison$interaction_without_ev_pvalue,
    primary_hit = comparison$modeled_hit,
    interaction_hit = comparison$interaction_without_ev_hit
  ),
  data.frame(
    variant = "with EV",
    primary_effect = comparison$evc_huber_effect,
    interaction_effect = comparison$interaction_with_ev_effect,
    primary_pvalue = comparison$evc_huber_pvalue,
    interaction_pvalue = comparison$interaction_with_ev_pvalue,
    primary_hit = comparison$evc_huber_hit,
    interaction_hit = comparison$interaction_with_ev_hit
  )
)
plot_data$status <- ifelse(
  plot_data$primary_hit & plot_data$interaction_hit,
  "Both",
  ifelse(plot_data$primary_hit, "Primary only", ifelse(plot_data$interaction_hit, "Interaction only", "Neither"))
)

effect_plot <- ggplot2$ggplot(
  plot_data,
  ggplot2$aes(primary_effect, interaction_effect, color = status)
) +
  ggplot2$geom_hline(yintercept = 0, color = "#d1d5db", linewidth = 0.3) +
  ggplot2$geom_vline(xintercept = 0, color = "#d1d5db", linewidth = 0.3) +
  ggplot2$geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#6b7280") +
  ggplot2$geom_point(alpha = 0.75, size = 1.6) +
  ggplot2$facet_wrap(~ variant, scales = "free") +
  ggplot2$scale_color_manual(values = c(
    "Both" = "#1f77b4",
    "Primary only" = "#d55e00",
    "Interaction only" = "#009e73",
    "Neither" = "#9ca3af"
  )) +
  ggplot2$labs(
    x = "Primary DStressR specific effect",
    y = "Interaction-model specific effect",
    color = "Hit status"
  ) +
  ggplot2$theme_bw(base_size = 11) +
  ggplot2$theme(
    legend.position = "bottom",
    strip.background = ggplot2$element_rect(fill = "#f3f4f6", color = "#d1d5db")
  )

ggplot2$ggsave(
  file.path(out_dir, "interaction_effect_comparison.png"),
  effect_plot,
  width = 7.2,
  height = 4.6,
  dpi = 300
)
ggplot2$ggsave(
  file.path(out_dir, "interaction_effect_comparison.pdf"),
  effect_plot,
  width = 7.2,
  height = 4.6
)

pvalue_plot <- ggplot2$ggplot(
  plot_data,
  ggplot2$aes(safe_neglog10(primary_pvalue), safe_neglog10(interaction_pvalue), color = status)
) +
  ggplot2$geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#6b7280") +
  ggplot2$geom_point(alpha = 0.75, size = 1.6) +
  ggplot2$facet_wrap(~ variant, scales = "free") +
  ggplot2$scale_color_manual(values = c(
    "Both" = "#1f77b4",
    "Primary only" = "#d55e00",
    "Interaction only" = "#009e73",
    "Neither" = "#9ca3af"
  )) +
  ggplot2$labs(
    x = "Primary DStressR -log10 raw p-value",
    y = "Interaction-model -log10 raw p-value",
    color = "Hit status"
  ) +
  ggplot2$theme_bw(base_size = 11) +
  ggplot2$theme(
    legend.position = "bottom",
    strip.background = ggplot2$element_rect(fill = "#f3f4f6", color = "#d1d5db")
  )

ggplot2$ggsave(
  file.path(out_dir, "interaction_pvalue_comparison.png"),
  pvalue_plot,
  width = 7.2,
  height = 4.6,
  dpi = 300
)
ggplot2$ggsave(
  file.path(out_dir, "interaction_pvalue_comparison.pdf"),
  pvalue_plot,
  width = 7.2,
  height = 4.6
)

message("Interaction sensitivity complete: ", out_dir)
