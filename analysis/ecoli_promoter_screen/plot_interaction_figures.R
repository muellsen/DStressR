source(file.path("analysis", "_helpers.R"))
source(file.path("analysis", "_paper_artifacts.R"))
load_destress_package()
set.seed(1) # Reproducible label placement only.
ggplot2 <- asNamespace("ggplot2")
target_arg <- grep("^--target=", commandArgs(trailingOnly = TRUE), value = TRUE)
target <- if (length(target_arg)) sub("^--target=", "", target_arg) else "total"
if (length(target) != 1L || !target %in% c("total", "relative")) {
  stop("Use --target=total or --target=relative.", call. = FALSE)
}
out_dir <- analysis_output_dir(if (target == "total") {
  "ecoli_interaction_total_effects"
} else {
  "ecoli_interaction_model"
})
artifact <- read_paper_artifact(file.path(out_dir, "binsfeld_interaction_paper.rds"))
comparison <- artifact$annotations$comparison
hit_cols <- artifact$annotations$hit_cols
safe_neglog10 <- function(x) -log10(pmax(as.numeric(x), .Machine$double.xmin))
set_labels <- c(
  binsfeld = "Binsfeld reference",
  dstressr_without_ev = "DStressR without EVC",
  dstressr_with_ev = "DStressR with EVC",
  interaction_without_ev = "Interaction without EVC",
  interaction_with_ev = "Interaction with EVC"
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
  "DStressR\nno EVC",
  "Interaction\nno EVC",
  "Interaction\nwith EVC",
  "DStressR\nwith EVC"
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
  grid::grid.text("DStressR\nwithout EVC", x = 0.075, y = 0.565, gp = grid::gpar(col = venn_outline[2], fontsize = 13, fontface = "bold", fontfamily = "sans", lineheight = 0.9))
  grid::grid.text("Interaction\nwithout EVC", x = 0.235, y = 0.065, gp = grid::gpar(col = venn_outline[3], fontsize = 15, fontface = "bold", fontfamily = "sans", lineheight = 0.9))
  grid::grid.text("Interaction\nwith EVC", x = 0.805, y = 0.065, gp = grid::gpar(col = venn_outline[4], fontsize = 15, fontface = "bold", fontfamily = "sans", lineheight = 0.9))
  grid::grid.text("DStressR\nwith EVC", x = 0.925, y = 0.605, gp = grid::gpar(col = venn_outline[5], fontsize = 13, fontface = "bold", fontfamily = "sans", lineheight = 0.9))
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
if (target == "total") {
  key_pairs <- data.frame(
    promoter = c("soxSp", "micFp", "micFp", "ompFp", "marRABp", "robp", "acrABp",
                 "micFp", "micFp", "soxSp"),
    compound = c("Paraquat", "Caffeine", "Procaine", "Procaine", "Clarithromycin",
                 "Clarithromycin", "Clarithromycin", "Paraquat", "Vanillin", "Caffeine")
  )
}
key_pairs <- merge(key_pairs, comparison, by = c("promoter", "compound"), all.x = TRUE, sort = FALSE)
key_cols <- c(
  "promoter", "compound",
  "binsfeld_hit", "estimated_exponent_hit", "evc_huber_hit",
  "interaction_without_ev_hit", "interaction_with_ev_hit",
  "mean_z", "estimated_exponent_effect", "evc_huber_effect",
  "interaction_without_ev_effect", "interaction_with_ev_effect",
  "binsfeld_padj", "estimated_exponent_padj_by_reporter", "evc_huber_padj_by_reporter",
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
    variant = "Without EVC",
    primary_effect = comparison$estimated_exponent_effect,
    interaction_effect = comparison$interaction_without_ev_effect,
    primary_pvalue = comparison$estimated_exponent_pvalue,
    interaction_pvalue = comparison$interaction_without_ev_pvalue,
    primary_hit = comparison$estimated_exponent_hit,
    interaction_hit = comparison$interaction_without_ev_hit
  ),
  data.frame(
    variant = "With EVC",
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
  ifelse(plot_data$primary_hit, "Reporter-wise only", ifelse(plot_data$interaction_hit, "Interaction only", "Neither"))
)
plot_data$variant <- factor(plot_data$variant, levels = c("Without EVC", "With EVC"))

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
    "Reporter-wise only" = "#d55e00",
    "Interaction only" = "#009e73",
    "Neither" = "#9ca3af"
  )) +
  ggplot2$labs(
    x = paste("Reporter-wise", target, "effect"),
    y = paste("Interaction-model", target, "effect"),
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
    "Reporter-wise only" = "#d55e00",
    "Interaction only" = "#009e73",
    "Neither" = "#9ca3af"
  )) +
  ggplot2$labs(
    x = "Reporter-wise -log10 raw p-value",
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

message("Wrote interaction-model figures: ", out_dir)
