#!/usr/bin/env Rscript

# Paper-style figures for the PFIbox promoter-activity subset workflow.
# This script assumes that run_pfibox_promoter_activity_subset.R has already
# created pair-level DStressR results for y_i^{PA}(T).

source(file.path("analysis", "_helpers.R"))
load_destress_package()

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package `ggplot2` is required for PFIbox PA figures.", call. = FALSE)
}
if (!requireNamespace("gridExtra", quietly = TRUE)) {
  stop("Package `gridExtra` is required for PFIbox PA figures.", call. = FALSE)
}
if (!requireNamespace("scales", quietly = TRUE)) {
  stop("Package `scales` is required for PFIbox PA figures.", call. = FALSE)
}

ggplot2 <- asNamespace("ggplot2")
gridExtra <- asNamespace("gridExtra")
scales <- asNamespace("scales")

in_dir <- analysis_output_dir("french_pfibox", "promoter_activity_subset")
run_tag <- Sys.getenv("PFIBOX_RUN_TAG", unset = "")
if (!nzchar(run_tag)) {
  candidates <- list.files(
    in_dir,
    pattern = "^pfibox_pa_destress_pair_results_.*[.]tsv$",
    full.names = TRUE
  )
  if (!length(candidates)) {
    stop(
      "No PFIbox PA result files found in ", in_dir,
      ". Run analysis/french_pfibox/run_pfibox_promoter_activity_subset.R first.",
      call. = FALSE
    )
  }
  candidates <- candidates[order(file.info(candidates)$mtime, decreasing = TRUE)]
  result_file <- candidates[1]
  run_tag <- sub("^pfibox_pa_destress_pair_results_", "", basename(result_file))
  run_tag <- sub("[.]tsv$", "", run_tag)
} else {
  result_file <- file.path(in_dir, paste0("pfibox_pa_destress_pair_results_", run_tag, ".tsv"))
}

response_file <- file.path(in_dir, paste0("pfibox_pa_response_table_", run_tag, ".tsv"))
diagnostic_file <- file.path(in_dir, paste0("pfibox_pa_perturbation_diagnostics_", run_tag, ".tsv"))
metadata_file <- file.path(in_dir, paste0("pfibox_pa_run_metadata_", run_tag, ".tsv"))

for (path in c(result_file, response_file, diagnostic_file, metadata_file)) {
  if (!file.exists(path)) {
    stop("Missing required PFIbox PA workflow output: ", path, call. = FALSE)
  }
}

res <- read_tsv_base(result_file)
response_tab <- read_tsv_base(response_file)
diagnostics <- read_tsv_base(diagnostic_file)
metadata <- read_tsv_base(metadata_file)

required_res <- c(
  "reporter", "perturbation", "total_effect", "specific_effect",
  "specific_pvalue", "specific_padj_by_reporter"
)
missing_res <- setdiff(required_res, names(res))
if (length(missing_res)) {
  stop("Missing result columns: ", paste(missing_res, collapse = ", "), call. = FALSE)
}

t_label <- paste0("[", metadata$realized_t_start_h[1], ", ", metadata$realized_t_end_h[1], "] h")

paper_theme <- function(base_size = 10) {
  ggplot2$theme_light(base_size = base_size) +
    ggplot2$theme(
      panel.grid.minor = ggplot2$element_blank(),
      panel.grid.major = ggplot2$element_line(linewidth = 0.24, color = "#E5E7EB"),
      axis.title = ggplot2$element_text(color = "#111827"),
      axis.text = ggplot2$element_text(color = "#4B5563"),
      legend.background = ggplot2$element_blank(),
      legend.key = ggplot2$element_blank(),
      legend.title = ggplot2$element_text(size = base_size - 1),
      legend.text = ggplot2$element_text(size = base_size - 2),
      plot.title = ggplot2$element_blank(),
      plot.subtitle = ggplot2$element_blank()
    )
}

panel_label <- function(label, plot) {
  gridExtra$arrangeGrob(
    grid::textGrob(
      label,
      x = grid::unit(0, "npc"),
      y = grid::unit(0.25, "npc"),
      hjust = 0,
      vjust = 0,
      gp = grid::gpar(fontsize = 14, fontface = "bold", col = "#111827")
    ),
    plot,
    ncol = 1,
    heights = c(0.06, 1)
  )
}

matrix_order <- function(tab, value) {
  wide <- stats::xtabs(stats::as.formula(paste(value, "~ reporter + perturbation")), data = tab)
  row_score <- apply(abs(wide), 1, max, na.rm = TRUE)
  rows <- names(sort(row_score, decreasing = TRUE))
  col_score <- colSums(abs(wide), na.rm = TRUE)
  cols <- names(sort(col_score, decreasing = TRUE))
  list(rows = rows, cols = cols)
}

make_effect_matrix_plot <- function(tab, value, legend_title, rows = NULL, cols = NULL, show_y = FALSE) {
  d <- tab[, c("reporter", "perturbation", value), drop = FALSE]
  names(d)[3] <- "effect"
  d$effect <- as.numeric(d$effect)
  d <- d[is.finite(d$effect), , drop = FALSE]
  if (is.null(rows) || is.null(cols)) {
    ord <- matrix_order(tab, value)
    rows <- ord$rows
    cols <- ord$cols
  }
  d$reporter <- factor(d$reporter, levels = rev(rows))
  d$perturbation <- factor(d$perturbation, levels = cols)
  limit <- stats::quantile(abs(d$effect), 0.98, na.rm = TRUE, names = FALSE)
  if (!is.finite(limit) || limit <= 0) {
    limit <- max(abs(d$effect), na.rm = TRUE)
  }
  ggplot2$ggplot(d, ggplot2$aes(perturbation, reporter, fill = effect)) +
    ggplot2$geom_tile(width = 0.98, height = 0.98, color = "white", linewidth = 0.02) +
    ggplot2$scale_fill_gradient2(
      low = "#3f6fb5",
      mid = "white",
      high = "#c93f3f",
      midpoint = 0,
      limits = c(-limit, limit),
      oob = scales$squish,
      name = legend_title
    ) +
    ggplot2$labs(x = NULL, y = NULL) +
    ggplot2$theme_minimal(base_size = 9) +
    ggplot2$theme(
      panel.grid = ggplot2$element_blank(),
      axis.text.x = ggplot2$element_text(angle = 50, hjust = 1, vjust = 1, size = 8.5),
      axis.text.y = if (show_y) ggplot2$element_text(size = 3.2) else ggplot2$element_blank(),
      axis.ticks.y = ggplot2$element_blank(),
      legend.position = "right",
      legend.title = ggplot2$element_text(size = 8),
      legend.text = ggplot2$element_text(size = 7),
      plot.margin = ggplot2$margin(0, 6, 2, 4)
    ) +
    ggplot2$guides(
      fill = ggplot2$guide_colorbar(
        title.position = "top",
        title.hjust = 0.5,
        barheight = grid::unit(20, "mm"),
        barwidth = grid::unit(3.2, "mm")
      )
    )
}

ord <- matrix_order(res, "total_effect")
p_total_matrix <- make_effect_matrix_plot(
  res,
  "total_effect",
  expression(hat(Delta*y)[aj]^"tot, PA"),
  rows = ord$rows,
  cols = ord$cols,
  show_y = FALSE
)
p_specific_matrix <- make_effect_matrix_plot(
  res,
  "specific_effect",
  expression(hat(delta)[aj]^PA),
  rows = ord$rows,
  cols = ord$cols,
  show_y = FALSE
)

label_data <- diagnostics[order(-diagnostics$variance_residual, -diagnostics$effect_variance), , drop = FALSE]
label_data <- utils::head(label_data[is.finite(label_data$effect_variance), , drop = FALSE], 8)
color_abs <- abs(diagnostics$mean_effect[is.finite(diagnostics$mean_effect)])
color_limit <- stats::quantile(color_abs, 0.92, names = FALSE, na.rm = TRUE)
if (!is.finite(color_limit) || color_limit <= 0) {
  color_limit <- max(color_abs, na.rm = TRUE)
}
diagnostics$color_effect <- pmax(pmin(diagnostics$mean_effect, color_limit), -color_limit)

p_rank_variance <- ggplot2$ggplot(
  diagnostics,
  ggplot2$aes(rank_abs_mean_effect, effect_variance)
) +
  ggplot2$geom_point(ggplot2$aes(color = color_effect), size = 1.7, alpha = 0.84) +
  ggplot2$geom_line(
    data = diagnostics[is.finite(diagnostics$variance_trend), , drop = FALSE],
    ggplot2$aes(y = variance_trend),
    color = "#008A5B",
    linewidth = 0.7
  ) +
  ggplot2$scale_color_gradient2(
    low = "#1D4ED8",
    mid = "#9CA3AF",
    high = "#B91C1C",
    midpoint = 0,
    limits = c(-color_limit, color_limit),
    name = expression(hat(Delta * bar(y))[j]^"tot, PA")
  ) +
  ggplot2$labs(
    x = expression("Rank of " * abs(hat(Delta * bar(y))[j]^"tot, PA")),
    y = expression("Variance across reporters of " * hat(Delta * y)["\u00b7" * j]^"tot, PA")
  ) +
  ggplot2$scale_y_continuous(expand = ggplot2$expansion(mult = c(0.08, 0.20))) +
  paper_theme(10) +
  ggplot2$theme(
    legend.position = "right",
    legend.justification = "center",
    legend.background = ggplot2$element_blank(),
    legend.key.height = grid::unit(10, "mm"),
    legend.key.width = grid::unit(3, "mm")
  )
if (nrow(label_data) > 0 && requireNamespace("ggrepel", quietly = TRUE)) {
  p_rank_variance <- p_rank_variance +
    ggrepel::geom_text_repel(
      data = label_data,
      ggplot2$aes(label = perturbation_label),
      size = 3.1,
      fontface = "bold",
      min.segment.length = 0,
      box.padding = 0.2,
      point.padding = 0.12,
      seed = 4,
      max.overlaps = Inf,
      ylim = c(NA, Inf),
      show.legend = FALSE
    )
}

effect_matrix_figure <- gridExtra$arrangeGrob(
  panel_label("a", p_total_matrix),
  panel_label("b", p_specific_matrix),
  panel_label("c", p_rank_variance),
  ncol = 1,
  heights = c(1.05, 1.05, 0.85)
)

effect_base <- file.path(in_dir, paste0("pfibox_pa_effect_matrix_and_variance_", run_tag))
ggplot2$ggsave(paste0(effect_base, ".png"), effect_matrix_figure, width = 10.6, height = 11.8, dpi = 300, bg = "white")
ggplot2$ggsave(paste0(effect_base, ".pdf"), effect_matrix_figure, width = 10.6, height = 11.8, bg = "white")

p_hist <- ggplot2$ggplot(res[is.finite(res$specific_pvalue), , drop = FALSE],
                         ggplot2$aes(specific_pvalue)) +
  ggplot2$geom_histogram(bins = 45, fill = "#4E79A7", color = "white", linewidth = 0.18) +
  ggplot2$coord_cartesian(xlim = c(0, 1)) +
  ggplot2$labs(
    x = "Specific-effect p-value",
    y = "Reporter-perturbation pairs"
  ) +
  paper_theme(10)

p_volcano <- plot_volcano(
  res,
  effect = "specific_effect",
  pvalue = "specific_pvalue",
  padj = "specific_padj_by_reporter",
  y = "pvalue",
  reporter = "reporter",
  perturbation = "perturbation",
  fdr = 0.05,
  top_n = 14,
  title = NULL,
  subtitle = NULL,
  xlab = expression(hat(delta)[aj]^PA),
  ylab = expression(-log[10]("p-value")),
  call_labels = c("Not called", "FDR < 0.05"),
  max_label_chars = 36,
  point_size = 1.15,
  point_alpha = 0.74,
  hit_color = "#D95F02",
  background_color = "#CBD5E1"
) +
  ggplot2$theme(
    plot.title = ggplot2$element_blank(),
    plot.subtitle = ggplot2$element_blank(),
    legend.position = "right",
    legend.justification = "center",
    legend.background = ggplot2$element_blank(),
    legend.key = ggplot2$element_blank()
  ) +
  ggplot2$guides(
    color = ggplot2$guide_legend(
      title = NULL,
      override.aes = list(size = 2.2, alpha = 1)
    )
  )

hit_table <- res
hit_table$hit <- is.finite(hit_table$specific_padj_by_reporter) &
  hit_table$specific_padj_by_reporter < 0.05
hit_table <- hit_table[hit_table$hit, , drop = FALSE]
hit_table <- hit_table[order(hit_table$specific_padj_by_reporter, -abs(hit_table$specific_effect)), , drop = FALSE]
hit_table <- utils::head(hit_table, 24)
hit_table$pair <- paste(hit_table$reporter, hit_table$perturbation, sep = " / ")
hit_table$pair <- factor(hit_table$pair, levels = rev(hit_table$pair))
p_hits <- ggplot2$ggplot(
  hit_table,
  ggplot2$aes(specific_effect, pair)
) +
  ggplot2$geom_vline(xintercept = 0, color = "#9CA3AF", linewidth = 0.3) +
  ggplot2$geom_point(ggplot2$aes(color = specific_effect), size = 2.2, alpha = 0.92) +
  ggplot2$scale_color_gradient2(
    low = "#3f6fb5",
    mid = "#9CA3AF",
    high = "#c93f3f",
    midpoint = 0,
    name = expression(hat(delta)[aj]^PA)
  ) +
  ggplot2$labs(
    x = expression(hat(delta)[aj]^PA),
    y = NULL
  ) +
  paper_theme(9) +
  ggplot2$theme(
    axis.text.y = ggplot2$element_text(size = 7),
    legend.position = "right",
    legend.background = ggplot2$element_blank(),
    legend.key = ggplot2$element_blank()
  ) +
  ggplot2$guides(
    color = ggplot2$guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barheight = grid::unit(18, "mm"),
      barwidth = grid::unit(3.5, "mm")
    )
  )

discovery_figure <- gridExtra$arrangeGrob(
  panel_label("a", p_hist),
  panel_label("b", p_volcano),
  panel_label("c", p_hits),
  ncol = 1,
  heights = c(0.7, 1.05, 1.0)
)

discovery_base <- file.path(in_dir, paste0("pfibox_pa_discovery_summary_", run_tag))
ggplot2$ggsave(paste0(discovery_base, ".png"), discovery_figure, width = 9.2, height = 10.2, dpi = 300, bg = "white")
ggplot2$ggsave(paste0(discovery_base, ".pdf"), discovery_figure, width = 9.2, height = 10.2, bg = "white")

utils::write.table(
  hit_table[order(hit_table$specific_padj_by_reporter), ],
  file.path(in_dir, paste0("pfibox_pa_top_specific_hits_", run_tag, ".tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Using PFIbox PA run tag: ", run_tag)
message("Time window T: ", t_label)
message("Wrote effect-matrix figure: ", paste0(effect_base, ".png"))
message("Wrote discovery figure: ", paste0(discovery_base, ".png"))
