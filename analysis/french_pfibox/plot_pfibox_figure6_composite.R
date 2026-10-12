#!/usr/bin/env Rscript

# Four-panel manuscript figure for the PFIbox DStressR analysis.
# The figure combines the DStressR total-effect clustering, the mapping to
# French et al. time-course clusters, and the corresponding inferential
# diagnostics.

if (!requireNamespace("gridExtra", quietly = TRUE)) {
  stop("Package `gridExtra` is required for the composite PFIbox figure.", call. = FALSE)
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package `ggplot2` is required for the composite PFIbox figure.", call. = FALSE)
}

out_dir <- file.path(
  "analysis", "outputs", "french_pfibox", "promoter_activity_subset"
)

Sys.setenv(
  PFIBOX_CLUSTER_TOP_REPORTERS = "2000",
  PFIBOX_CLUSTER_LINKAGE = "ward.D2",
  PFIBOX_CLUSTER_DISTANCE = "euclidean",
  PFIBOX_CLUSTER_OUTPUT_SUFFIX = "figure6_ward",
  PFIBOX_VOLCANO_LABELS = "18"
)

cluster_env <- new.env(parent = globalenv())
source(
  file.path("analysis", "french_pfibox", "plot_pfibox_model_cluster_and_discovery.R"),
  local = cluster_env
)

alluvial_env <- new.env(parent = globalenv())
source(
  file.path("analysis", "french_pfibox", "plot_french_cluster_alluvial.R"),
  local = alluvial_env
)

panel_label <- function(label, grob) {
  gridExtra::arrangeGrob(
    grid::textGrob(
      label,
      x = grid::unit(0, "npc"),
      y = grid::unit(0.06, "npc"),
      hjust = 0,
      vjust = 0,
      gp = grid::gpar(fontsize = 15, fontface = "bold", col = "#111827")
    ),
    grob,
    ncol = 1,
    heights = c(0.060, 1)
  )
}

strip_facet <- function(plot) {
  plot +
    ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_blank()
    )
}

reporter_map <- alluvial_env$reporter_cluster_order[, c("reporter_cluster", "display_cluster")]
names(reporter_map) <- c("raw_cluster", "display_cluster")
reporter_map$raw_cluster <- as.character(reporter_map$raw_cluster)
reporter_map$display_cluster <- as.character(reporter_map$display_cluster)

row_df <- data.frame(
  reporter = cluster_env$row_order,
  row_top = seq_along(cluster_env$row_order),
  raw_cluster = paste0("R", cluster_env$row_cluster[cluster_env$row_order]),
  stringsAsFactors = FALSE
)
row_df$y <- length(cluster_env$row_order) - row_df$row_top + 1
row_df <- merge(row_df, reporter_map, by = "raw_cluster", all.x = TRUE, sort = FALSE)
row_df$display_cluster <- factor(
  row_df$display_cluster,
  levels = alluvial_env$reporter_cluster_order$display_cluster
)

cluster_blocks <- stats::aggregate(
  y ~ display_cluster,
  row_df,
  function(z) c(ymin = min(z) - 0.5, ymax = max(z) + 0.5, ymid = mean(range(z)))
)
cluster_blocks <- data.frame(
  display_cluster = cluster_blocks$display_cluster,
  ymin = cluster_blocks$y[, "ymin"],
  ymax = cluster_blocks$y[, "ymax"],
  ymid = cluster_blocks$y[, "ymid"],
  stringsAsFactors = FALSE
)
cluster_blocks$display_cluster <- factor(
  cluster_blocks$display_cluster,
  levels = alluvial_env$reporter_cluster_order$display_cluster
)
cluster_blocks <- cluster_blocks[order(cluster_blocks$display_cluster), ]
reporter_cols <- alluvial_env$reporter_cols
cluster_blocks$fill_color <- reporter_cols[as.character(cluster_blocks$display_cluster)]

plot_long <- merge(
  cluster_env$plot_long,
  row_df[, c("reporter", "display_cluster")],
  by = "reporter",
  all.x = TRUE,
  sort = FALSE
)
plot_long$reporter <- as.character(plot_long$reporter)
plot_long$perturbation <- as.character(plot_long$perturbation)

cluster_label_df <- merge(
  cluster_blocks,
  alluvial_env$reporter_cluster_order[, c("display_cluster", "total_n")],
  by = "display_cluster",
  all.x = TRUE,
  sort = FALSE
)
cluster_label_df$label <- as.character(cluster_label_df$display_cluster)

p_cluster <- ggplot2::ggplot(plot_long, ggplot2::aes(col_index, row_index, fill = effect)) +
  ggplot2::geom_tile(width = 1, height = 1, color = NA) +
  ggplot2::geom_rect(
    data = cluster_blocks,
    ggplot2::aes(
      xmin = -1.10,
      xmax = -0.25,
      ymin = ymin,
      ymax = ymax
    ),
    inherit.aes = FALSE,
    fill = cluster_blocks$fill_color,
    color = "white",
    linewidth = 0.45
  ) +
  ggplot2::geom_text(
    data = cluster_label_df,
    ggplot2::aes(x = -0.675, y = ymid, label = label),
    inherit.aes = FALSE,
    hjust = 0.5,
    size = 2.55,
    fontface = "bold",
    color = "white"
  ) +
  ggplot2::geom_hline(
    yintercept = cluster_blocks$ymin[-1],
    color = "#374151",
    linewidth = 0.20,
    alpha = 0.62
  ) +
  ggplot2::geom_vline(
    xintercept = cluster_env$col_breaks,
    color = "white",
    linewidth = 0.22
  ) +
  ggplot2::scale_fill_gradient2(
    low = "#3f6fb5",
    mid = "white",
    high = "#c93f3f",
    midpoint = 0,
    limits = c(-cluster_env$clip_limit, cluster_env$clip_limit),
    oob = scales::squish,
    na.value = "#F3F4F6",
    name = expression(hat(Delta * r)[aj]^"tot")
  ) +
  ggplot2::scale_x_continuous(
    breaks = seq_along(cluster_env$col_order),
    labels = cluster_env$col_order,
    limits = c(-1.45, length(cluster_env$col_order) + 0.5),
    expand = c(0, 0)
  ) +
  ggplot2::scale_y_continuous(
    limits = c(0.5, length(cluster_env$row_order) + 0.5),
    expand = c(0, 0)
  ) +
  ggplot2::labs(x = NULL, y = "Reporters clustered") +
  ggplot2::theme_minimal(base_size = 10) +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(angle = 58, hjust = 1, vjust = 1, size = 6.4),
    axis.text.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_text(size = 8.8, color = "#111827"),
    legend.position = "right",
    legend.title = ggplot2::element_text(size = 8.5),
    legend.text = ggplot2::element_text(size = 7.6),
    plot.margin = ggplot2::margin(0, 5, 0, 0)
  ) +
  ggplot2::guides(
    fill = ggplot2::guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barheight = grid::unit(30, "mm"),
      barwidth = grid::unit(4, "mm")
    )
  )

p_alluvial <- alluvial_env$p +
  ggplot2::theme(
    plot.margin = ggplot2::margin(-2, 2, 44, 2)
  )

p_hist <- cluster_env$p_hist +
  ggplot2::theme(
    axis.title = ggplot2::element_text(size = 9),
    axis.text = ggplot2::element_text(size = 8),
    plot.margin = ggplot2::margin(2, 6, 2, 2)
  )

p_volcano <- strip_facet(cluster_env$p_volcano) +
  ggplot2::theme(
    axis.title = ggplot2::element_text(size = 9),
    axis.text = ggplot2::element_text(size = 8),
    plot.margin = ggplot2::margin(2, 2, 2, 6)
  )

top_row <- gridExtra::arrangeGrob(
  panel_label("a", p_cluster),
  panel_label("b", p_alluvial),
  ncol = 2,
  widths = c(1.18, 1)
)

bottom_row <- gridExtra::arrangeGrob(
  panel_label("c", p_hist),
  panel_label("d", p_volcano),
  ncol = 2,
  widths = c(0.82, 1.18)
)

figure6 <- gridExtra::arrangeGrob(
  top_row,
  bottom_row,
  ncol = 1,
  heights = c(1.48, 1)
)

png_path <- file.path(out_dir, "pfibox_figure6_composite.png")
pdf_path <- file.path(out_dir, "pfibox_figure6_composite.pdf")

ggplot2::ggsave(png_path, figure6, width = 13.6, height = 11.2, dpi = 320, bg = "white")
ggplot2::ggsave(pdf_path, figure6, width = 13.6, height = 11.2, device = grDevices::pdf, bg = "white")

message("Wrote composite Figure 6:")
message("  ", png_path)
message("  ", pdf_path)
