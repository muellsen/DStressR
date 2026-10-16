source(file.path("analysis", "_helpers.R"))
load_destress_package()

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package `ggplot2` is required for effect matrix plots.", call. = FALSE)
}
ggplot2 <- asNamespace("ggplot2")
scales <- asNamespace("scales")

out_dir <- analysis_output_dir("binsfeld_effect_matrix_views")
source(file.path("analysis", "_paper_artifacts.R"))
artifact <- read_paper_artifact(file.path(out_dir, "binsfeld_effect_paper.rds"))
res <- artifact$results$with_EVC
res$promoter <- res$reporter
res$compound <- res$perturbation
shared_color_limit <- max(abs(as.numeric(res$specific_effect)), na.rm = TRUE)
style_matrix_plot <- function(plot) {
  plot +
    ggplot2$theme(
      legend.position = "right",
      legend.justification = "center",
      legend.background = ggplot2$element_blank(),
      legend.key.height = grid::unit(14, "mm"),
      legend.key.width = grid::unit(3, "mm"),
      legend.title = ggplot2$element_text(size = 7),
      legend.text = ggplot2$element_text(size = 6),
      legend.margin = ggplot2$margin(0, 0, 0, 4),
      legend.spacing = grid::unit(1, "mm"),
      axis.text.x = ggplot2$element_text(size = 6, margin = ggplot2$margin(t = 1)),
      axis.text.y = ggplot2$element_text(size = 8),
      axis.title.x = ggplot2$element_text(size = 9, margin = ggplot2$margin(t = 4))
    ) +
    ggplot2$guides(
      fill = ggplot2$guide_colorbar(
        title.position = "top",
        title.hjust = 0.5,
        label.position = "right",
        barwidth = grid::unit(3, "mm"),
        barheight = grid::unit(14, "mm")
      )
    )
}

plots <- list(
  raw_effect_matrix = plot_response_heatmap(
    res,
    value = "specific_effect",
    reporter = "promoter",
    perturbation = "compound",
    show_perturbation_ids = FALSE,
    top_n_perturbations = Inf,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    clip_quantile = 1,
    color_limit = shared_color_limit,
    title = NULL,
    subtitle = "",
    xlab = "Compounds",
    ylab = NULL,
    legend_title = "Specific effect"
  ),
  clustered_effect_matrix = plot_response_heatmap(
    res,
    value = "specific_effect",
    reporter = "promoter",
    perturbation = "compound",
    show_perturbation_ids = FALSE,
    top_n_perturbations = Inf,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    clip_quantile = 1,
    color_limit = shared_color_limit,
    title = NULL,
    subtitle = "",
    xlab = "Compounds",
    ylab = NULL,
    legend_title = "Specific effect"
  ),
  significant_hit_matrix = plot_hit_heatmap(
    res,
    effect = "specific_effect",
    padj = "specific_padj_by_reporter",
    reporter = "promoter",
    perturbation = "compound",
    show_perturbation_ids = FALSE,
    top_n_perturbations = Inf,
    order_rows = "global",
    order_cols = "input",
    clip_quantile = 1,
    color_limit = shared_color_limit,
    title = NULL,
    subtitle = "",
    xlab = "Compounds",
    ylab = NULL,
    legend_title = "Specific effect"
  ),
  clustered_significant_hit_matrix = plot_hit_heatmap(
    res,
    effect = "specific_effect",
    padj = "specific_padj_by_reporter",
    reporter = "promoter",
    perturbation = "compound",
    show_perturbation_ids = FALSE,
    top_n_perturbations = Inf,
    order_rows = "cluster",
    order_cols = "cluster",
    clip_quantile = 1,
    color_limit = shared_color_limit,
    title = NULL,
    subtitle = "",
    xlab = "Compounds",
    ylab = NULL,
    legend_title = "Specific effect"
  )
)
plots <- lapply(plots, style_matrix_plot)

for (name in names(plots)) {
  ggplot2$ggsave(
    file.path(out_dir, paste0("binsfeld_evc_specific_", name, ".png")),
    plots[[name]],
    width = 10.4,
    height = 4.4,
    dpi = 320,
    bg = "white"
  )
  ggplot2$ggsave(
    file.path(out_dir, paste0("binsfeld_evc_specific_", name, ".pdf")),
    plots[[name]],
    width = 10.4,
    height = 4.4
  )
}

message("Wrote optional exploratory matrix views.")
