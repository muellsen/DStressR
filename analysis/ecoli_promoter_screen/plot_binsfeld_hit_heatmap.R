source(file.path("analysis", "_helpers.R"))
load_destress_package()

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package `ggplot2` is required for hit heatmap plots.", call. = FALSE)
}
ggplot2 <- asNamespace("ggplot2")

out_dir <- analysis_output_dir("binsfeld_hit_diagnostics")
result_file <- file.path(
  analysis_project_root(),
  "analysis", "outputs", "binsfeld_evc_calibrated", "evc_huber_pair_results.tsv"
)
if (!file.exists(result_file)) {
  stop("Missing Binsfeld EVC result table: ", result_file, call. = FALSE)
}

res <- read_tsv_base(result_file)
required <- c("promoter", "compound", "specific_effect", "specific_padj_by_reporter")
missing_cols <- setdiff(required, names(res))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
}

hit_summary <- summarize_hits(
  res,
  effect = "specific_effect",
  padj = "specific_padj_by_reporter",
  reporter = "promoter",
  perturbation = "compound",
  fdr = 0.05
)
utils::write.table(
  hit_summary$pairs,
  file.path(out_dir, "binsfeld_evc_hit_pairs.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
utils::write.table(
  hit_summary$promoters,
  file.path(out_dir, "binsfeld_evc_hit_summary_by_reporter.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
utils::write.table(
  hit_summary$compounds,
  file.path(out_dir, "binsfeld_evc_hit_summary_by_compound.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

base_hit_plot <- plot_hit_heatmap(
  res,
  effect = "specific_effect",
  padj = "specific_padj_by_reporter",
  reporter = "promoter",
  perturbation = "compound",
  show_perturbation_ids = FALSE,
  top_n_perturbations = Inf,
  order_cols = "frequency",
  title = NULL,
  subtitle = "",
  xlab = "Compounds",
  ylab = NULL,
  legend_title = "Specific effect"
)

hit_plot <- base_hit_plot +
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

ggplot2$ggsave(
  file.path(out_dir, "binsfeld_evc_specific_hit_heatmap.png"),
  hit_plot,
  width = 10.4,
  height = 4.4,
  dpi = 320,
  bg = "white"
)
ggplot2$ggsave(
  file.path(out_dir, "binsfeld_evc_specific_hit_heatmap.pdf"),
  hit_plot,
  width = 10.4,
  height = 4.4
)

message("Wrote Binsfeld hit diagnostics to: ", out_dir)
