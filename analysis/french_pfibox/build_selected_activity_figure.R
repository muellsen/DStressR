source(file.path("analysis", "_helpers.R"))
source(file.path("analysis", "french_pfibox", "activity_plotting.R"))
parent <- analysis_output_dir("french_pfibox", "temporal_unscaled_integral_analysis")
out <- analysis_output_dir("french_pfibox", "temporal_unscaled_integral_analysis", "selected_four_figure6")
d <- readRDS(file.path(parent, "timecourses", "top10_each_direction_by_volcano_diverse_timecourses.rds"))
keys <- data.frame(direction = c("negative", "positive", "negative", "positive"),
  direction_rank = c(1, 2, 9, 9))
selected <- d$top_pairs[match(paste(keys$direction, keys$direction_rank),
  paste(d$top_pairs$direction, d$top_pairs$direction_rank)), ]
stopifnot(identical(as.character(selected$reporter), c("ompC", "recN", "msrB", "rmf")))
selected_file <- file.path(out, "selected_four_pairs.tsv")
write.table(selected, selected_file, sep = "\t", quote = FALSE, row.names = FALSE)
long <- d$trajectories[d$trajectories$rank %in% selected$rank, ]
labels <- paste(selected$reporter, selected$perturbation, sep = " | ")
long$panel <- factor(long$rank, levels = selected$rank, labels = labels)
annotations <- data.frame(panel = factor(labels, levels = labels),
  label = vapply(seq_len(nrow(selected)), function(i) {
    exponent <- floor(log10(selected$total_pvalue[i]))
    mantissa <- selected$total_pvalue[i] / 10^exponent
    sprintf("atop(hat(Delta*y)[aj]^plain(tot)==%.0f,p[aj]==%.2f%%*%%10^{%d})",
      selected$total_effect[i], mantissa, exponent)
  }, character(1)))
annotations$.panel <- annotations$panel
p <- plot_activity_trajectories(long) +
  ggplot2::geom_text(data = annotations, ggplot2::aes(x = Inf, y = Inf, label = label),
    inherit.aes = FALSE, parse = TRUE, hjust = 1.08, vjust = 1.2, size = 3.1) +
  ggplot2::labs(y = "Unscaled promoter activity (per h)")
grob <- ggplot2::ggplotGrob(p)
panel_file <- file.path(out, "selected_four_timecourse_panel.rds")
saveRDS(grob, panel_file)
ggplot2::ggsave(file.path(out, "selected_four_timecourses.png"), grob,
  width = 10, height = 7, dpi = 250, bg = "white")
Sys.setenv(PFIBOX_ANALYSIS_DIR = parent, PFIBOX_FIGURE_OUTPUT_DIR = out,
  PFIBOX_RUN_TAG = "reporters2000_perturbations100_Tall_available_treatment_replicate_scalable_temporal_unscaled_integral",
  PFIBOX_CLUSTER_LINKAGE = "complete", PFIBOX_CLUSTER_DISTANCE = "correlation",
  PFIBOX_CLUSTER_OUTPUT_SUFFIX = "total_effect_correlation_complete",
  PFIBOX_ADJUSTMENT_SCOPE = "global", PFIBOX_VOLCANO_LABEL_REGION = "near",
  PFIBOX_VOLCANO_LABEL_PAIRS = selected_file, PFIBOX_TIMECOURSE_PANEL_RDS = panel_file,
  PFIBOX_HEATMAP_MARK_PAIRS = selected_file, PFIBOX_RASTER_PDF = "1")
source(file.path("analysis", "french_pfibox", "plot_clustering_and_inference.R"),
  local = new.env(parent = globalenv()))
