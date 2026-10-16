source(file.path("analysis", "_helpers.R"))
source(file.path("analysis", "_paper_artifacts.R"))
load_destress_package()
set.seed(1) # Reproducible label placement only.

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package `ggplot2` is required for rank-1 adjustment plots.", call. = FALSE)
}
if (!requireNamespace("gridExtra", quietly = TRUE)) {
  stop("Package `gridExtra` is required for rank-1 adjustment plots.", call. = FALSE)
}

ggplot2 <- asNamespace("ggplot2")
gridExtra <- asNamespace("gridExtra")

panel_label <- function(label, size = 15) {
  grid::textGrob(
    label,
    x = grid::unit(0, "npc"),
    y = grid::unit(0.08, "npc"),
    hjust = 0,
    vjust = 0,
    gp = grid::gpar(fontsize = size, fontface = "bold", col = "#111827")
  )
}
out_dir <- analysis_output_dir("dryad_global_regulators")
out_dir <- analysis_output_dir("dryad_global_regulators")
artifact <- read_paper_artifact(file.path(out_dir, "dash_rank_paper.rds"))
res_rank0 <- artifact$results$rank0
res_rank1 <- artifact$results$rank1
rank_results <- rbind(res_rank0, res_rank1)
rank_diagnostics <- artifact$diagnostics$rank_diagnostics
invisible(list2env(artifact$annotations, environment()))
p_destress_rank_scree <- plot_background_rank_diagnostics(
  diagnostics = rank_diagnostics,
  xlab = "Component",
  ylab = "Singular value"
)
ggplot2::ggsave(
  file.path(out_dir, "dryad_weak_stress_windows_calibrated_alpha_destress_rank_scree.png"),
  p_destress_rank_scree,
  width = 5.2,
  height = 3.8,
  dpi = 300,
  bg = "white"
)
ggplot2::ggsave(
  file.path(out_dir, "dryad_weak_stress_windows_calibrated_alpha_destress_rank_scree.pdf"),
  p_destress_rank_scree,
  width = 5.2,
  height = 3.8
)

effect_long <- rbind(
  data.frame(
    res_rank1[, c("promoter", "base_promoter", "window", "compound")],
    effect_type = "Total effect",
    effect = res_rank1$total_effect,
    stringsAsFactors = FALSE
  ),
  data.frame(
    res_rank1[, c("promoter", "base_promoter", "window", "compound")],
    effect_type = "Specific effect",
    effect = res_rank0$specific_effect,
    stringsAsFactors = FALSE
  ),
  data.frame(
    res_rank1[, c("promoter", "base_promoter", "window", "compound")],
    effect_type = "Rank-1 component",
    effect = res_rank1$low_rank_effect,
    stringsAsFactors = FALSE
  ),
  data.frame(
    res_rank1[, c("promoter", "base_promoter", "window", "compound")],
    effect_type = "Rank-adjusted total effect",
    effect = res_rank1$rank_adjusted_total_effect,
    stringsAsFactors = FALSE
  )
)
effect_long$effect_type <- factor(
  effect_long$effect_type,
  levels = c("Total effect", "Specific effect", "Rank-1 component", "Rank-adjusted total effect")
)
effect_long$compound <- factor(effect_long$compound, levels = compound_levels)
effect_long$promoter <- factor(effect_long$promoter, levels = rev(pseudo_reporter_levels))
effect_limit <- max(abs(effect_long$effect), na.rm = TRUE)
expected_boxes <- merge(
  expand.grid(
    base_promoter = promoter_levels,
    window = window_levels,
    effect_type = c("Total effect", "Specific effect", "Rank-adjusted total effect"),
    stringsAsFactors = FALSE
  ),
  expected[, c("base_promoter", "compound")],
  by = "base_promoter",
  all.x = TRUE,
  sort = FALSE
)
expected_boxes$promoter <- factor(
  paste(expected_boxes$base_promoter, expected_boxes$window, sep = " | "),
  levels = rev(pseudo_reporter_levels)
)
expected_boxes$compound <- factor(expected_boxes$compound, levels = compound_levels)
expected_boxes$effect_type <- factor(expected_boxes$effect_type, levels = levels(effect_long$effect_type))

p_rank_heatmaps <- ggplot2::ggplot(
  effect_long,
  ggplot2::aes(compound, promoter, fill = pmax(pmin(effect, effect_limit), -effect_limit))
) +
  ggplot2::geom_tile(color = "white", linewidth = 0.35) +
  ggplot2::geom_tile(
    data = expected_boxes,
    ggplot2::aes(compound, promoter),
    inherit.aes = FALSE,
    fill = NA,
    color = "black",
    linewidth = 0.85
  ) +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", effect)), size = 2.2, color = "white") +
  ggplot2::facet_wrap(ggplot2::vars(effect_type), nrow = 1) +
  ggplot2::scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-effect_limit, effect_limit)
  ) +
  ggplot2::theme_light(base_size = 8) +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(angle = 30, hjust = 1),
    strip.background = ggplot2::element_rect(fill = "grey95", color = "grey75"),
    strip.text = ggplot2::element_text(face = "bold", size = 8, color = "grey15")
  ) +
  ggplot2::labs(
    x = "Weak-stress condition",
    y = "Reporter-window unit",
    fill = "Estimated effect"
  )
p_rank_heatmaps_labeled <- gridExtra$arrangeGrob(
  gridExtra$arrangeGrob(
    grid::nullGrob(),
    panel_label("a"),
    panel_label("b"),
    panel_label("c"),
    panel_label("d"),
    grid::nullGrob(),
    ncol = 6,
    widths = c(0.24, 1, 1, 1, 1, 0.30)
  ),
  p_rank_heatmaps,
  ncol = 1,
  heights = c(0.06, 1)
)
ggplot2::ggsave(
  file.path(out_dir, "dryad_weak_stress_windows_calibrated_alpha_rank1_effect_decomposition.png"),
  p_rank_heatmaps_labeled,
  width = 13,
  height = 6.0,
  dpi = 300
)
ggplot2::ggsave(
  file.path(out_dir, "dryad_weak_stress_windows_calibrated_alpha_rank1_effect_decomposition.pdf"),
  p_rank_heatmaps_labeled,
  width = 13,
  height = 6.0
)

p_hist <- ggplot2::ggplot(
  rank_results,
  ggplot2::aes(x = rank_target_pvalue, fill = factor(background_rank))
) +
  ggplot2::geom_histogram(
    breaks = seq(0, 1, by = 0.05),
    color = "white",
    linewidth = 0.25
  ) +
  ggplot2::facet_wrap(
    ggplot2::vars(background_rank),
    ncol = 1,
    labeller = ggplot2::as_labeller(c(
      "0" = "Default specific effect",
      "1" = "Rank-adjusted total effect"
    ))
  ) +
  ggplot2::scale_fill_manual(values = c("0" = "#009E73", "1" = "#0072B2"), guide = "none") +
  ggplot2::theme_light(base_size = 10) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    strip.background = ggplot2::element_rect(fill = "grey95", color = "grey75"),
    strip.text = ggplot2::element_text(face = "bold", size = 9, color = "grey15")
  ) +
  ggplot2::labs(
    x = "Raw p-value",
    y = "Number of tests"
  )
ggplot2::ggsave(
  file.path(out_dir, "dryad_weak_stress_windows_calibrated_alpha_rank1_pvalue_histograms.png"),
  p_hist,
  width = 7,
  height = 6,
  dpi = 300
)
ggplot2::ggsave(
  file.path(out_dir, "dryad_weak_stress_windows_calibrated_alpha_rank1_pvalue_histograms.pdf"),
  p_hist,
  width = 7,
  height = 6
)

p_volcano <- plot_volcano(
  res_rank1,
  effect = "rank_adjusted_total_effect",
  padj = "rank_adjusted_total_padj_by_reporter",
  title = NULL,
  label_by = "pair",
  top_n = 15,
  top_reporters = 8,
  xlab = "Estimated rank-adjusted total effect",
  ylab = "-log10 adjusted p-value"
)
ggplot2::ggsave(
  file.path(out_dir, "dryad_weak_stress_windows_calibrated_alpha_rank1_volcano.png"),
  p_volcano,
  width = 8,
  height = 5.5,
  dpi = 300
)
ggplot2::ggsave(
  file.path(out_dir, "dryad_weak_stress_windows_calibrated_alpha_rank1_volcano.pdf"),
  p_volcano,
  width = 8,
  height = 5.5
)

message("Wrote Dryad rank-1 adjustment outputs to: ", out_dir)
