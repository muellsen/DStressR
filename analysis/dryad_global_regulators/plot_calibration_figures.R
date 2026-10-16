source(file.path("analysis", "_helpers.R"))
source(file.path("analysis", "_paper_artifacts.R"))
load_destress_package()
set.seed(1) # Reproducible label placement only.
ggplot2 <- asNamespace("ggplot2")
out_dir <- analysis_output_dir("dryad_global_regulators")
artifact <- read_paper_artifact(file.path(out_dir, "dash_calibration_paper.rds"))
res <- artifact$results$default
invisible(list2env(artifact$parameters, environment()))
invisible(list2env(artifact$annotations, environment()))
alpha_promoter_long <- rbind(
  data.frame(
    promoter = alpha_fit$reporter,
    panel = "Raw growth exponent",
    estimate = alpha_fit$alpha_raw,
    se = alpha_fit$alpha_raw_se,
    control_n = alpha_fit$control_n,
    stringsAsFactors = FALSE
  ),
  data.frame(
    promoter = alpha_fit$reporter,
    panel = "Shrunken growth exponent",
    estimate = alpha_fit$alpha_shrunk,
    se = alpha_fit$alpha_shrunk_se,
    control_n = alpha_fit$control_n,
    stringsAsFactors = FALSE
  )
)
alpha_promoter_long$panel <- factor(
  alpha_promoter_long$panel,
  levels = c("Raw growth exponent", "Shrunken growth exponent")
)
promoter_alpha_order <- alpha_fit[order(alpha_fit$alpha_shrunk), "reporter"]
alpha_promoter_long$promoter <- factor(alpha_promoter_long$promoter, levels = promoter_alpha_order)
alpha_global <- unique(alpha_fit$alpha_global)
alpha_global <- alpha_global[is.finite(alpha_global)][1]

p_promoter_alpha <- ggplot2::ggplot(
  alpha_promoter_long,
  ggplot2::aes(x = estimate, y = promoter)
) +
  ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "#94A3B8", linewidth = 0.35) +
  ggplot2::geom_vline(xintercept = alpha_global, linetype = "dotted", color = "#475569", linewidth = 0.4) +
  ggplot2::geom_errorbar(
    ggplot2::aes(xmin = estimate - 1.96 * se, xmax = estimate + 1.96 * se),
    width = 0,
    color = "#64748B",
    linewidth = 0.45
  ) +
  ggplot2::geom_point(color = "#1D4ED8", fill = "#93C5FD", size = 2.2) +
  ggplot2::facet_wrap(
    ggplot2::vars(panel),
    nrow = 1
  ) +
  ggplot2::theme_light(base_size = 10) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    strip.background = ggplot2::element_rect(fill = "#F1F5F9", color = "#CBD5E1"),
    strip.text = ggplot2::element_text(color = "#111827"),
    legend.position = "none"
  ) +
  ggplot2::labs(
    x = expression("Growth-response exponent estimate " * hat(alpha)[a]),
    y = "Reporter"
  )
ggplot2::ggsave(file.path(out_dir, "dryad_growth_transition_matching_windows_alpha_shrinkage_estimates.png"), p_promoter_alpha, width = 8, height = 5.8, dpi = 300)
ggplot2::ggsave(file.path(out_dir, "dryad_growth_transition_matching_windows_alpha_shrinkage_estimates.pdf"), p_promoter_alpha, width = 8, height = 5.8)

alpha_window_long <- rbind(
  data.frame(
    base_promoter = alpha_window_fit$base_promoter,
    window = alpha_window_fit$window,
    estimate_type = "Raw",
    alpha = alpha_window_fit$alpha_raw,
    se = alpha_window_fit$alpha_raw_se,
    stringsAsFactors = FALSE
  ),
  data.frame(
    base_promoter = alpha_window_fit$base_promoter,
    window = alpha_window_fit$window,
    estimate_type = "Shrunken",
    alpha = alpha_window_fit$alpha_shrunk,
    se = alpha_window_fit$alpha_shrunk_se,
    stringsAsFactors = FALSE
  )
)
alpha_window_long$estimate_type <- factor(alpha_window_long$estimate_type, levels = c("Raw", "Shrunken"))
reporter_order <- stats::aggregate(alpha_shrunk ~ base_promoter, alpha_window_fit, median)
reporter_order <- reporter_order[order(reporter_order$alpha_shrunk), "base_promoter"]
alpha_window_long$base_promoter <- factor(alpha_window_long$base_promoter, levels = reporter_order)
alpha_limit <- max(abs(alpha_window_long$alpha), na.rm = TRUE)

p_alpha_heatmap <- ggplot2::ggplot(
  alpha_window_long,
  ggplot2::aes(window, base_promoter, fill = pmax(pmin(alpha, alpha_limit), -alpha_limit))
) +
  ggplot2::geom_tile(color = "white", linewidth = 0.3) +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", alpha)), size = 2.6) +
  ggplot2::facet_wrap(ggplot2::vars(estimate_type), nrow = 1) +
  ggplot2::scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-alpha_limit, alpha_limit)
  ) +
  ggplot2::theme_light(base_size = 9) +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    strip.background = ggplot2::element_rect(fill = "grey95", color = "grey75"),
    strip.text = ggplot2::element_text(face = "bold", size = 9, color = "grey15")
  ) +
  ggplot2::labs(
    x = "Time window",
    y = "Reporter",
    fill = expression(hat(alpha)[a])
  )
ggplot2::ggsave(file.path(out_dir, "dryad_growth_transition_matching_windows_alpha_heatmap.png"), p_alpha_heatmap, width = 7.5, height = 7, dpi = 300)
ggplot2::ggsave(file.path(out_dir, "dryad_growth_transition_matching_windows_alpha_heatmap.pdf"), p_alpha_heatmap, width = 7.5, height = 7)

alpha_curve <- alpha_window_fit
alpha_curve$base_promoter <- factor(alpha_curve$base_promoter, levels = reporter_order)
p_alpha_curves <- ggplot2::ggplot(alpha_curve, ggplot2::aes(window, alpha_raw, group = base_promoter)) +
  ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "grey55", linewidth = 0.25) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = alpha_raw - alpha_raw_se, ymax = alpha_raw + alpha_raw_se),
    width = 0.12,
    color = "grey65",
    linewidth = 0.25
  ) +
  ggplot2::geom_point(color = "grey30", size = 1.4) +
  ggplot2::geom_line(ggplot2::aes(y = alpha_shrunk), color = "#009E73", linewidth = 0.55) +
  ggplot2::geom_point(ggplot2::aes(y = alpha_shrunk), color = "#009E73", size = 1.4) +
  ggplot2::facet_wrap(ggplot2::vars(base_promoter), ncol = 3) +
  ggplot2::coord_cartesian(ylim = stats::quantile(c(alpha_curve$alpha_raw, alpha_curve$alpha_shrunk), c(0.02, 0.98), na.rm = TRUE)) +
  ggplot2::theme_light(base_size = 8) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(angle = 30, hjust = 1),
    strip.background = ggplot2::element_rect(fill = "grey95", color = "grey75"),
    strip.text = ggplot2::element_text(face = "bold", size = 8, color = "grey15")
  ) +
  ggplot2::labs(
    x = "Time window",
    y = expression("Growth-response exponent estimate " * hat(alpha)[a])
  )
ggplot2::ggsave(file.path(out_dir, "dryad_growth_transition_matching_windows_alpha_curves.png"), p_alpha_curves, width = 8, height = 8, dpi = 300)
ggplot2::ggsave(file.path(out_dir, "dryad_growth_transition_matching_windows_alpha_curves.pdf"), p_alpha_curves, width = 8, height = 8)

plot_window_heatmap <- function(tab, effect_col, fill_label) {
  d <- tab
  d$compound <- factor(d$compound, levels = c("Iron", "Tetracycline", "H2O2", "Kanamycin"))
  d$pseudo_reporter <- factor(d$promoter, levels = rev(pseudo_reporter_levels))
  limit <- max(abs(d[[effect_col]]), na.rm = TRUE)
  ggplot2::ggplot(d, ggplot2::aes(compound, pseudo_reporter, fill = .data[[effect_col]])) +
    ggplot2::geom_tile(color = "white", linewidth = 0.35) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", .data[[effect_col]])), size = 2.4) +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-limit, limit)
    ) +
    ggplot2::theme_light(base_size = 9) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 30, hjust = 1)
    ) +
    ggplot2::labs(x = "Weak-stress condition", y = "Reporter-window unit", fill = fill_label)
}

p_total <- plot_window_heatmap(
  res,
  "total_effect",
  "Estimated total effect"
)
ggplot2::ggsave(file.path(out_dir, "dryad_weak_stress_windows_calibrated_alpha_total_effect_heatmap.png"), p_total, width = 7.5, height = 6, dpi = 300)
ggplot2::ggsave(file.path(out_dir, "dryad_weak_stress_windows_calibrated_alpha_total_effect_heatmap.pdf"), p_total, width = 7.5, height = 6)

p_specific <- plot_window_heatmap(
  res,
  "specific_effect",
  "Estimated specific effect"
)
ggplot2::ggsave(file.path(out_dir, "dryad_weak_stress_windows_calibrated_alpha_specific_effect_heatmap.png"), p_specific, width = 7.5, height = 6, dpi = 300)
ggplot2::ggsave(file.path(out_dir, "dryad_weak_stress_windows_calibrated_alpha_specific_effect_heatmap.pdf"), p_specific, width = 7.5, height = 6)

p_hist <- ggplot2::ggplot(res, ggplot2::aes(x = specific_pvalue)) +
  ggplot2::geom_histogram(breaks = seq(0, 1, by = 0.05), fill = "#009E73", color = "white", linewidth = 0.25) +
  ggplot2::theme_light(base_size = 10) +
  ggplot2::labs(
    x = "Raw p-value",
    y = "Number of tests"
  )
ggplot2::ggsave(file.path(out_dir, "dryad_weak_stress_windows_calibrated_alpha_pvalue_histogram_combined.png"), p_hist, width = 7, height = 4.5, dpi = 300)
ggplot2::ggsave(file.path(out_dir, "dryad_weak_stress_windows_calibrated_alpha_pvalue_histogram_combined.pdf"), p_hist, width = 7, height = 4.5)

p_volcano <- plot_volcano(
  res,
  effect = "specific_effect",
  padj = "specific_padj_by_reporter",
  title = NULL,
  label_by = "pair",
  top_n = 15,
  top_reporters = 8,
  xlab = "Estimated specific effect",
  ylab = "-log10 adjusted p-value"
)
ggplot2::ggsave(file.path(out_dir, "dryad_weak_stress_windows_calibrated_alpha_volcano.png"), p_volcano, width = 8, height = 5.5, dpi = 300)
ggplot2::ggsave(file.path(out_dir, "dryad_weak_stress_windows_calibrated_alpha_volcano.pdf"), p_volcano, width = 8, height = 5.5)

