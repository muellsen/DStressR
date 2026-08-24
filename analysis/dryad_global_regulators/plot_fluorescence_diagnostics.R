source(file.path("analysis", "_helpers.R"))
load_destress_package()

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package `ggplot2` is required for fluorescence diagnostics.", call. = FALSE)
}
if (!requireNamespace("patchwork", quietly = TRUE)) {
  stop("Package `patchwork` is required for fluorescence diagnostics.", call. = FALSE)
}

ggplot2 <- asNamespace("ggplot2")
patchwork <- asNamespace("patchwork")

out_dir <- analysis_output_dir("dryad_global_regulators")
res_file <- file.path(out_dir, "dryad_weak_stress_windows_calibrated_alpha_pair_results.tsv")
if (!file.exists(res_file)) {
  stop("Missing fluorescence result table: ", res_file, call. = FALSE)
}

res <- read_tsv_base(res_file)
required <- c("promoter", "compound", "total_effect", "specific_effect", "specific_pvalue")
missing_cols <- setdiff(required, names(res))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
}

res$compound <- factor(
  res$compound,
  levels = c("Iron", "Tetracycline", "H2O2", "Kanamycin")
)
res$promoter <- factor(
  res$promoter,
  levels = c(
    "Fur | Early", "Fur | Middle", "Fur | Late",
    "MarA | Early", "MarA | Middle", "MarA | Late",
    "SoxS | Early", "SoxS | Middle", "SoxS | Late",
    "LexA | Early", "LexA | Middle", "LexA | Late"
  )
)

compound_diag <- perturbation_diagnostics(
  res,
  mean_effect = "total_effect",
  variance_effect = "total_effect",
  reporter = "promoter",
  perturbation = "compound",
  perturbation_label = "compound",
  min_reporters = 2,
  trend_span = 0.75
)

reporter_diag <- perturbation_diagnostics(
  res,
  mean_effect = "total_effect",
  variance_effect = "specific_effect",
  reporter = "compound",
  perturbation = "promoter",
  perturbation_label = "promoter",
  min_reporters = 2,
  trend_span = 0.75
)

utils::write.table(
  compound_diag,
  file.path(out_dir, "dryad_fluorescence_compound_diagnostics.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
utils::write.table(
  reporter_diag,
  file.path(out_dir, "dryad_fluorescence_reporter_window_diagnostics.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

compound_plot <- plot_mean_variance_diagnostic(
  diagnostics = compound_diag,
  mean_effect = "total_effect",
  variance_effect = "total_effect",
  label_by = "abs_mean",
  top_n = 4,
  add_trend = FALSE,
  point_size = 2.2,
  title = NULL,
  subtitle = NULL,
  xlab = expression("Rank of " * abs(hat(Delta * bar(y))[j]^"tot")),
  ylab = expression("Variance across reporter-windows of " * hat(Delta * y)["\u00b7" * j]^"tot"),
  legend_position = "right"
)

reporter_plot <- plot_mean_variance_diagnostic(
  diagnostics = reporter_diag,
  mean_effect = "total_effect",
  variance_effect = "specific_effect",
  label_by = "variance",
  top_n = 6,
  add_trend = TRUE,
  point_size = 2.0,
  title = NULL,
  subtitle = NULL,
  xlab = expression("Rank of " * abs(hat(Delta * bar(y))[a]^"tot")),
  ylab = expression("Variance across perturbations of " * hat(Delta * y)[a * "\u00b7"]^"spec"),
  legend_position = "right"
)

p_hist <- plot_pvalue_histogram(
  res,
  pvalue = "specific_pvalue",
  promoter = "promoter",
  by = "all",
  bins = 20,
  title = NULL,
  subtitle = NULL,
  xlab = "Raw p-value"
) +
  ggplot2$theme(legend.position = "none")

p_volcano <- plot_volcano(
  res,
  effect = "specific_effect",
  pvalue = "specific_pvalue",
  padj = "specific_padj_by_promoter",
  promoter = "promoter",
  compound = "compound",
  title = NULL,
  subtitle = NULL,
  xlab = expression(hat(Delta * y)[aj]^"spec"),
  ylab = expression(-log[10] * " raw p-value"),
  top_n = 8,
  top_promoters = 12
) +
  ggplot2$theme(legend.position = "right")

panel <- (compound_plot + reporter_plot) / (p_hist + p_volcano) +
  patchwork$plot_layout(guides = "collect") &
  ggplot2$theme(
    legend.position = "right",
    plot.margin = ggplot2$margin(4, 4, 4, 4)
  )

ggplot2$ggsave(
  file.path(out_dir, "dryad_fluorescence_diagnostic_panel.png"),
  panel,
  width = 12.2,
  height = 8.4,
  dpi = 320,
  bg = "white"
)
ggplot2$ggsave(
  file.path(out_dir, "dryad_fluorescence_diagnostic_panel.pdf"),
  panel,
  width = 12.2,
  height = 8.4
)
ggplot2$ggsave(
  file.path(out_dir, "dryad_fluorescence_compound_mean_variance.png"),
  compound_plot,
  width = 6.2,
  height = 4.4,
  dpi = 320,
  bg = "white"
)
ggplot2$ggsave(
  file.path(out_dir, "dryad_fluorescence_reporter_window_mean_variance.png"),
  reporter_plot,
  width = 6.2,
  height = 4.4,
  dpi = 320,
  bg = "white"
)

message("Wrote fluorescence diagnostics to: ", out_dir)
