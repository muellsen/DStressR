source(file.path("analysis", "_helpers.R"))
load_destress_package()

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package `ggplot2` is required for Binsfeld mean-variance diagnostics.", call. = FALSE)
}
ggplot2 <- asNamespace("ggplot2")

out_dir <- analysis_output_dir("binsfeld_variance_diagnostics")
diagnostics_file <- file.path(out_dir, "ecoli_raw_vs_moderated_variance.tsv")
if (!file.exists(diagnostics_file)) {
  stop(
    "Missing variance diagnostics table: ", diagnostics_file,
    "\nRun analysis/ecoli_promoter_screen/plot_variance_moderation_diagnostics.R first.",
    call. = FALSE
  )
}

diagnostics <- read_tsv_base(diagnostics_file)
required <- c("workflow", "promoter", "compound", "mean_response", "effect", "moderated_variance")
missing_cols <- setdiff(required, names(diagnostics))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
}

diagnostics <- diagnostics[
  diagnostics$workflow %in% c("DStressR without EV control", "DStressR with EV control") &
    is.finite(diagnostics$mean_response) &
    is.finite(diagnostics$effect) &
    is.finite(diagnostics$moderated_variance),
  ,
  drop = FALSE
]
diagnostics$workflow <- sub("without EV control", "without EVC", diagnostics$workflow, fixed = TRUE)
diagnostics$workflow <- sub("with EV control", "with EVC", diagnostics$workflow, fixed = TRUE)

make_workflow_diag <- function(d, variance_effect) {
  out <- do.call(rbind, lapply(split(d, d$workflow), function(wd) {
    tab <- perturbation_diagnostics(
      wd,
      mean_effect = "mean_response",
      variance_effect = variance_effect,
      reporter = "promoter",
      perturbation = "compound",
      trend_span = 0.45
    )
    tab$workflow <- unique(wd$workflow)
    tab
  }))
  rownames(out) <- NULL
  out
}

total_diag <- make_workflow_diag(diagnostics, variance_effect = "mean_response")
specific_diag <- make_workflow_diag(diagnostics, variance_effect = "effect")

model_variance <- stats::aggregate(
  moderated_variance ~ workflow + compound,
  diagnostics,
  function(x) stats::median(x, na.rm = TRUE)
)
names(model_variance)[3] <- "median_moderated_model_variance"

compound_diag <- merge(
  total_diag[, c(
    "workflow", "perturbation", "mean_effect", "abs_mean_effect",
    "effect_variance", "n_reporters", "rank_abs_mean_effect",
    "variance_trend", "variance_residual", "rank_variance_residual"
  )],
  specific_diag[, c("workflow", "perturbation", "effect_variance")],
  by = c("workflow", "perturbation"),
  all = TRUE,
  sort = FALSE
)
names(compound_diag) <- c(
  "workflow", "compound", "mean_total_response", "abs_mean_total_response",
  "total_effect_variance", "n_promoters", "rank_abs_mean_total_response",
  "variance_trend", "variance_residual", "rank_variance_residual",
  "specific_effect_variance"
)
compound_diag <- merge(compound_diag, model_variance, by = c("workflow", "compound"), all.x = TRUE, sort = FALSE)
compound_diag <- compound_diag[
  order(compound_diag$workflow, compound_diag$rank_abs_mean_total_response),
  ,
  drop = FALSE
]

utils::write.table(
  compound_diag,
  file.path(out_dir, "ecoli_compound_mean_specific_variance_diagnostic.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

total_diag$workflow <- factor(
  total_diag$workflow,
  levels = c("DStressR without EVC", "DStressR with EVC")
)
finite_color <- abs(total_diag$mean_effect[is.finite(total_diag$mean_effect)])
color_max <- stats::quantile(finite_color, probs = 0.92, names = FALSE, na.rm = TRUE)
if (!is.finite(color_max) || color_max <= 0) {
  color_max <- max(finite_color, na.rm = TRUE)
}
if (!is.finite(color_max) || color_max <= 0) {
  color_max <- 1
}
total_diag$color_effect <- pmax(pmin(total_diag$mean_effect, color_max), -color_max)
label_data <- do.call(rbind, lapply(split(total_diag, total_diag$workflow), function(d) {
  d <- d[order(-d$variance_residual, -d$effect_variance), , drop = FALSE]
  d <- d[is.finite(d$variance_residual) | is.finite(d$effect_variance), , drop = FALSE]
  utils::head(d, 7)
}))

mean_variance_plot <- ggplot2$ggplot(
  total_diag,
  ggplot2$aes(rank_abs_mean_effect, effect_variance)
) +
  ggplot2$geom_point(
    ggplot2$aes(color = color_effect),
    size = 1.65,
    alpha = 0.82
  ) +
  ggplot2$geom_line(
    data = total_diag[is.finite(total_diag$variance_trend), , drop = FALSE],
    ggplot2$aes(y = variance_trend),
    color = "#008A5B",
    linewidth = 0.65
  ) +
  ggplot2$facet_wrap(ggplot2$vars(workflow), nrow = 1) +
  ggplot2$scale_color_gradient2(
    low = "#1D4ED8",
    mid = "#9CA3AF",
    high = "#B91C1C",
    midpoint = 0,
    limits = c(-color_max, color_max),
    name = expression(hat(Delta * bar(y))[j]^"tot")
  ) +
  ggplot2$theme_light(base_size = 10) +
  ggplot2$theme(
    panel.grid.minor = ggplot2$element_blank(),
    panel.grid.major = ggplot2$element_line(linewidth = 0.25, color = "#E5E7EB"),
    strip.text = ggplot2$element_text(face = "bold", color = "#111827"),
    axis.title = ggplot2$element_text(color = "#111827"),
    axis.text = ggplot2$element_text(color = "#4B5563"),
    legend.position = "right",
    legend.justification = "center",
    legend.background = ggplot2$element_blank(),
    legend.key.height = grid::unit(14, "mm"),
    legend.key.width = grid::unit(3, "mm"),
    legend.title = ggplot2$element_text(size = 8),
    legend.text = ggplot2$element_text(size = 7),
    plot.title = ggplot2$element_blank(),
    plot.subtitle = ggplot2$element_blank()
  ) +
  ggplot2$guides(
    color = ggplot2$guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barwidth = grid::unit(3, "mm"),
      barheight = grid::unit(18, "mm")
    )
  ) +
  ggplot2$labs(
    x = expression("Rank of " * abs(hat(Delta * bar(y))[j]^"tot")),
    y = expression("Variance across promoters of " * hat(Delta * y)["\u00b7" * j]^"tot")
  )

if (nrow(label_data) > 0 && requireNamespace("ggrepel", quietly = TRUE)) {
  mean_variance_plot <- mean_variance_plot +
    ggrepel::geom_text_repel(
      data = label_data,
      ggplot2$aes(label = perturbation_label),
      size = 2.35,
      min.segment.length = 0,
      box.padding = 0.2,
      point.padding = 0.1,
      max.overlaps = Inf,
      seed = 1,
      show.legend = FALSE
    )
}

ggplot2$ggsave(
  file.path(out_dir, "dstressr_rank_variance_diagnostic.png"),
  mean_variance_plot,
  width = 9.2,
  height = 4.6,
  dpi = 320,
  bg = "white"
)
ggplot2$ggsave(
  file.path(out_dir, "dstressr_rank_variance_diagnostic.pdf"),
  mean_variance_plot,
  width = 9.2,
  height = 4.6
)

fit_summary <- fit_variance_distribution(
  compound_diag,
  variance = "specific_effect_variance"
)
utils::write.table(
  fit_summary,
  file.path(out_dir, "ecoli_compound_variance_distribution_fits.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

variance_density_plot <- plot_variance_distribution(
  compound_diag,
  fit = fit_summary,
  variance = "specific_effect_variance",
  title = NULL,
  subtitle = ""
)

ggplot2$ggsave(
  file.path(out_dir, "ecoli_compound_variance_density.png"),
  variance_density_plot,
  width = 8.0,
  height = 4.2,
  dpi = 320,
  bg = "white"
)
ggplot2$ggsave(
  file.path(out_dir, "ecoli_compound_variance_density.pdf"),
  variance_density_plot,
  width = 8.0,
  height = 4.2
)

message("Wrote DStressR rank-variance diagnostic to: ", out_dir)
