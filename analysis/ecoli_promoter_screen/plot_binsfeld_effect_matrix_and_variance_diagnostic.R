source(file.path("analysis", "_helpers.R"))
load_destress_package()

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package `ggplot2` is required for effect-matrix diagnostics.", call. = FALSE)
}
if (!requireNamespace("gridExtra", quietly = TRUE)) {
  stop("Package `gridExtra` is required for combined figures.", call. = FALSE)
}
if (!requireNamespace("scales", quietly = TRUE)) {
  stop("Package `scales` is required for color clipping.", call. = FALSE)
}

ggplot2 <- asNamespace("ggplot2")
gridExtra <- asNamespace("gridExtra")
scales <- asNamespace("scales")

panel_label <- function(label, size = 14) {
  grid::textGrob(
    label,
    x = grid::unit(0, "npc"),
    y = grid::unit(0.06, "npc"),
    hjust = 0,
    vjust = 0,
    gp = grid::gpar(fontsize = size, fontface = "bold", col = "#111827")
  )
}

effect_dir <- analysis_output_dir("binsfeld_effect_matrix_views")
variance_dir <- analysis_output_dir("binsfeld_variance_diagnostics")
comparison_file <- file.path(
  analysis_project_root(),
  "analysis", "outputs", "binsfeld_evc_calibrated",
  "evc_huber_comparison_to_binsfeld_and_default.tsv"
)
diagnostics_file <- file.path(variance_dir, "ecoli_raw_vs_moderated_variance.tsv")

if (!file.exists(comparison_file)) {
  stop("Missing comparison table: ", comparison_file, call. = FALSE)
}
if (!file.exists(diagnostics_file)) {
  stop("Missing variance diagnostics table: ", diagnostics_file, call. = FALSE)
}

comparison <- read_tsv_base(comparison_file)
diagnostics <- read_tsv_base(diagnostics_file)
diagnostics$workflow <- sub("without EV control", "without EVC", diagnostics$workflow, fixed = TRUE)
diagnostics$workflow <- sub("with EV control", "with EVC", diagnostics$workflow, fixed = TRUE)

required_comparison <- c("promoter", "compound", "modeled_hit", "evc_huber_hit")
required_diagnostics <- c("workflow", "promoter", "compound", "mean_response")
missing_comparison <- setdiff(required_comparison, names(comparison))
missing_diagnostics <- setdiff(required_diagnostics, names(diagnostics))
if (length(missing_comparison) > 0) {
  stop("Missing required comparison columns: ", paste(missing_comparison, collapse = ", "), call. = FALSE)
}
if (length(missing_diagnostics) > 0) {
  stop("Missing required diagnostics columns: ", paste(missing_diagnostics, collapse = ", "), call. = FALSE)
}

reporter_order <- c("acrABp", "marRABp", "micFp", "ompFp", "robp", "soxSp", "tolCp")
total_no_evc <- diagnostics[
  diagnostics$workflow == "DStressR without EVC",
  c("promoter", "compound", "mean_response"),
  drop = FALSE
]
names(total_no_evc)[3] <- "modeled_total_effect"
total_evc <- diagnostics[
  diagnostics$workflow == "DStressR with EVC",
  c("promoter", "compound", "mean_response"),
  drop = FALSE
]
names(total_evc)[3] <- "evc_huber_total_effect"

comparison <- merge(comparison, total_no_evc, by = c("promoter", "compound"), all.x = TRUE, sort = FALSE)
comparison <- merge(comparison, total_evc, by = c("promoter", "compound"), all.x = TRUE, sort = FALSE)
comparison$hit_count <- as.integer(comparison$modeled_hit) + as.integer(comparison$evc_huber_hit)
comparison$abs_sum <- abs(comparison$modeled_total_effect) + abs(comparison$evc_huber_total_effect)

compound_score <- stats::aggregate(
  cbind(hit_count, abs_sum) ~ compound,
  data = comparison,
  FUN = sum,
  na.rm = TRUE
)
compound_score <- compound_score[
  order(-compound_score$hit_count, -compound_score$abs_sum, compound_score$compound),
  ,
  drop = FALSE
]
compound_order <- compound_score$compound

total_long <- rbind(
  data.frame(
    workflow = "DStressR without EVC",
    promoter = comparison$promoter,
    compound = comparison$compound,
    effect = comparison$modeled_total_effect,
    stringsAsFactors = FALSE
  ),
  data.frame(
    workflow = "DStressR with EVC",
    promoter = comparison$promoter,
    compound = comparison$compound,
    effect = comparison$evc_huber_total_effect,
    stringsAsFactors = FALSE
  )
)
total_long$promoter <- factor(total_long$promoter, levels = rev(reporter_order))
total_long$compound <- factor(total_long$compound, levels = compound_order)
total_long$workflow <- factor(
  total_long$workflow,
  levels = c("DStressR without EVC", "DStressR with EVC"),
  labels = c("a", "b")
)
total_limit <- stats::quantile(abs(total_long$effect), 0.98, na.rm = TRUE)

make_workflow_diag <- function(d) {
  out <- do.call(rbind, lapply(split(d, d$workflow), function(wd) {
    tab <- perturbation_diagnostics(
      wd,
      mean_effect = "mean_response",
      variance_effect = "mean_response",
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

total_diag <- make_workflow_diag(diagnostics)
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
highlighted_compounds <- unique(label_data$perturbation)
compound_axis_labels <- function(x) {
  label_text <- vapply(x, function(label) {
    escaped <- gsub("\\\\", "\\\\\\\\", label)
    escaped <- gsub("'", "\\\\'", escaped)
    if (label %in% highlighted_compounds) {
      paste0("bold('", escaped, "')")
    } else {
      paste0("'", escaped, "'")
    }
  }, character(1))
  parse(text = label_text)
}

total_plot <- ggplot2$ggplot(total_long, ggplot2$aes(compound, promoter, fill = effect)) +
  ggplot2$geom_tile(width = 0.98, height = 0.95, color = "white", linewidth = 0.07) +
  ggplot2$facet_wrap(ggplot2$vars(workflow), ncol = 1, strip.position = "top") +
  ggplot2$scale_x_discrete(labels = compound_axis_labels) +
  ggplot2$scale_fill_gradient2(
    low = "#3f6fb5",
    mid = "white",
    high = "#c93f3f",
    midpoint = 0,
    limits = c(-total_limit, total_limit),
    oob = scales$squish,
    name = "Total effect"
  ) +
  ggplot2$labs(x = NULL, y = NULL) +
  ggplot2$theme_minimal(base_size = 11) +
  ggplot2$theme(
    panel.grid = ggplot2$element_blank(),
    strip.background = ggplot2$element_blank(),
    strip.clip = "off",
    strip.text = ggplot2$element_text(hjust = -0.055, face = "bold", size = 14, margin = ggplot2$margin(b = 1)),
    axis.text.x = ggplot2$element_text(angle = 63, hjust = 1, vjust = 1, size = 5.8, margin = ggplot2$margin(t = 1)),
    axis.text.y = ggplot2$element_text(size = 10.5),
    legend.position = "right",
    legend.title = ggplot2$element_text(size = 9),
    legend.text = ggplot2$element_text(size = 8),
    plot.margin = ggplot2$margin(0, 8, 0, 6)
  ) +
  ggplot2$guides(
    fill = ggplot2$guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barheight = grid::unit(28, "mm"),
      barwidth = grid::unit(4, "mm")
    )
  )

plot_rank_variance <- function(workflow, show_y = TRUE, show_legend = TRUE) {
  plot_data <- total_diag[total_diag$workflow == workflow, , drop = FALSE]
  labels <- label_data[label_data$workflow == workflow, , drop = FALSE]
  p <- ggplot2$ggplot(plot_data, ggplot2$aes(rank_abs_mean_effect, effect_variance)) +
    ggplot2$geom_point(ggplot2$aes(color = color_effect), size = 1.55, alpha = 0.82) +
    ggplot2$geom_line(
      data = plot_data[is.finite(plot_data$variance_trend), , drop = FALSE],
      ggplot2$aes(y = variance_trend),
      color = "#008A5B",
      linewidth = 0.65
    ) +
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
      axis.title = ggplot2$element_text(color = "#111827", size = 9),
      axis.text = ggplot2$element_text(color = "#4B5563", size = 8),
      legend.position = if (show_legend) "right" else "none",
      legend.justification = "center",
      legend.background = ggplot2$element_blank(),
      legend.key.height = grid::unit(12, "mm"),
      legend.key.width = grid::unit(3, "mm"),
      legend.title = ggplot2$element_text(size = 8),
      legend.text = ggplot2$element_text(size = 7),
      plot.margin = ggplot2$margin(0, 8, 4, 6)
    ) +
    ggplot2$guides(
      color = ggplot2$guide_colorbar(
        title.position = "top",
        title.hjust = 0.5,
        barwidth = grid::unit(3, "mm"),
        barheight = grid::unit(16, "mm")
      )
    ) +
    ggplot2$labs(
      x = expression("Rank of " * abs(hat(Delta * bar(y))[j]^"tot")),
      y = expression("Variance across promoters of " * hat(Delta * y)["\u00b7" * j]^"tot")
    )
  if (!show_y) {
    p <- p + ggplot2$labs(y = NULL)
  }
  if (nrow(labels) > 0 && requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p +
      ggrepel::geom_text_repel(
        data = labels,
        ggplot2$aes(label = perturbation_label),
        size = 3.05,
        fontface = "bold",
        min.segment.length = 0,
        box.padding = 0.2,
        point.padding = 0.1,
        max.overlaps = Inf,
        seed = 1,
        show.legend = FALSE
      )
  }
  p
}

combined <- gridExtra$arrangeGrob(
  total_plot,
  gridExtra$arrangeGrob(panel_label("c"), panel_label("d"), ncol = 2),
  gridExtra$arrangeGrob(
    plot_rank_variance("DStressR without EVC", show_y = TRUE, show_legend = FALSE),
    plot_rank_variance("DStressR with EVC", show_y = FALSE, show_legend = TRUE),
    ncol = 2
  ),
  ncol = 1,
  heights = c(2.2, 0.045, 1.0)
)

ggplot2$ggsave(
  file.path(effect_dir, "binsfeld_dstressr_effect_matrix_and_variance_diagnostic.png"),
  combined,
  width = 16.5,
  height = 11.2,
  dpi = 260,
  bg = "white",
  limitsize = FALSE
)
ggplot2$ggsave(
  file.path(effect_dir, "binsfeld_dstressr_effect_matrix_and_variance_diagnostic.pdf"),
  combined,
  width = 16.5,
  height = 11.2,
  limitsize = FALSE
)

utils::write.table(
  compound_score,
  file.path(effect_dir, "binsfeld_dstressr_effect_matrix_and_variance_diagnostic_compound_order.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Wrote combined DStressR effect-matrix diagnostic to: ", effect_dir)
