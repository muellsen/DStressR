source(file.path("analysis", "_helpers.R"))
load_destress_package()

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package `ggplot2` is required for effect matrix plots.", call. = FALSE)
}
ggplot2 <- asNamespace("ggplot2")
scales <- asNamespace("scales")

out_dir <- analysis_output_dir("binsfeld_effect_matrix_views")
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

comparison_file <- file.path(
  analysis_project_root(),
  "analysis", "outputs", "binsfeld_evc_calibrated",
  "evc_huber_comparison_to_binsfeld_and_default.tsv"
)
diagnostics_file <- file.path(
  analysis_project_root(),
  "analysis", "outputs", "binsfeld_variance_diagnostics",
  "ecoli_raw_vs_moderated_variance.tsv"
)
if (file.exists(comparison_file) && file.exists(diagnostics_file)) {
  comparison <- read_tsv_base(comparison_file)
  diagnostics <- read_tsv_base(diagnostics_file)
  required_comparison <- c(
    "promoter", "compound", "mean_z", "binsfeld_hit",
    "modeled_hit", "evc_huber_hit"
  )
  required_diagnostics <- c("workflow", "promoter", "compound", "mean_response")
  missing_comparison <- setdiff(required_comparison, names(comparison))
  if (length(missing_comparison) > 0) {
    stop("Missing required comparison columns: ", paste(missing_comparison, collapse = ", "), call. = FALSE)
  }
  missing_diagnostics <- setdiff(required_diagnostics, names(diagnostics))
  if (length(missing_diagnostics) > 0) {
    stop("Missing required diagnostics columns: ", paste(missing_diagnostics, collapse = ", "), call. = FALSE)
  }

  reporter_order <- c("acrABp", "marRABp", "micFp", "ompFp", "robp", "soxSp", "tolCp")
  diagnostics$workflow <- sub("without EV control", "without EVC", diagnostics$workflow, fixed = TRUE)
  diagnostics$workflow <- sub("with EV control", "with EVC", diagnostics$workflow, fixed = TRUE)
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

  comparison$hit_count <- as.integer(comparison$modeled_hit) +
    as.integer(comparison$evc_huber_hit)
  comparison$abs_sum <- abs(comparison$modeled_total_effect) +
    abs(comparison$evc_huber_total_effect)
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

  effect_long <- rbind(
    data.frame(
      method = "DStressR without EVC",
      promoter = comparison$promoter,
      compound = comparison$compound,
      effect = comparison$modeled_total_effect,
      hit = comparison$modeled_hit,
      stringsAsFactors = FALSE
    ),
    data.frame(
      method = "DStressR with EVC",
      promoter = comparison$promoter,
      compound = comparison$compound,
      effect = comparison$evc_huber_total_effect,
      hit = comparison$evc_huber_hit,
      stringsAsFactors = FALSE
    )
  )
  effect_long$promoter <- factor(effect_long$promoter, levels = rev(reporter_order))
  effect_long$compound <- factor(effect_long$compound, levels = compound_score$compound)
  effect_long$method <- factor(
    effect_long$method,
    levels = c("DStressR without EVC", "DStressR with EVC")
  )
  effect_limit <- stats::quantile(abs(effect_long$effect), 0.98, na.rm = TRUE)

  three_approach_plot <- ggplot2$ggplot(
    effect_long,
    ggplot2$aes(compound, promoter, fill = effect)
  ) +
    ggplot2$geom_tile(width = 0.98, height = 0.95, color = "white", linewidth = 0.07) +
    ggplot2$facet_grid(ggplot2$vars(method), switch = "y") +
    ggplot2$scale_fill_gradient2(
      low = "#3f6fb5",
      mid = "white",
      high = "#c93f3f",
      midpoint = 0,
      limits = c(-effect_limit, effect_limit),
      oob = scales$squish,
      name = "Total effect"
    ) +
    ggplot2$labs(x = NULL, y = NULL) +
    ggplot2$theme_minimal(base_size = 11) +
    ggplot2$theme(
      panel.grid = ggplot2$element_blank(),
      strip.text.y.left = ggplot2$element_text(
        angle = 0,
        face = "bold",
        size = 10.5,
        margin = ggplot2$margin(r = 6)
      ),
      strip.placement = "outside",
      axis.text.x = ggplot2$element_text(
        angle = 63,
        hjust = 1,
        vjust = 1,
        size = 5.1,
        margin = ggplot2$margin(t = 1)
      ),
      axis.text.y = ggplot2$element_text(size = 10.5),
      legend.position = "right",
      legend.title = ggplot2$element_text(size = 9),
      legend.text = ggplot2$element_text(size = 8),
      plot.margin = ggplot2$margin(6, 8, 5, 6)
    ) +
    ggplot2$guides(
      fill = ggplot2$guide_colorbar(
        title.position = "top",
        title.hjust = 0.5,
        barheight = grid::unit(36, "mm"),
        barwidth = grid::unit(4, "mm")
      )
    )

  ggplot2$ggsave(
    file.path(out_dir, "binsfeld_dstressr_total_effect_matrix.png"),
    three_approach_plot,
    width = 16.5,
    height = 6.2,
    dpi = 260,
    bg = "white",
    limitsize = FALSE
  )
  ggplot2$ggsave(
    file.path(out_dir, "binsfeld_dstressr_total_effect_matrix.pdf"),
    three_approach_plot,
    width = 16.5,
    height = 6.2,
    limitsize = FALSE
  )

  specific_long <- rbind(
    data.frame(
      method = "Original workflow",
      promoter = comparison$promoter,
      compound = comparison$compound,
      effect = comparison$mean_z,
      stringsAsFactors = FALSE
    ),
    data.frame(
      method = "DStressR without EVC",
      promoter = comparison$promoter,
      compound = comparison$compound,
      effect = comparison$modeled_effect,
      stringsAsFactors = FALSE
    ),
    data.frame(
      method = "DStressR with EVC",
      promoter = comparison$promoter,
      compound = comparison$compound,
      effect = comparison$evc_huber_effect,
      stringsAsFactors = FALSE
    )
  )
  specific_long$promoter <- factor(specific_long$promoter, levels = rev(reporter_order))
  specific_long$compound <- factor(specific_long$compound, levels = compound_score$compound)
  specific_long$method <- factor(
    specific_long$method,
    levels = c("Original workflow", "DStressR without EVC", "DStressR with EVC")
  )
  specific_limit <- stats::quantile(abs(specific_long$effect), 0.98, na.rm = TRUE)

  specific_plot <- ggplot2$ggplot(
    specific_long,
    ggplot2$aes(compound, promoter, fill = effect)
  ) +
    ggplot2$geom_tile(width = 0.98, height = 0.95, color = "white", linewidth = 0.07) +
    ggplot2$facet_grid(ggplot2$vars(method), switch = "y") +
    ggplot2$scale_fill_gradient2(
      low = "#3f6fb5",
      mid = "white",
      high = "#c93f3f",
      midpoint = 0,
      limits = c(-specific_limit, specific_limit),
      oob = scales$squish,
      name = "Matrix value"
    ) +
    ggplot2$labs(x = NULL, y = NULL) +
    ggplot2$theme_minimal(base_size = 11) +
    ggplot2$theme(
      panel.grid = ggplot2$element_blank(),
      strip.text.y.left = ggplot2$element_text(
        angle = 0,
        face = "bold",
        size = 10.5,
        margin = ggplot2$margin(r = 6)
      ),
      strip.placement = "outside",
      axis.text.x = ggplot2$element_text(
        angle = 63,
        hjust = 1,
        vjust = 1,
        size = 5.1,
        margin = ggplot2$margin(t = 1)
      ),
      axis.text.y = ggplot2$element_text(size = 10.5),
      legend.position = "right",
      legend.title = ggplot2$element_text(size = 9),
      legend.text = ggplot2$element_text(size = 8),
      plot.margin = ggplot2$margin(6, 8, 5, 6)
    ) +
    ggplot2$guides(
      fill = ggplot2$guide_colorbar(
        title.position = "top",
        title.hjust = 0.5,
        barheight = grid::unit(36, "mm"),
        barwidth = grid::unit(4, "mm")
      )
    )

  ggplot2$ggsave(
    file.path(out_dir, "binsfeld_three_specific_effect_matrices.png"),
    specific_plot,
    width = 16.5,
    height = 8.5,
    dpi = 260,
    bg = "white",
    limitsize = FALSE
  )
  ggplot2$ggsave(
    file.path(out_dir, "binsfeld_three_specific_effect_matrices.pdf"),
    specific_plot,
    width = 16.5,
    height = 8.5,
    limitsize = FALSE
  )
  utils::write.table(
    compound_score,
    file.path(out_dir, "binsfeld_dstressr_total_effect_matrix_compound_order.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

message("Wrote Binsfeld effect matrix views to: ", out_dir)
