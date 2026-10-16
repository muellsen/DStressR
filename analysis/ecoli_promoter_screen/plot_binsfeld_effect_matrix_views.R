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
comparison <- artifact$annotations$comparison
diagnostics <- artifact$diagnostics$display
if (nrow(comparison)) {
  required_comparison <- c(
    "promoter", "compound", "mean_z", "binsfeld_hit",
    "estimated_exponent_hit", "evc_huber_hit"
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
  names(total_no_evc)[3] <- "estimated_exponent_total_effect"
  total_evc <- diagnostics[
    diagnostics$workflow == "DStressR with EVC",
    c("promoter", "compound", "mean_response"),
    drop = FALSE
  ]
  names(total_evc)[3] <- "evc_huber_total_effect"
  comparison <- merge(comparison, total_no_evc, by = c("promoter", "compound"), all.x = TRUE, sort = FALSE)
  comparison <- merge(comparison, total_evc, by = c("promoter", "compound"), all.x = TRUE, sort = FALSE)

  comparison$hit_count <- as.integer(comparison$estimated_exponent_hit) +
    as.integer(comparison$evc_huber_hit)
  comparison$abs_sum <- abs(comparison$estimated_exponent_total_effect) +
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
      effect = comparison$estimated_exponent_total_effect,
      hit = comparison$estimated_exponent_hit,
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
      effect = comparison$estimated_exponent_effect,
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
