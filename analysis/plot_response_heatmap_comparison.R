#!/usr/bin/env Rscript

source(file.path("analysis", "_helpers.R"))

suppressPackageStartupMessages({
  library(ggplot2)
})

out_dir <- comparison_results_dir("response_heatmaps")
top_n <- as.integer(Sys.getenv("DSTRESSR_RESPONSE_HEATMAP_TOP_N", unset = "220"))

read_response <- function(method) {
  path <- package_results_dir(paste0(method, "_response.tsv"))
  if (!file.exists(path)) {
    stop(
      "Missing package response output: ", path,
      "\nRun the package export that writes legacy_ratio_response.tsv and modeled_response.tsv first.",
      call. = FALSE
    )
  }
  tab <- read_tsv_base(path)
  required <- c("promoter", "compound", "compound_label", "response_centered")
  missing <- setdiff(required, names(tab))
  if (length(missing) > 0) {
    stop("Response output is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  tab
}

make_response_matrix <- function(tab, promoters, compounds) {
  d <- stats::aggregate(
    response_centered ~ promoter + compound,
    tab,
    mean,
    na.rm = TRUE
  )
  mat <- matrix(
    NA_real_,
    nrow = length(promoters),
    ncol = length(compounds),
    dimnames = list(promoters, compounds)
  )
  keep <- d$promoter %in% promoters & d$compound %in% compounds
  d <- d[keep, , drop = FALSE]
  mat[cbind(d$promoter, d$compound)] <- d$response_centered
  mat
}

cluster_order <- function(mat, margin) {
  x <- mat
  x[!is.finite(x)] <- 0
  if (margin == 1) {
    if (nrow(x) < 2) {
      return(seq_len(nrow(x)))
    }
    return(stats::hclust(stats::dist(x))$order)
  }
  if (ncol(x) < 2) {
    return(seq_len(ncol(x)))
  }
  stats::hclust(stats::dist(t(x)))$order
}

matrix_to_long <- function(mat, method) {
  out <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
  names(out) <- c("promoter", "compound", "response_centered")
  out$method <- method
  out
}

plot_heatmap <- function(plot_data, title, subtitle, stem, width = 13.5) {
  limit <- stats::quantile(abs(plot_data$response_centered), 0.98, na.rm = TRUE)
  if (!is.finite(limit) || limit <= 0) {
    limit <- max(abs(plot_data$response_centered), na.rm = TRUE)
  }

  p <- ggplot(plot_data, aes(compound, promoter, fill = response_centered)) +
    geom_raster() +
    scale_fill_gradient2(
      low = "#2563eb",
      mid = "#f8fafc",
      high = "#b91c1c",
      midpoint = 0,
      limits = c(-limit, limit),
      na.value = "#e5e7eb",
      name = "Centered response"
    ) +
    theme_light(base_size = 8) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid = element_blank(),
      legend.position = "bottom",
      plot.title = element_text(face = "bold"),
      plot.title.position = "plot",
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      x = paste0("Shared top ", top_n, " compounds by mean absolute centered response"),
      y = "Promoter"
    )

  ggsave(file.path(out_dir, paste0(stem, ".png")), p, width = width, height = 7.8, dpi = 300, bg = "white")
  ggsave(file.path(out_dir, paste0(stem, ".pdf")), p, width = width, height = 7.8, bg = "white")
  p
}

write_matrix <- function(mat, stem) {
  write.table(
    data.frame(promoter = rownames(mat), mat, check.names = FALSE),
    file.path(out_dir, paste0(stem, "_matrix.tsv")),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
}

ratio <- read_response("legacy_ratio")
modeled <- read_response("modeled")

combined_for_selection <- rbind(
  ratio[, c("compound", "compound_label", "response_centered"), drop = FALSE],
  modeled[, c("compound", "compound_label", "response_centered"), drop = FALSE]
)
compound_summary <- stats::aggregate(
  abs(response_centered) ~ compound + compound_label,
  combined_for_selection,
  mean,
  na.rm = TRUE
)
names(compound_summary)[3] <- "mean_abs_centered_response"
compound_summary <- compound_summary[order(-compound_summary$mean_abs_centered_response), , drop = FALSE]
selected_compounds <- utils::head(compound_summary$compound, top_n)

promoters <- sort(intersect(unique(ratio$promoter), unique(modeled$promoter)))
ratio_mat <- make_response_matrix(ratio, promoters, selected_compounds)
modeled_mat <- make_response_matrix(modeled, promoters, selected_compounds)

row_order <- rownames(modeled_mat)[cluster_order(modeled_mat, 1)]
col_order <- colnames(modeled_mat)[cluster_order(modeled_mat, 2)]
ratio_mat <- ratio_mat[row_order, col_order, drop = FALSE]
modeled_mat <- modeled_mat[row_order, col_order, drop = FALSE]
diff_mat <- modeled_mat - ratio_mat

write_matrix(ratio_mat, "legacy_ratio_response_heatmap")
write_matrix(modeled_mat, "modeled_response_heatmap")
write_matrix(diff_mat, "modeled_minus_ratio_response_heatmap")
write.table(
  data.frame(promoter = row_order, promoter_order = seq_along(row_order)),
  file.path(out_dir, "shared_promoter_order.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
write.table(
  data.frame(compound = col_order, compound_order = seq_along(col_order)),
  file.path(out_dir, "shared_compound_order.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

ratio_long <- matrix_to_long(ratio_mat, "Standard LUX/growth ratio")
modeled_long <- matrix_to_long(modeled_mat, "Modeled growth response")
diff_long <- matrix_to_long(diff_mat, "Modeled - baseline alpha=1")
plot_data <- rbind(ratio_long, modeled_long)
plot_data$promoter <- factor(plot_data$promoter, levels = rev(row_order))
plot_data$compound <- factor(plot_data$compound, levels = col_order)
plot_data$method <- factor(plot_data$method, levels = c(
  "Standard LUX/growth ratio",
  "Modeled growth response"
))
diff_long$promoter <- factor(diff_long$promoter, levels = rev(row_order))
diff_long$compound <- factor(diff_long$compound, levels = col_order)

ratio_plot <- plot_heatmap(
  subset(plot_data, method == "Standard LUX/growth ratio"),
  title = "Screen response heatmap: standard LUX/growth ratio",
  subtitle = "DMSO-centered log2(LUX / growth), matching the legacy response definition",
  stem = "legacy_ratio_response_heatmap"
)
modeled_plot <- plot_heatmap(
  subset(plot_data, method == "Modeled growth response"),
  title = "Screen response heatmap: modeled growth response",
  subtitle = "DMSO-centered log2(LUX) - alpha[g] log2(growth), using promoter-specific shrunken alpha[g]",
  stem = "modeled_response_heatmap"
)
diff_plot <- plot_heatmap(
  diff_long,
  title = "Response difference heatmap: modeled minus baseline alpha = 1",
  subtitle = "DMSO-centered modeled response minus DMSO-centered log2(LUX / growth)",
  stem = "modeled_minus_ratio_response_heatmap"
)

combined_limit <- stats::quantile(abs(plot_data$response_centered), 0.98, na.rm = TRUE)
combined_plot <- ggplot(plot_data, aes(compound, promoter, fill = response_centered)) +
  geom_raster() +
  facet_grid(method ~ .) +
  scale_fill_gradient2(
    low = "#2563eb",
    mid = "#f8fafc",
    high = "#b91c1c",
    midpoint = 0,
    limits = c(-combined_limit, combined_limit),
    na.value = "#e5e7eb",
    name = "Centered response"
  ) +
  theme_light(base_size = 8) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "bottom",
    strip.background = element_rect(fill = "#f8fafc", color = "#cbd5e1"),
    strip.text = element_text(face = "bold", color = "#0f172a"),
    plot.title = element_text(face = "bold"),
    plot.title.position = "plot",
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  labs(
    title = "Screen response heatmaps with matched row and column ordering",
    subtitle = "Rows and columns ordered once from the modeled response; both panels use the same compound set and color scale",
    x = paste0("Shared top ", top_n, " compounds by mean absolute centered response"),
    y = "Promoter"
  )
ggsave(file.path(out_dir, "matched_response_heatmaps.png"), combined_plot, width = 13.5, height = 11.5, dpi = 300, bg = "white")
ggsave(file.path(out_dir, "matched_response_heatmaps.pdf"), combined_plot, width = 13.5, height = 11.5, bg = "white")

diff_limit <- stats::quantile(abs(diff_long$response_centered), 0.98, na.rm = TRUE)
if (!is.finite(diff_limit) || diff_limit <= 0) {
  diff_limit <- max(abs(diff_long$response_centered), na.rm = TRUE)
}
diff_panel <- ggplot(diff_long, aes(compound, promoter, fill = response_centered)) +
  geom_raster() +
  scale_fill_gradient2(
    low = "#2563eb",
    mid = "#f8fafc",
    high = "#b91c1c",
    midpoint = 0,
    limits = c(-diff_limit, diff_limit),
    na.value = "#e5e7eb",
    name = "Modeled - baseline"
  ) +
  theme_light(base_size = 8) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(face = "bold"),
    plot.title.position = "plot",
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  labs(
    title = "Response difference heatmap: modeled growth response minus baseline alpha = 1",
    subtitle = "Same row and column ordering as the matched response heatmaps",
    x = paste0("Shared top ", top_n, " compounds by mean absolute centered response"),
    y = "Promoter"
  )
ggsave(file.path(out_dir, "modeled_minus_ratio_response_heatmap.png"), diff_panel, width = 13.5, height = 7.8, dpi = 300, bg = "white")
ggsave(file.path(out_dir, "modeled_minus_ratio_response_heatmap.pdf"), diff_panel, width = 13.5, height = 7.8, bg = "white")

summary <- data.frame(
  method = c("legacy_ratio", "modeled_response"),
  rows = c(nrow(ratio), nrow(modeled)),
  promoters = c(length(unique(ratio$promoter)), length(unique(modeled$promoter))),
  compounds = c(length(unique(ratio$compound)), length(unique(modeled$compound))),
  top_n_compounds = top_n,
  stringsAsFactors = FALSE
)
write.table(
  summary,
  file.path(out_dir, "response_heatmap_summary.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

message("Wrote response heatmaps to: ", out_dir)
print(summary)
