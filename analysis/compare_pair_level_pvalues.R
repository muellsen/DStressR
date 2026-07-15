#!/usr/bin/env Rscript

source(file.path("analysis", "_helpers.R"))

suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)
})

methods <- c("median_polish", "destress_standard", "destress_moderated")
out_dir <- comparison_results_dir("pair_level")
adjustment <- comparison_adjustment()

pair_table <- merge_package_pair_results(methods)
pair_table <- add_hit_columns(pair_table, methods, fdr = 0.05, adjustment = adjustment)

write.table(
  pair_table,
  file.path(out_dir, "pair_level_pvalue_comparison.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

summary_table <- do.call(rbind, lapply(methods, function(method) {
  data.frame(
    method = method,
    label = method_label(method),
    common_pair_rows = nrow(pair_table),
    finite_pvalues = sum(is.finite(pair_table[[paste0(method, "_pvalue")]])),
    adjustment = adjustment,
    finite_adjusted_pvalues = sum(is.finite(pair_table[[padj_column(method, adjustment)]])),
    raw_p_lt_0.05 = sum(pair_table[[paste0(method, "_pvalue")]] < 0.05, na.rm = TRUE),
    adjusted_p_lt_0.05 = sum(pair_table[[method_hit_column(method)]], na.rm = TRUE),
    median_pvalue = stats::median(pair_table[[paste0(method, "_pvalue")]], na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))

method_pairs <- utils::combn(methods, 2, simplify = FALSE)
cor_table <- do.call(rbind, lapply(method_pairs, function(method_pair) {
  x <- paste0(method_pair[1], "_pvalue")
  y <- paste0(method_pair[2], "_pvalue")
  keep <- is.finite(pair_table[[x]]) & is.finite(pair_table[[y]])
  data.frame(
    method_x = method_pair[1],
    method_y = method_pair[2],
    label_x = method_label(method_pair[1]),
    label_y = method_label(method_pair[2]),
    n = sum(keep),
    spearman_neglog10p = stats::cor(
      safe_neglog10(pair_table[[x]][keep]),
      safe_neglog10(pair_table[[y]][keep]),
      method = "spearman"
    ),
    stringsAsFactors = FALSE
  )
}))

membership <- pair_table[, c("promoter", "compound", "pair_id"), drop = FALSE]
for (method in methods) {
  membership[[method]] <- pair_table[[method_hit_column(method)]]
}
membership$hit_region <- paste(
  ifelse(membership$median_polish, "median_polish", ""),
  ifelse(membership$destress_standard, "destress_standard", ""),
  ifelse(membership$destress_moderated, "destress_moderated", ""),
  sep = ";"
)
membership$hit_region <- gsub(";+", ";", membership$hit_region)
membership$hit_region <- gsub("^;+|;+$", "", membership$hit_region)
membership$hit_region[membership$hit_region == ""] <- "not_hit"

region_counts <- as.data.frame(table(membership$hit_region), stringsAsFactors = FALSE)
names(region_counts) <- c("hit_region", "n")
region_counts <- region_counts[order(region_counts$hit_region), , drop = FALSE]

write.table(summary_table, file.path(out_dir, "pair_level_method_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(cor_table, file.path(out_dir, "pair_level_method_correlations.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(membership, file.path(out_dir, "pair_level_hit_membership.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(region_counts, file.path(out_dir, "pair_level_hit_region_counts.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

scatter_df <- pair_table
scatter_df$median_polish_neglog10p <- safe_neglog10(scatter_df$median_polish_pvalue)
scatter_df$destress_standard_neglog10p <- safe_neglog10(scatter_df$destress_standard_pvalue)
scatter_df$destress_moderated_neglog10p <- safe_neglog10(scatter_df$destress_moderated_pvalue)

standard_vs_median <- scatter_df[
  is.finite(scatter_df$median_polish_pvalue) &
    is.finite(scatter_df$destress_standard_pvalue),
  ,
  drop = FALSE
]
standard_vs_median$both_global_hits <-
  standard_vs_median[[method_hit_column("median_polish")]] &
  standard_vs_median[[method_hit_column("destress_standard")]]
standard_vs_median$hit_class <- "Not significant in both"
standard_vs_median$hit_class[
  standard_vs_median[[method_hit_column("median_polish")]] &
    !standard_vs_median[[method_hit_column("destress_standard")]]
] <- "Median-polish max-p only"
standard_vs_median$hit_class[
  !standard_vs_median[[method_hit_column("median_polish")]] &
    standard_vs_median[[method_hit_column("destress_standard")]]
] <- "DStressR standard only"
standard_vs_median$hit_class[standard_vs_median$both_global_hits] <- "Both"
standard_vs_median$hit_class <- factor(
  standard_vs_median$hit_class,
  levels = c("Not significant in both", "Median-polish max-p only", "DStressR standard only", "Both")
)

pair_pvalue_scatter <- function(d, x_col, y_col, value_label, zoom = FALSE) {
  title <- if (zoom) {
    paste0(value_label, ", zoomed to [0, 0.1]")
  } else {
    paste0(value_label, ", full range")
  }
  ggplot(d, aes(.data[[x_col]], .data[[y_col]])) +
    geom_point(aes(color = hit_class), alpha = 0.5, size = 1.05, stroke = 0) +
    geom_abline(slope = 1, intercept = 0, color = "#0f172a", linewidth = 0.35) +
    scale_color_manual(
      values = c(
        "Not significant in both" = "#94a3b8",
        "Median-polish max-p only" = "#f97316",
        "DStressR standard only" = "#2563eb",
        "Both" = "#b91c1c"
      ),
      name = paste0(adjustment, " BH hit")
    ) +
    coord_equal(xlim = if (zoom) c(0, 0.1) else c(0, 1),
                ylim = if (zoom) c(0, 0.1) else c(0, 1),
                expand = FALSE) +
    theme_bw(base_size = 9) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      plot.title = element_text(face = "bold"),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    ) +
    labs(
      title = title,
      subtitle = "One point per promoter-compound pair",
      x = paste0("Median-polish max-p ", value_label),
      y = paste0("DStressR standard ", value_label)
    )
}

raw_full <- pair_pvalue_scatter(
  standard_vs_median,
  "median_polish_pvalue",
  "destress_standard_pvalue",
  "raw p-values",
  zoom = FALSE
)
raw_zoom <- pair_pvalue_scatter(
  standard_vs_median,
  "median_polish_pvalue",
  "destress_standard_pvalue",
  "raw p-values",
  zoom = TRUE
)
raw_combined <- arrangeGrob(
  raw_full + theme(legend.position = "none"),
  raw_zoom + theme(legend.position = "none"),
  ncol = 2,
  top = "Median-polish max-p model vs DStressR standard model raw p-values"
)

adjusted_x <- padj_column("median_polish", adjustment)
adjusted_y <- padj_column("destress_standard", adjustment)
adjusted_label <- paste0(adjustment, " BH adjusted p-values")
adjusted_df <- standard_vs_median[
  is.finite(standard_vs_median[[adjusted_x]]) &
    is.finite(standard_vs_median[[adjusted_y]]),
  ,
  drop = FALSE
]
adjusted_full <- pair_pvalue_scatter(
  adjusted_df,
  adjusted_x,
  adjusted_y,
  adjusted_label,
  zoom = FALSE
)
adjusted_zoom <- pair_pvalue_scatter(
  adjusted_df,
  adjusted_x,
  adjusted_y,
  adjusted_label,
  zoom = TRUE
)
adjusted_combined <- arrangeGrob(
  adjusted_full + theme(legend.position = "none"),
  adjusted_zoom + theme(legend.position = "none"),
  ncol = 2,
  top = paste0("Median-polish max-p model vs DStressR standard model ", adjusted_label)
)

p1 <- ggplot(scatter_df, aes(median_polish_neglog10p, destress_moderated_neglog10p)) +
  geom_point(alpha = 0.12, size = 0.35, color = "#1f2937") +
  geom_abline(slope = 1, intercept = 0, color = "#b91c1c", linewidth = 0.4) +
  coord_equal() +
  theme_bw(base_size = 10) +
  labs(
    title = "Pair-level package p-values",
    subtitle = "One row per promoter-compound pair",
    x = "Median-polish max-p model -log10(p)",
    y = "DStressR moderated model -log10(p)"
  )

p2 <- ggplot(scatter_df, aes(destress_standard_neglog10p, destress_moderated_neglog10p)) +
  geom_point(alpha = 0.12, size = 0.35, color = "#1f2937") +
  geom_abline(slope = 1, intercept = 0, color = "#b91c1c", linewidth = 0.4) +
  coord_equal() +
  theme_bw(base_size = 10) +
  labs(
    title = "Standard vs moderated DStressR package p-values",
    subtitle = "One row per promoter-compound pair",
    x = "DStressR standard -log10(p)",
    y = "DStressR moderated model -log10(p)"
  )

ggsave(file.path(out_dir, "median_polish_vs_destress_moderated_pair_pvalues.png"),
       p1, width = 6.5, height = 6, dpi = 300, bg = "white")
ggsave(file.path(out_dir, "median_polish_vs_destress_moderated_pair_pvalues.pdf"),
       p1, width = 6.5, height = 6, bg = "white")
ggsave(file.path(out_dir, "destress_standard_vs_moderated_pair_pvalues.png"),
       p2, width = 6.5, height = 6, dpi = 300, bg = "white")
ggsave(file.path(out_dir, "destress_standard_vs_moderated_pair_pvalues.pdf"),
       p2, width = 6.5, height = 6, bg = "white")

ggsave(file.path(out_dir, "median_polish_vs_destress_standard_raw_pvalues_full.png"),
       raw_full, width = 6.5, height = 6.2, dpi = 300, bg = "white")
ggsave(file.path(out_dir, "median_polish_vs_destress_standard_raw_pvalues_full.pdf"),
       raw_full, width = 6.5, height = 6.2, bg = "white")
ggsave(file.path(out_dir, "median_polish_vs_destress_standard_raw_pvalues_zoom_0_0.1.png"),
       raw_zoom, width = 6.5, height = 6.2, dpi = 300, bg = "white")
ggsave(file.path(out_dir, "median_polish_vs_destress_standard_raw_pvalues_zoom_0_0.1.pdf"),
       raw_zoom, width = 6.5, height = 6.2, bg = "white")
ggsave(file.path(out_dir, "median_polish_vs_destress_standard_raw_pvalues_full_and_zoom.png"),
       raw_combined, width = 12, height = 5.8, dpi = 300, bg = "white")
ggsave(file.path(out_dir, "median_polish_vs_destress_standard_raw_pvalues_full_and_zoom.pdf"),
       raw_combined, width = 12, height = 5.8, bg = "white")

ggsave(file.path(out_dir, "median_polish_vs_destress_standard_adjusted_pvalues_full.png"),
       adjusted_full, width = 6.5, height = 6.2, dpi = 300, bg = "white")
ggsave(file.path(out_dir, "median_polish_vs_destress_standard_adjusted_pvalues_full.pdf"),
       adjusted_full, width = 6.5, height = 6.2, bg = "white")
ggsave(file.path(out_dir, "median_polish_vs_destress_standard_adjusted_pvalues_zoom_0_0.1.png"),
       adjusted_zoom, width = 6.5, height = 6.2, dpi = 300, bg = "white")
ggsave(file.path(out_dir, "median_polish_vs_destress_standard_adjusted_pvalues_zoom_0_0.1.pdf"),
       adjusted_zoom, width = 6.5, height = 6.2, bg = "white")
ggsave(file.path(out_dir, "median_polish_vs_destress_standard_adjusted_pvalues_full_and_zoom.png"),
       adjusted_combined, width = 12, height = 5.8, dpi = 300, bg = "white")
ggsave(file.path(out_dir, "median_polish_vs_destress_standard_adjusted_pvalues_full_and_zoom.pdf"),
       adjusted_combined, width = 12, height = 5.8, bg = "white")

print(summary_table)
print(cor_table)
print(region_counts)
message("Wrote pair-level comparison outputs to: ", out_dir)
