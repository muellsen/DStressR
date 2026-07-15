#!/usr/bin/env Rscript

source(file.path("analysis", "_helpers.R"))

suppressPackageStartupMessages({
  library(ggplot2)
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

print(summary_table)
print(cor_table)
print(region_counts)
message("Wrote pair-level comparison outputs to: ", out_dir)
