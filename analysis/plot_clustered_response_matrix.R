#!/usr/bin/env Rscript

source(file.path("analysis", "_helpers.R"))

suppressPackageStartupMessages({
  library(ggplot2)
})

if (requireNamespace("DStressR", quietly = TRUE)) {
  suppressPackageStartupMessages(library(DStressR))
} else {
  load_destress_package()
}

out_dir <- file.path(getwd(), "analysis", "outputs", "normalized_matrix", "clusters")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

matrix_file <- file.path(
  getwd(),
  "analysis",
  "outputs",
  "normalized_matrix",
  "normalized_promoter_by_compound_matrix_destress_eb_effect.tsv"
)

mat_df <- read.delim(matrix_file, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
mat <- as.matrix(mat_df[, setdiff(names(mat_df), "promoter"), drop = FALSE])
mode(mat) <- "numeric"
rownames(mat) <- mat_df$promoter

long_df <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
names(long_df) <- c("promoter", "compound", "specific_effect")
long_df <- long_df[is.finite(long_df$specific_effect), , drop = FALSE]

promoter_k <- 6
compound_k <- 14
plot_response_clustered_heatmap(
  long_df,
  value = "specific_effect",
  promoter = "promoter",
  compound = "compound",
  compound_label = "compound",
  show_compound_ids = FALSE,
  top_n_compounds = 400,
  n_promoter_clusters = promoter_k,
  n_compound_clusters = compound_k,
  file = file.path(out_dir, "pretty_clustered_heatmap_top400.png"),
  width = 14,
  height = 8,
  title = "Hierarchically clustered promoter-by-compound response heatmap",
  subtitle = "Top 400 compounds by mean absolute DStressR EB effect",
  show_rownames = TRUE,
  show_colnames = FALSE
)
plot_response_clustered_heatmap(
  long_df,
  value = "specific_effect",
  promoter = "promoter",
  compound = "compound",
  compound_label = "compound",
  show_compound_ids = FALSE,
  top_n_compounds = 400,
  n_promoter_clusters = promoter_k,
  n_compound_clusters = compound_k,
  file = file.path(out_dir, "pretty_clustered_heatmap_top400.pdf"),
  width = 14,
  height = 8,
  title = "Hierarchically clustered promoter-by-compound response heatmap",
  subtitle = "Top 400 compounds by mean absolute DStressR EB effect",
  show_rownames = TRUE,
  show_colnames = FALSE
)
plot_response_clustered_heatmap(
  long_df,
  value = "specific_effect",
  promoter = "promoter",
  compound = "compound",
  compound_label = "compound",
  show_compound_ids = FALSE,
  top_n_compounds = Inf,
  n_promoter_clusters = promoter_k,
  n_compound_clusters = compound_k,
  file = file.path(out_dir, "pretty_clustered_heatmap_full.png"),
  width = 16,
  height = 8,
  title = "Full hierarchically clustered response heatmap",
  subtitle = "All compounds; column labels suppressed",
  show_rownames = TRUE,
  show_colnames = FALSE
)

p_block <- plot_response_cluster_blocks(
  long_df,
  value = "specific_effect",
  promoter = "promoter",
  compound = "compound",
  compound_label = "compound",
  show_compound_ids = FALSE,
  n_promoter_clusters = promoter_k,
  n_compound_clusters = compound_k,
  title = "Clustered promoter-by-compound response map",
  subtitle = paste0(
    promoter_k, " promoter clusters x ", compound_k,
    " compound clusters; tile labels show compounds per cluster"
  )
)

promoter_assignments <- attr(p_block, "promoter_clusters")
compound_assignments <- attr(p_block, "compound_clusters")
block_df <- attr(p_block, "block_summary")
names(block_df) <- sub("^\\.", "", names(block_df))

write.table(
  promoter_assignments,
  file.path(out_dir, "promoter_cluster_assignments.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
write.table(
  compound_assignments,
  file.path(out_dir, "compound_cluster_assignments.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
write.table(
  block_df,
  file.path(out_dir, "promoter_compound_cluster_block_summary.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

top_compounds <- do.call(
  rbind,
  lapply(split(compound_assignments, compound_assignments$compound_cluster), function(d) {
    d <- d[order(-d$mean_abs_effect), ]
    utils::head(d, 12)
  })
)
write.table(
  top_compounds,
  file.path(out_dir, "top_compounds_per_cluster.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

ggsave(
  file.path(out_dir, "promoter_compound_cluster_block_heatmap.png"),
  p_block,
  width = 9,
  height = 5.5,
  dpi = 300,
  bg = "white"
)
ggsave(
  file.path(out_dir, "promoter_compound_cluster_block_heatmap.pdf"),
  p_block,
  width = 9,
  height = 5.5,
  bg = "white"
)

row_hc <- attr(p_block, "row_hclust")
col_hc <- attr(p_block, "col_hclust")
ordered_mat <- mat[row_hc$order, col_hc$order, drop = FALSE]
full_df <- as.data.frame(as.table(ordered_mat), stringsAsFactors = FALSE)
names(full_df) <- c("promoter", "compound_display", "effect")
full_df$promoter <- factor(full_df$promoter, levels = rownames(ordered_mat))
full_df$compound_display <- factor(full_df$compound_display, levels = colnames(ordered_mat))

limit <- stats::quantile(abs(full_df$effect), 0.98, na.rm = TRUE)
if (!is.finite(limit) || limit <= 0) {
  limit <- max(abs(full_df$effect), na.rm = TRUE)
}
full_df$plot_effect <- pmax(pmin(full_df$effect, limit), -limit)

row_cluster <- stats::setNames(
  promoter_assignments$promoter_cluster,
  promoter_assignments$promoter
)
col_cluster <- stats::setNames(
  compound_assignments$compound_cluster,
  compound_assignments$compound_display
)
ordered_row_clusters <- row_cluster[rownames(ordered_mat)]
ordered_col_clusters <- col_cluster[colnames(ordered_mat)]
row_boundaries <- which(ordered_row_clusters[-1] != ordered_row_clusters[-length(ordered_row_clusters)]) + 0.5
col_boundaries <- which(ordered_col_clusters[-1] != ordered_col_clusters[-length(ordered_col_clusters)]) + 0.5

p_ordered <- ggplot(full_df, aes(compound_display, promoter, fill = plot_effect)) +
  geom_raster() +
  geom_hline(yintercept = row_boundaries, color = "#334155", linewidth = 0.25) +
  geom_vline(xintercept = col_boundaries, color = "#334155", linewidth = 0.18) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-limit, limit),
    name = "DStressR\nEB effect"
  ) +
  theme_light(base_size = 8) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "bottom"
  ) +
  labs(
    title = "Hierarchically clustered full response matrix",
    subtitle = "Lines mark promoter and compound cluster boundaries",
    x = "Compounds",
    y = "Promoters"
  )

ggsave(
  file.path(out_dir, "clustered_full_response_matrix_with_boundaries.png"),
  p_ordered,
  width = 14,
  height = 6.5,
  dpi = 300,
  bg = "white"
)
ggsave(
  file.path(out_dir, "clustered_full_response_matrix_with_boundaries.pdf"),
  p_ordered,
  width = 14,
  height = 6.5,
  bg = "white"
)

message("Wrote clustered response representations to: ", out_dir)
message("Promoter clusters: ", promoter_k)
message("Compound clusters: ", compound_k)
