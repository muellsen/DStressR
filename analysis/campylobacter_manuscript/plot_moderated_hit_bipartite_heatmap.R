#!/usr/bin/env Rscript

source(file.path("analysis", "_helpers.R"))

suppressPackageStartupMessages({
  library(ggplot2)
})

out_dir <- comparison_results_dir("hit_bipartite_heatmap")
adjustment <- comparison_adjustment()
dstressr_method <- "destress_moderated"
methods <- c(dstressr_method, "median_polish")

libmap <- read_tsv_base(libmap_path())
libmap$libplate <- paste0("lp", libmap[["Library plate"]])
libmap$compound <- paste(libmap$libplate, libmap[["Well"]], sep = "_")
libmap$ProductName <- ifelse(
  is.na(libmap$ProductName) | libmap$ProductName == "NA" | libmap$ProductName == "",
  libmap[["Catalog Number"]],
  libmap$ProductName
)

moderated <- read_package_pair_results(dstressr_method)
median_polish <- read_package_pair_results("median_polish")

moderated_padj <- padj_column(dstressr_method, adjustment)
median_padj <- padj_column("median_polish", adjustment)

moderated$hit <- is.finite(moderated[[moderated_padj]]) & moderated[[moderated_padj]] < 0.05
median_polish$hit <- is.finite(median_polish[[median_padj]]) & median_polish[[median_padj]] < 0.05

moderated_hits <- moderated[moderated$hit, c("promoter", "compound"), drop = FALSE]
if (nrow(moderated_hits) == 0) {
  stop("No DStressR default moderated hits available for bipartite heatmap.", call. = FALSE)
}

compound_promoters <- split(moderated_hits$promoter, moderated_hits$compound)
pcmeA_unique_compounds <- names(compound_promoters)[vapply(
  compound_promoters,
  function(x) identical(sort(unique(x)), "PcmeA"),
  logical(1)
)]
moderated_hits <- moderated_hits[!moderated_hits$compound %in% pcmeA_unique_compounds, , drop = FALSE]
if (nrow(moderated_hits) == 0) {
  stop("No hits remain after removing PcmeA-unique compounds.", call. = FALSE)
}

promoters <- sort(unique(moderated$promoter))
compounds <- sort(unique(moderated_hits$compound))

moderated_matrix <- xtabs(hit ~ promoter + compound, moderated)
moderated_matrix <- moderated_matrix[promoters, compounds, drop = FALSE] > 0

cluster_order <- function(mat, margin) {
  if (margin == 1) {
    if (nrow(mat) < 2) {
      return(rownames(mat))
    }
    d <- stats::dist(as.matrix(mat), method = "binary")
    return(rownames(mat)[stats::hclust(d, method = "average")$order])
  }
  if (ncol(mat) < 2) {
    return(colnames(mat))
  }
  d <- stats::dist(t(as.matrix(mat)), method = "binary")
  colnames(mat)[stats::hclust(d, method = "average")$order]
}

row_order <- cluster_order(moderated_matrix, 1)
compound_degree <- colSums(moderated_matrix)
col_order <- names(sort(compound_degree, decreasing = TRUE))
compound_lookup <- unique(libmap[, c("compound", "ProductName", "Catalog Number", "Target"), drop = FALSE])

moderated_matrix <- moderated_matrix[row_order, col_order, drop = FALSE]
compound_labels <- data.frame(compound = col_order, stringsAsFactors = FALSE)
compound_labels <- merge(compound_labels, compound_lookup, by = "compound", all.x = TRUE, sort = FALSE)
compound_labels <- compound_labels[match(col_order, compound_labels$compound), , drop = FALSE]
compound_labels$label <- ifelse(
  is.na(compound_labels$ProductName) | compound_labels$ProductName == "" | compound_labels$ProductName == "NA",
  compound_labels$compound,
  compound_labels$ProductName
)
compound_labels$label <- iconv(compound_labels$label, from = "", to = "ASCII//TRANSLIT")
compound_labels$label[is.na(compound_labels$label)] <- compound_labels$compound[is.na(compound_labels$label)]
compound_labels$label <- gsub("\u03b3", "gamma", compound_labels$label, fixed = TRUE)
too_long <- nchar(compound_labels$label) > 36
compound_labels$label[too_long] <- paste0(substr(compound_labels$label[too_long], 1, 33), "...")
compound_axis_labels <- stats::setNames(compound_labels$label, compound_labels$compound)

make_matrix <- function(tab, method) {
  m <- xtabs(hit ~ promoter + compound, tab)
  m <- m[row_order, col_order, drop = FALSE] > 0
  effects <- xtabs(tab[[paste0(method, "_effect")]] ~ promoter + compound, tab)
  effects <- effects[row_order, col_order, drop = FALSE]
  grid <- expand.grid(
    promoter = row_order,
    compound = col_order,
    stringsAsFactors = FALSE
  )
  grid$hit <- as.logical(m[cbind(grid$promoter, grid$compound)])
  grid$effect <- as.numeric(effects[cbind(grid$promoter, grid$compound)])
  grid$hit_effect <- ifelse(grid$hit, grid$effect, NA_real_)
  grid$method <- method_label(method)
  grid
}

plot_data <- rbind(
  make_matrix(moderated, dstressr_method),
  make_matrix(median_polish, "median_polish")
)
plot_data$promoter <- factor(plot_data$promoter, levels = rev(row_order))
plot_data$compound <- factor(plot_data$compound, levels = col_order)
plot_data$method <- factor(plot_data$method, levels = c(
  method_label(dstressr_method),
  method_label("median_polish")
))

compound_order <- data.frame(
  compound = col_order,
  compound_order = seq_along(col_order),
  moderated_promoter_hits = colSums(moderated_matrix[row_order, col_order, drop = FALSE]),
  stringsAsFactors = FALSE
)
compound_order <- merge(compound_order, compound_lookup, by = "compound", all.x = TRUE, sort = FALSE)
compound_order <- compound_order[order(compound_order$compound_order), , drop = FALSE]

promoter_order <- data.frame(
  promoter = row_order,
  promoter_order = seq_along(row_order),
  moderated_compound_hits = rowSums(moderated_matrix[row_order, col_order, drop = FALSE]),
  stringsAsFactors = FALSE
)

summary <- aggregate(hit ~ method, plot_data, sum)
names(summary)[2] <- "hit_pairs"
summary$n_promoters <- length(row_order)
summary$n_compounds <- length(col_order)
summary$adjustment <- adjustment
effect_limit <- max(abs(plot_data$hit_effect), na.rm = TRUE)
if (!is.finite(effect_limit) || effect_limit == 0) {
  effect_limit <- 1
}

p <- ggplot(plot_data, aes(compound, promoter, fill = hit_effect)) +
  geom_tile(color = "white", linewidth = 0.04) +
  facet_grid(method ~ ., switch = "y") +
  scale_fill_gradient2(
    low = "#2563eb",
    mid = "#f8fafc",
    high = "#b91c1c",
    midpoint = 0,
    limits = c(-effect_limit, effect_limit),
    na.value = "#f1f5f9",
    name = "Effect size"
  ) +
  scale_x_discrete(labels = compound_axis_labels, drop = FALSE) +
  theme_bw(base_size = 9) +
  theme(
    axis.text.x = element_text(size = 7, angle = 55, hjust = 1, vjust = 1),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 7),
    strip.background = element_rect(fill = "#f8fafc", color = "#cbd5e1"),
    strip.text.y.left = element_text(face = "bold", angle = 90),
    panel.grid = element_blank(),
    panel.spacing.y = unit(0.25, "lines"),
    legend.position = "bottom",
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(8, 12, 70, 12)
  ) +
  labs(
    title = "Bipartite hit-effect heatmap",
    subtitle = paste0(
      "Compounds uniquely connected to PcmeA removed; columns sorted by DStressR default moderated in-degree; ",
      "tile color shows signed effect for ", adjustment, " BH hits"
    ),
    x = paste0("DStressR default moderated significant compounds after removing PcmeA-unique hits (n = ", length(col_order), ")"),
    y = "Promoter"
  )

write.table(
  plot_data,
  file.path(out_dir, "moderated_order_hit_matrix_long.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
write.table(
  promoter_order,
  file.path(out_dir, "moderated_order_promoters.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
write.table(
  compound_order,
  file.path(out_dir, "moderated_order_compounds.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
write.table(
  summary,
  file.path(out_dir, "moderated_order_hit_matrix_summary.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

ggsave(
  file.path(out_dir, "dstressr_default_moderated_order_hit_bipartite_heatmap.png"),
  p,
  width = 30,
  height = 12,
  dpi = 300,
  bg = "white"
)
ggsave(
  file.path(out_dir, "dstressr_default_moderated_order_hit_bipartite_heatmap.pdf"),
  p,
  width = 30,
  height = 12,
  bg = "white"
)
ggsave(
  file.path(out_dir, "moderated_order_hit_bipartite_heatmap.png"),
  p,
  width = 30,
  height = 12,
  dpi = 300,
  bg = "white"
)
ggsave(
  file.path(out_dir, "moderated_order_hit_bipartite_heatmap.pdf"),
  p,
  width = 30,
  height = 12,
  bg = "white"
)

message("Wrote bipartite hit heatmap to: ", out_dir)
print(summary)
