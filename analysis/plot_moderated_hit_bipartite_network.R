#!/usr/bin/env Rscript

source(file.path("analysis", "_helpers.R"))

suppressPackageStartupMessages({
  library(ggplot2)
})

out_dir <- comparison_results_dir("hit_bipartite_network")
adjustment <- comparison_adjustment()
method <- "destress_standard"
padj_col <- padj_column(method, adjustment)

libmap <- read_tsv_base(libmap_path())
libmap$libplate <- paste0("lp", libmap[["Library plate"]])
libmap$compound <- paste(libmap$libplate, libmap[["Well"]], sep = "_")
libmap$ProductName <- ifelse(
  is.na(libmap$ProductName) | libmap$ProductName == "NA" | libmap$ProductName == "",
  libmap[["Catalog Number"]],
  libmap$ProductName
)
compound_lookup <- unique(libmap[, c("compound", "ProductName", "Catalog Number", "Target"), drop = FALSE])

tab <- read_package_pair_results(method)
tab$hit <- is.finite(tab[[padj_col]]) & tab[[padj_col]] < 0.05
hits <- tab[tab$hit, , drop = FALSE]
if (nrow(hits) == 0) {
  stop("No DStressR standard hits available for bipartite network.", call. = FALSE)
}

hits <- merge(hits, compound_lookup, by = "compound", all.x = TRUE, sort = FALSE)
hits$compound_label <- ifelse(
  is.na(hits$ProductName) | hits$ProductName == "" | hits$ProductName == "NA",
  hits$compound,
  hits$ProductName
)

promoters <- sort(unique(tab$promoter))
compounds <- sort(unique(hits$compound))
incidence <- xtabs(hit ~ promoter + compound, tab)
incidence <- incidence[promoters, compounds, drop = FALSE] > 0

compound_degree <- colSums(incidence)
col_order <- names(sort(compound_degree, decreasing = TRUE))
compound_rank <- stats::setNames(seq_along(col_order), col_order)
promoter_score <- vapply(rownames(incidence), function(promoter) {
  connected <- colnames(incidence)[incidence[promoter, ]]
  if (length(connected) == 0) {
    return(Inf)
  }
  stats::median(compound_rank[connected], na.rm = TRUE)
}, numeric(1))
row_order <- names(sort(promoter_score, decreasing = FALSE))

promoter_nodes <- data.frame(
  node_id = row_order,
  label = row_order,
  type = "promoter",
  x = 0,
  y = seq(1, 0, length.out = length(row_order)),
  stringsAsFactors = FALSE
)
promoter_nodes$n_hits <- rowSums(incidence[row_order, col_order, drop = FALSE])

compound_nodes <- data.frame(
  node_id = col_order,
  type = "compound",
  x = 1,
  y = seq(1, 0, length.out = length(col_order)),
  stringsAsFactors = FALSE
)
compound_nodes$n_hits <- colSums(incidence[row_order, col_order, drop = FALSE])
compound_nodes <- merge(compound_nodes, compound_lookup, by.x = "node_id", by.y = "compound", all.x = TRUE, sort = FALSE)
compound_nodes$label <- ifelse(
  is.na(compound_nodes$ProductName) | compound_nodes$ProductName == "" | compound_nodes$ProductName == "NA",
  compound_nodes$node_id,
  compound_nodes$ProductName
)

hit_promoters_by_compound <- split(hits$promoter, hits$compound)
compound_nodes$promoter_signature <- vapply(
  compound_nodes$node_id,
  function(compound) paste(sort(unique(hit_promoters_by_compound[[compound]])), collapse = ";"),
  character(1)
)
compound_nodes$pcmeA_only <- compound_nodes$promoter_signature == "PcmeA"

top_labeled_compounds <- head(col_order, 25)
compound_nodes$label_compound <- compound_nodes$node_id %in% top_labeled_compounds
compound_nodes$label_reason <- ifelse(
  compound_nodes$label_compound,
  "top_25_by_indegree",
  "unlabeled"
)

edge_base <- hits[, c("promoter", "compound", paste0(method, "_effect"), padj_col), drop = FALSE]
names(edge_base)[names(edge_base) == paste0(method, "_effect")] <- "effect"
names(edge_base)[names(edge_base) == padj_col] <- "padj"
edge_base <- merge(edge_base, promoter_nodes[, c("node_id", "x", "y")],
                   by.x = "promoter", by.y = "node_id", all.x = TRUE, sort = FALSE)
names(edge_base)[names(edge_base) == "x"] <- "x_promoter"
names(edge_base)[names(edge_base) == "y"] <- "y_promoter"
edge_base <- merge(edge_base, compound_nodes[, c("node_id", "x", "y", "label")],
                   by.x = "compound", by.y = "node_id", all.x = TRUE, sort = FALSE)
names(edge_base)[names(edge_base) == "x"] <- "x_compound"
names(edge_base)[names(edge_base) == "y"] <- "y_compound"
edge_base$direction <- ifelse(edge_base$effect >= 0, "Up", "Down")

edge_plot <- edge_base
edge_plot$x <- edge_plot$x_promoter + 0.02
edge_plot$xend <- edge_plot$x_compound - 0.02
edge_plot$alpha <- pmin(0.6, pmax(0.09, -log10(pmax(edge_plot$padj, .Machine$double.xmin)) / 15))

compound_label_nodes <- compound_nodes[compound_nodes$label_compound, , drop = FALSE]
compound_label_nodes$label_short <- iconv(compound_label_nodes$label, from = "", to = "ASCII//TRANSLIT")
compound_label_nodes$label_short[is.na(compound_label_nodes$label_short)] <- compound_label_nodes$label[is.na(compound_label_nodes$label_short)]
compound_label_nodes$label_short <- gsub("\u03b3", "gamma", compound_label_nodes$label_short, fixed = TRUE)
too_long <- nchar(compound_label_nodes$label_short) > 42
compound_label_nodes$label_short[too_long] <- paste0(substr(compound_label_nodes$label_short[too_long], 1, 39), "...")
compound_label_nodes$label_x <- 1.12
compound_label_nodes$label_y <- seq(0.985, 0.025, length.out = nrow(compound_label_nodes))
compound_label_nodes$rank_label <- paste0(seq_len(nrow(compound_label_nodes)), ". ", compound_label_nodes$label_short)

p <- ggplot() +
  geom_segment(
    data = edge_plot,
    aes(x = x, y = y_promoter, xend = xend, yend = y_compound, color = direction, alpha = alpha),
    linewidth = 0.22,
    lineend = "round"
  ) +
  geom_segment(
    data = compound_label_nodes,
    aes(x = x + 0.012, y = y, xend = label_x - 0.012, yend = label_y),
    color = "#64748b",
    linewidth = 0.18,
    alpha = 0.55
  ) +
  geom_point(
    data = promoter_nodes,
    aes(x = x, y = y, size = n_hits),
    shape = 21,
    fill = "#fee2e2",
    color = "#991b1b",
    stroke = 0.5
  ) +
  geom_point(
    data = compound_nodes,
    aes(x = x, y = y, size = n_hits),
    shape = 21,
    fill = "#e0f2fe",
    color = "#0369a1",
    stroke = 0.25
  ) +
  geom_text(
    data = promoter_nodes,
    aes(x = x - 0.02, y = y, label = label),
    hjust = 1,
    size = 2.6,
    fontface = "bold"
  ) +
  geom_text(
    data = compound_label_nodes,
    aes(x = label_x, y = label_y, label = rank_label),
    hjust = 0,
    size = 2.45,
    color = "#0f172a"
  ) +
  scale_color_manual(values = c("Down" = "#2563eb", "Up" = "#b91c1c")) +
  scale_alpha_identity() +
  scale_size_continuous(range = c(1.2, 8), guide = "none") +
  coord_cartesian(xlim = c(-0.18, 1.62), ylim = c(-0.02, 1.02), clip = "off") +
  theme_void(base_size = 10) +
  theme(
    legend.position = "bottom",
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(12, 230, 12, 72)
  ) +
  labs(
    title = "DStressR standard model bipartite hit network",
    subtitle = paste0(
      "Compounds ranked top-to-bottom by in-degree; promoters ordered by connected compound ranks; ",
      adjustment, " BH, FDR < 0.05. Top 25 compound names shown."
    ),
    color = "Effect"
  )

write.table(
  edge_base,
  file.path(out_dir, "moderated_bipartite_network_edges.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
write.table(
  promoter_nodes,
  file.path(out_dir, "moderated_bipartite_network_promoters.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
write.table(
  compound_nodes,
  file.path(out_dir, "moderated_bipartite_network_compounds.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

ggsave(
  file.path(out_dir, "dstressr_standard_bipartite_network.png"),
  p,
  width = 16,
  height = 12,
  dpi = 300,
  bg = "white",
  limitsize = FALSE
)
ggsave(
  file.path(out_dir, "dstressr_standard_bipartite_network.pdf"),
  p,
  width = 16,
  height = 12,
  bg = "white",
  limitsize = FALSE
)
ggsave(
  file.path(out_dir, "moderated_bipartite_network.png"),
  p,
  width = 16,
  height = 12,
  dpi = 300,
  bg = "white",
  limitsize = FALSE
)
ggsave(
  file.path(out_dir, "moderated_bipartite_network.pdf"),
  p,
  width = 16,
  height = 12,
  bg = "white",
  limitsize = FALSE
)

message("Wrote DStressR standard model bipartite network to: ", out_dir)
message("Promoters: ", nrow(promoter_nodes), "; compounds: ", nrow(compound_nodes), "; edges: ", nrow(edge_base))
message("Labeled compounds: ", sum(compound_nodes$label_compound), " / ", nrow(compound_nodes))
