#!/usr/bin/env Rscript

source(file.path("analysis", "_helpers.R"))

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(scales)
})

out_dir <- comparison_results_dir("hit_network")
adjustment <- comparison_adjustment()
method <- "destress_moderated"
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
  stop("No DStressR default moderated hits available for network.", call. = FALSE)
}

compound_promoters <- split(hits$promoter, hits$compound)
pcmeA_unique_compounds <- names(compound_promoters)[vapply(
  compound_promoters,
  function(x) identical(sort(unique(x)), "PcmeA"),
  logical(1)
)]
hits <- hits[!hits$compound %in% pcmeA_unique_compounds, , drop = FALSE]
if (nrow(hits) == 0) {
  stop("No hits remain after removing PcmeA-unique compounds.", call. = FALSE)
}

hits <- merge(hits, compound_lookup, by = "compound", all.x = TRUE, sort = FALSE)
hits$compound_label <- ifelse(
  is.na(hits$ProductName) | hits$ProductName == "" | hits$ProductName == "NA",
  hits$compound,
  hits$ProductName
)
hits$direction <- ifelse(hits[[paste0(method, "_effect")]] >= 0, "Up", "Down")

promoters <- sort(unique(hits$promoter))
compounds <- sort(unique(hits$compound))
incidence <- xtabs(~ promoter + compound, hits)
incidence <- incidence[promoters, compounds, drop = FALSE]

promoter_order <- if (nrow(incidence) > 2 && ncol(incidence) > 1) {
  rownames(incidence)[stats::hclust(stats::dist(as.matrix(incidence > 0), method = "binary"), method = "average")$order]
} else {
  promoters
}

theta <- seq(0, 2 * pi, length.out = length(promoter_order) + 1)[-length(promoter_order) - 1]
theta <- theta + pi / 2
promoter_nodes <- data.frame(
  node_id = promoter_order,
  label = promoter_order,
  type = "promoter",
  x = cos(theta),
  y = sin(theta),
  stringsAsFactors = FALSE
)
promoter_nodes$degree <- rowSums(incidence[promoter_nodes$node_id, , drop = FALSE] > 0)

compound_summary <- do.call(rbind, lapply(split(hits, hits$compound), function(d) {
  ps <- sort(unique(d$promoter))
  data.frame(
    node_id = d$compound[1],
    label = d$compound_label[1],
    type = "compound",
    degree = length(ps),
    promoter_signature = paste(ps, collapse = ";"),
    min_padj = min(d[[padj_col]], na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))

set.seed(17)
compound_coords <- do.call(rbind, lapply(seq_len(nrow(compound_summary)), function(i) {
  row <- compound_summary[i, ]
  ps <- unlist(strsplit(row$promoter_signature, ";", fixed = TRUE))
  pxy <- promoter_nodes[match(ps, promoter_nodes$node_id), c("x", "y"), drop = FALSE]
  center_x <- mean(pxy$x)
  center_y <- mean(pxy$y)
  len <- sqrt(center_x^2 + center_y^2)
  if (!is.finite(len) || len == 0) {
    center_x <- stats::runif(1, -0.15, 0.15)
    center_y <- stats::runif(1, -0.15, 0.15)
    len <- sqrt(center_x^2 + center_y^2)
  }
  ux <- center_x / len
  uy <- center_y / len
  radius <- if (row$degree == 1) {
    0.72
  } else {
    scales::rescale(row$degree, to = c(0.52, 0.18), from = range(compound_summary$degree))
  }
  if (!is.finite(radius)) {
    radius <- 0.42
  }
  jitter <- stats::runif(2, -0.035, 0.035)
  data.frame(
    node_id = row$node_id,
    x = radius * ux + jitter[1],
    y = radius * uy + jitter[2]
  )
}))

compound_nodes <- merge(compound_summary, compound_coords, by = "node_id", all.x = TRUE, sort = FALSE)
compound_nodes$label_plot <- iconv(compound_nodes$label, from = "", to = "ASCII//TRANSLIT")
compound_nodes$label_plot[is.na(compound_nodes$label_plot)] <- compound_nodes$label[is.na(compound_nodes$label_plot)]
compound_nodes$label_plot <- gsub("\u03b3", "gamma", compound_nodes$label_plot, fixed = TRUE)
too_long <- nchar(compound_nodes$label_plot) > 36
compound_nodes$label_plot[too_long] <- paste0(substr(compound_nodes$label_plot[too_long], 1, 33), "...")

top_compounds <- compound_nodes[order(-compound_nodes$degree, compound_nodes$min_padj, compound_nodes$label), ]
top_compounds <- head(top_compounds, 25)
top_compounds$label_side <- ifelse(top_compounds$x >= stats::median(top_compounds$x), "right", "left")
if (length(unique(top_compounds$label_side)) == 1 && nrow(top_compounds) > 1) {
  top_compounds$label_side[order(top_compounds$x)[seq_len(floor(nrow(top_compounds) / 2))]] <- "left"
  top_compounds$label_side[top_compounds$label_side != "left"] <- "right"
}
top_compounds <- do.call(rbind, lapply(split(top_compounds, top_compounds$label_side), function(d) {
  d <- d[order(-d$degree, d$min_padj, d$label_plot), , drop = FALSE]
  d$label_x <- ifelse(d$label_side[1] == "right", 1.36, -1.36)
  d$label_y <- seq(0.86, -0.86, length.out = nrow(d))
  d$label_hjust <- ifelse(d$label_side[1] == "right", 0, 1)
  d
}))

nodes <- rbind(
  promoter_nodes[, c("node_id", "label", "type", "x", "y", "degree")],
  data.frame(
    node_id = compound_nodes$node_id,
    label = compound_nodes$label,
    type = compound_nodes$type,
    x = compound_nodes$x,
    y = compound_nodes$y,
    degree = compound_nodes$degree,
    stringsAsFactors = FALSE
  )
)

edge_base <- merge(hits, promoter_nodes[, c("node_id", "x", "y")],
                   by.x = "promoter", by.y = "node_id", all.x = TRUE, sort = FALSE)
names(edge_base)[names(edge_base) == "x"] <- "x_promoter"
names(edge_base)[names(edge_base) == "y"] <- "y_promoter"
edge_base <- merge(edge_base, compound_nodes[, c("node_id", "x", "y")],
                   by.x = "compound", by.y = "node_id", all.x = TRUE, sort = FALSE)
names(edge_base)[names(edge_base) == "x"] <- "x_compound"
names(edge_base)[names(edge_base) == "y"] <- "y_compound"
edge_base$alpha <- pmin(0.55, pmax(0.12, -log10(pmax(edge_base[[padj_col]], .Machine$double.xmin)) / 15))

p <- ggplot() +
  geom_curve(
    data = edge_base,
    aes(x = x_promoter, y = y_promoter, xend = x_compound, yend = y_compound,
        color = direction, alpha = alpha),
    curvature = 0.07,
    linewidth = 0.26
  ) +
  geom_point(
    data = subset(nodes, type == "compound"),
    aes(x = x, y = y, size = degree),
    shape = 21,
    stroke = 0.25,
    fill = "#e0f2fe",
    color = "#0369a1"
  ) +
  geom_point(
    data = subset(nodes, type == "promoter"),
    aes(x = x, y = y, size = degree),
    shape = 21,
    stroke = 0.65,
    fill = "#fee2e2",
    color = "#991b1b"
  ) +
  geom_text_repel(
    data = promoter_nodes,
    aes(x = x * 1.08, y = y * 1.08, label = label),
    size = 3,
    fontface = "bold",
    segment.color = NA,
    max.overlaps = Inf,
    box.padding = 0.12,
    point.padding = 0.08
  ) +
  geom_segment(
    data = top_compounds,
    aes(x = x, y = y, xend = label_x, yend = label_y),
    color = "#475569",
    linewidth = 0.34,
    alpha = 0.95
  ) +
  geom_point(
    data = top_compounds,
    aes(x = x, y = y),
    shape = 21,
    size = 2.1,
    stroke = 0.55,
    fill = "white",
    color = "#0f172a"
  ) +
  geom_label(
    data = top_compounds,
    aes(x = label_x, y = label_y, label = label_plot, hjust = label_hjust),
    size = 2.75,
    color = "#0f172a",
    label.size = 0,
    fill = alpha("white", 0.82),
    label.padding = unit(0.08, "lines"),
    lineheight = 0.92
  ) +
  scale_color_manual(values = c("Down" = "#2563eb", "Up" = "#b91c1c")) +
  scale_alpha_identity() +
  scale_size_continuous(range = c(1.2, 13), guide = "none") +
  coord_equal(xlim = c(-1.65, 1.65), ylim = c(-1.18, 1.18), clip = "off") +
  theme_void(base_size = 10) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(12, 36, 12, 36)
  ) +
  guides(color = guide_legend(override.aes = list(linewidth = 1.2, alpha = 1))) +
  labs(
    title = "DStressR default (moderated) hit network",
    subtitle = paste0(
      "Compounds uniquely connected to PcmeA removed; node size proportional to degree; ",
      "top 25 compounds labeled"
    ),
    color = "Effect"
  )

write.table(
  edge_base,
  file.path(out_dir, "dstressr_default_moderated_network_no_pcmeA_unique_edges.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
write.table(
  nodes,
  file.path(out_dir, "dstressr_default_moderated_network_no_pcmeA_unique_nodes.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
write.table(
  top_compounds,
  file.path(out_dir, "dstressr_default_moderated_network_no_pcmeA_unique_top25_compounds.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
write.table(
  edge_base,
  file.path(out_dir, "moderated_network_no_pcmeA_unique_edges.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
write.table(
  nodes,
  file.path(out_dir, "moderated_network_no_pcmeA_unique_nodes.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
write.table(
  top_compounds,
  file.path(out_dir, "moderated_network_no_pcmeA_unique_top25_compounds.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

ggsave(
  file.path(out_dir, "dstressr_default_moderated_network_no_pcmeA_unique.png"),
  p,
  width = 12,
  height = 10,
  dpi = 300,
  bg = "white"
)
ggsave(
  file.path(out_dir, "dstressr_default_moderated_network_no_pcmeA_unique.pdf"),
  p,
  width = 12,
  height = 10,
  bg = "white"
)
ggsave(
  file.path(out_dir, "moderated_network_no_pcmeA_unique.png"),
  p,
  width = 12,
  height = 10,
  dpi = 300,
  bg = "white"
)
ggsave(
  file.path(out_dir, "moderated_network_no_pcmeA_unique.pdf"),
  p,
  width = 12,
  height = 10,
  bg = "white"
)

message("Wrote DStressR default moderated network without PcmeA-unique compounds to: ", out_dir)
message("Removed PcmeA-unique compounds: ", length(pcmeA_unique_compounds))
message("Promoters: ", sum(nodes$type == "promoter"), "; compounds: ", sum(nodes$type == "compound"), "; edges: ", nrow(edge_base))
