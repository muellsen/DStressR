#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(scales)
})

out_dir <- file.path(getwd(), "analysis", "outputs", "network")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

input_file <- file.path(
  getwd(),
  "analysis",
  "outputs",
  "venn",
  "bh_rejected_compound_promoter_pairs_explicit.tsv"
)

if (!file.exists(input_file)) {
  stop("Missing explicit rejected pair table. Run analysis/export_rejected_pair_list.R first.", call. = FALSE)
}

truthy <- function(x) {
  if (is.logical(x)) return(x)
  tolower(as.character(x)) %in% c("true", "t", "1", "yes")
}

unit_vector <- function(x, y) {
  len <- sqrt(x^2 + y^2)
  len[len == 0] <- 1
  data.frame(x = x / len, y = y / len)
}

compound_status <- function(median_any, eb_any) {
  ifelse(
    median_any & eb_any,
    "Both",
    ifelse(eb_any, "EB only", "Median-polish only")
  )
}

set.seed(71)

raw <- read.delim(input_file, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
raw$median_polish <- truthy(raw$median_polish)
raw$destress_eb <- truthy(raw$destress_eb)
raw$compound_label <- ifelse(
  is.na(raw$ProductName) | raw$ProductName == "" | raw$ProductName == "NA",
  raw$srn_code,
  raw$ProductName
)

hits <- raw[raw$median_polish | raw$destress_eb, ]
hits$edge_status <- ifelse(
  hits$median_polish & hits$destress_eb,
  "Both",
  ifelse(hits$destress_eb, "EB only", "Median-polish only")
)
hits$edge_status <- factor(hits$edge_status, levels = c("Both", "EB only", "Median-polish only"))

promoters <- sort(unique(hits$promoter))
compounds <- sort(unique(hits$srn_code))

incidence <- xtabs(~ promoter + srn_code, hits)
incidence <- incidence[promoters, compounds, drop = FALSE]
if (nrow(incidence) > 2 && ncol(incidence) > 1) {
  promoter_dist <- stats::dist(as.matrix(incidence > 0), method = "binary")
  promoter_order <- rownames(incidence)[stats::hclust(promoter_dist, method = "average")$order]
} else {
  promoter_order <- promoters
}

theta <- seq(0, 2 * pi, length.out = length(promoter_order) + 1)[-length(promoter_order) - 1]
theta <- theta + pi / 2

promoter_counts <- do.call(
  rbind,
  lapply(split(hits, hits$promoter), function(d) {
    data.frame(
      promoter = d$promoter[1],
      eb_compounds = length(unique(d$srn_code[d$destress_eb])),
      median_compounds = length(unique(d$srn_code[d$median_polish])),
      union_compounds = length(unique(d$srn_code)),
      stringsAsFactors = FALSE
    )
  })
)

promoter_nodes <- data.frame(
  node_id = promoter_order,
  label = promoter_order,
  type = "promoter",
  angle = theta,
  x = cos(theta),
  y = sin(theta),
  stringsAsFactors = FALSE
)
promoter_nodes <- merge(promoter_nodes, promoter_counts, by.x = "node_id", by.y = "promoter", all.x = TRUE, sort = FALSE)
promoter_nodes[is.na(promoter_nodes)] <- 0
promoter_nodes$degree <- promoter_nodes$union_compounds
promoter_nodes$node_size <- rescale(sqrt(promoter_nodes$union_compounds), to = c(5, 17))

compound_summary <- do.call(
  rbind,
  lapply(split(hits, hits$srn_code), function(d) {
    ps <- unique(d$promoter)
    data.frame(
      node_id = d$srn_code[1],
      label = d$compound_label[1],
      type = "compound",
      n_promoters = length(ps),
      n_edges = nrow(d),
      median_any = any(d$median_polish),
      eb_any = any(d$destress_eb),
      promoter_signature = paste(sort(ps), collapse = ";"),
      stringsAsFactors = FALSE
    )
  })
)
compound_summary$status <- compound_status(compound_summary$median_any, compound_summary$eb_any)
compound_summary$status <- factor(compound_summary$status, levels = c("Both", "EB only", "Median-polish only"))

compound_coords <- do.call(
  rbind,
  lapply(seq_len(nrow(compound_summary)), function(i) {
    row <- compound_summary[i, ]
    ps <- unlist(strsplit(row$promoter_signature, ";", fixed = TRUE))
    pxy <- promoter_nodes[match(ps, promoter_nodes$node_id), c("x", "y", "angle"), drop = FALSE]
    center_x <- mean(pxy$x)
    center_y <- mean(pxy$y)
    uv <- unit_vector(center_x, center_y)

    if (row$n_promoters == 1) {
      angle <- pxy$angle[1] + stats::runif(1, -0.1, 0.1)
      # Put singleton compounds on method-specific lanes around the promoter.
      radius <- switch(
        as.character(row$status),
        "Both" = stats::runif(1, 0.54, 0.64),
        "EB only" = stats::runif(1, 0.66, 0.82),
        "Median-polish only" = stats::runif(1, 0.42, 0.52)
      )
      x <- radius * cos(angle)
      y <- radius * sin(angle)
    } else {
      radius <- rescale(row$n_promoters, to = c(0.5, 0.18), from = range(compound_summary$n_promoters))
      if (!is.finite(radius)) radius <- 0.42
      status_shift <- switch(
        as.character(row$status),
        "Both" = 0,
        "EB only" = 0.04,
        "Median-polish only" = -0.04
      )
      jitter <- stats::runif(2, -0.018, 0.018)
      x <- (radius + status_shift) * uv$x + jitter[1]
      y <- (radius + status_shift) * uv$y + jitter[2]
    }
    data.frame(node_id = row$node_id, x = x, y = y)
  })
)

compound_nodes <- merge(compound_summary, compound_coords, by = "node_id", all.x = TRUE, sort = FALSE)
compound_nodes$degree <- compound_nodes$n_promoters
compound_nodes$node_size <- rescale(sqrt(compound_nodes$n_promoters), to = c(1.2, 5))

nodes <- rbind(
  promoter_nodes[, c("node_id", "label", "type", "x", "y", "degree", "node_size")],
  transform(
    compound_nodes[, c("node_id", "label", "type", "x", "y", "degree", "node_size")],
    label = as.character(label)
  )
)

edge_base <- merge(
  hits[, c("promoter", "srn_code", "edge_status", "median_polish", "destress_eb")],
  promoter_nodes[, c("node_id", "x", "y")],
  by.x = "promoter",
  by.y = "node_id",
  all.x = TRUE,
  sort = FALSE
)
names(edge_base)[names(edge_base) == "x"] <- "x_promoter"
names(edge_base)[names(edge_base) == "y"] <- "y_promoter"
edge_base <- merge(
  edge_base,
  compound_nodes[, c("node_id", "x", "y")],
  by.x = "srn_code",
  by.y = "node_id",
  all.x = TRUE,
  sort = FALSE
)
names(edge_base)[names(edge_base) == "x"] <- "x_compound"
names(edge_base)[names(edge_base) == "y"] <- "y_compound"

edge_colors <- c(
  "Both" = "#4c1d95",
  "EB only" = "#16a34a",
  "Median-polish only" = "#dc2626"
)

compound_fills <- c(
  "Both" = "#ede9fe",
  "EB only" = "#dcfce7",
  "Median-polish only" = "#fee2e2"
)

compound_label_nodes <- compound_nodes[
  compound_nodes$n_promoters >= 3 |
    compound_nodes$status == "Both" & compound_nodes$n_promoters >= 2 |
    compound_nodes$node_id %in% names(sort(table(hits$srn_code), decreasing = TRUE))[seq_len(min(30, length(unique(hits$srn_code))))],
]

status_counts <- as.data.frame(table(edge_base$edge_status), stringsAsFactors = FALSE)
names(status_counts) <- c("edge_status", "n_edges")

p <- ggplot() +
  geom_curve(
    data = edge_base,
    aes(x = x_promoter, y = y_promoter, xend = x_compound, yend = y_compound, color = edge_status),
    curvature = 0.08,
    linewidth = 0.4,
    alpha = 0.55
  ) +
  geom_point(
    data = compound_nodes,
    aes(x = x, y = y, size = node_size, fill = status),
    shape = 21,
    stroke = 0.25,
    color = "#475569",
    alpha = 0.92
  ) +
  geom_point(
    data = promoter_nodes,
    aes(x = x, y = y, size = node_size),
    shape = 21,
    stroke = 0.85,
    fill = "#f8fafc",
    color = "#111827"
  ) +
  geom_text_repel(
    data = promoter_nodes,
    aes(x = x * 1.06, y = y * 1.06, label = label),
    size = 3.2,
    segment.color = NA,
    max.overlaps = Inf,
    box.padding = 0.15,
    point.padding = 0.1
  ) +
  geom_text_repel(
    data = compound_label_nodes,
    aes(x = x, y = y, label = label),
    size = 2.0,
    color = "#334155",
    alpha = 0.9,
    max.overlaps = 90,
    box.padding = 0.12,
    point.padding = 0.08,
    min.segment.length = 0.02,
    segment.size = 0.15
  ) +
  scale_color_manual(values = edge_colors, drop = FALSE) +
  scale_fill_manual(values = compound_fills, drop = FALSE) +
  scale_size_identity() +
  coord_equal(xlim = c(-1.25, 1.25), ylim = c(-1.18, 1.18), clip = "off") +
  theme_void(base_size = 10) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  guides(
    color = guide_legend(override.aes = list(linewidth = 1.4, alpha = 1)),
    fill = guide_legend(override.aes = list(size = 4))
  ) +
  labs(
    title = "Median-polish vs DStressR EB BH-rejected network",
    subtitle = "Purple: agreement; green: EB additions; red: median-polish-only removals",
    color = "Pair status",
    fill = "Compound status"
  )

write.table(
  nodes,
  file.path(out_dir, "eb_vs_medianpolish_difference_network_nodes.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
write.table(
  edge_base[, c("promoter", "srn_code", "edge_status", "median_polish", "destress_eb", "x_promoter", "y_promoter", "x_compound", "y_compound")],
  file.path(out_dir, "eb_vs_medianpolish_difference_network_edges.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
write.table(
  status_counts,
  file.path(out_dir, "eb_vs_medianpolish_difference_network_edge_counts.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

ggsave(
  file.path(out_dir, "eb_vs_medianpolish_difference_network.png"),
  p,
  width = 14,
  height = 12,
  dpi = 300,
  bg = "white"
)
ggsave(
  file.path(out_dir, "eb_vs_medianpolish_difference_network.pdf"),
  p,
  width = 14,
  height = 12,
  bg = "white"
)

message("Wrote EB vs median-polish difference network outputs to: ", out_dir)
message(
  "Union pairs: ", nrow(edge_base),
  "; shared: ", sum(edge_base$edge_status == "Both"),
  "; EB only: ", sum(edge_base$edge_status == "EB only"),
  "; median-polish only: ", sum(edge_base$edge_status == "Median-polish only")
)
