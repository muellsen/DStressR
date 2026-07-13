#!/usr/bin/env Rscript

source(file.path("analysis", "_helpers.R"))

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
  stop(
    "Missing explicit rejected pair table: ", input_file,
    "\nRun analysis/export_rejected_pair_list.R first.",
    call. = FALSE
  )
}

truthy <- function(x) {
  if (is.logical(x)) {
    return(x)
  }
  tolower(as.character(x)) %in% c("true", "t", "1", "yes")
}

unit_vector <- function(x, y) {
  len <- sqrt(x^2 + y^2)
  len[len == 0] <- 1
  data.frame(x = x / len, y = y / len)
}

set.seed(17)

hits <- read.delim(input_file, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
hits$median_polish <- truthy(hits$median_polish)
hits$destress_gaussian <- truthy(hits$destress_gaussian)
hits$destress_eb <- truthy(hits$destress_eb)
hits$compound_label <- ifelse(
  is.na(hits$ProductName) | hits$ProductName == "" | hits$ProductName == "NA",
  hits$srn_code,
  hits$ProductName
)

promoters <- sort(unique(hits$promoter))
compounds <- sort(unique(hits$srn_code))

incidence <- xtabs(~ promoter + srn_code, hits)
incidence <- incidence[promoters, compounds, drop = FALSE]
if (nrow(incidence) > 2) {
  promoter_dist <- stats::dist(as.matrix(incidence > 0), method = "binary")
  promoter_order <- rownames(incidence)[stats::hclust(promoter_dist, method = "average")$order]
} else {
  promoter_order <- promoters
}

theta <- seq(0, 2 * pi, length.out = length(promoter_order) + 1)[-length(promoter_order) - 1]
theta <- theta + pi / 2
promoter_counts <- aggregate(srn_code ~ promoter, hits, function(x) length(unique(x)))
names(promoter_counts)[2] <- "n_compounds"

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
promoter_nodes$n_compounds[is.na(promoter_nodes$n_compounds)] <- 0
promoter_nodes$degree <- promoter_nodes$n_compounds
promoter_nodes$node_size <- rescale(sqrt(promoter_nodes$n_compounds), to = c(4.5, 16))

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
      promoter_signature = paste(sort(ps), collapse = ";"),
      stringsAsFactors = FALSE
    )
  })
)

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
      # Put singleton compounds in the promoter's sector, closer to the edge.
      angle <- pxy$angle[1] + stats::runif(1, -0.12, 0.12)
      radius <- stats::runif(1, 0.62, 0.88)
      x <- radius * cos(angle)
      y <- radius * sin(angle)
    } else {
      # Shared compounds sit near the barycenter of the promoters they connect.
      # The more promoters they hit, the closer they move to the center.
      radius <- rescale(row$n_promoters, to = c(0.52, 0.18), from = range(compound_summary$n_promoters))
      if (!is.finite(radius)) {
        radius <- 0.42
      }
      jitter <- stats::runif(2, -0.025, 0.025)
      x <- radius * uv$x + jitter[1]
      y <- radius * uv$y + jitter[2]
    }
    data.frame(node_id = row$node_id, x = x, y = y)
  })
)

compound_nodes <- merge(compound_summary, compound_coords, by = "node_id", all.x = TRUE, sort = FALSE)
compound_nodes$degree <- compound_nodes$n_promoters
compound_nodes$node_size <- rescale(sqrt(compound_nodes$n_promoters), to = c(1.1, 4.5))

nodes <- rbind(
  promoter_nodes[, c("node_id", "label", "type", "x", "y", "degree", "node_size")],
  compound_nodes[, c("node_id", "label", "type", "x", "y", "degree", "node_size")]
)

method_long <- rbind(
  data.frame(
    promoter = hits$promoter,
    srn_code = hits$srn_code,
    method = "Median-polish",
    significant = hits$median_polish,
    stringsAsFactors = FALSE
  ),
  data.frame(
    promoter = hits$promoter,
    srn_code = hits$srn_code,
    method = "DStressR Gaussian",
    significant = hits$destress_gaussian,
    stringsAsFactors = FALSE
  ),
  data.frame(
    promoter = hits$promoter,
    srn_code = hits$srn_code,
    method = "DStressR EB",
    significant = hits$destress_eb,
    stringsAsFactors = FALSE
  )
)
method_long <- method_long[method_long$significant, ]
method_long$method <- factor(method_long$method, levels = c("Median-polish", "DStressR Gaussian", "DStressR EB"))
method_long$method_offset <- c("Median-polish" = -0.012, "DStressR Gaussian" = 0, "DStressR EB" = 0.012)[as.character(method_long$method)]

edge_base <- merge(
  method_long,
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

dx <- edge_base$x_compound - edge_base$x_promoter
dy <- edge_base$y_compound - edge_base$y_promoter
perp <- unit_vector(-dy, dx)
edge_base$x <- edge_base$x_promoter + edge_base$method_offset * perp$x
edge_base$y <- edge_base$y_promoter + edge_base$method_offset * perp$y
edge_base$xend <- edge_base$x_compound + edge_base$method_offset * perp$x
edge_base$yend <- edge_base$y_compound + edge_base$method_offset * perp$y

compound_label_nodes <- compound_nodes[
  compound_nodes$n_promoters >= 3 |
    compound_nodes$node_id %in% names(sort(table(hits$srn_code), decreasing = TRUE))[seq_len(min(18, length(unique(hits$srn_code))))],
]

method_colors <- c(
  "Median-polish" = "#64748b",
  "DStressR Gaussian" = "#2563eb",
  "DStressR EB" = "#16a34a"
)

p <- ggplot() +
  geom_curve(
    data = edge_base,
    aes(x = x, y = y, xend = xend, yend = yend, color = method),
    curvature = 0.08,
    linewidth = 0.28,
    alpha = 0.5
  ) +
  geom_point(
    data = subset(nodes, type == "compound"),
    aes(x = x, y = y, size = node_size),
    shape = 21,
    stroke = 0.2,
    fill = "#f8fafc",
    color = "#475569",
    alpha = 0.85
  ) +
  geom_point(
    data = subset(nodes, type == "promoter"),
    aes(x = x, y = y, size = node_size),
    shape = 21,
    stroke = 0.7,
    fill = "#fee2e2",
    color = "#991b1b"
  ) +
  geom_text_repel(
    data = promoter_nodes,
    aes(x = x * 1.06, y = y * 1.06, label = label),
    size = 3.1,
    segment.color = NA,
    max.overlaps = Inf,
    box.padding = 0.15,
    point.padding = 0.1
  ) +
  geom_text_repel(
    data = compound_label_nodes,
    aes(x = x, y = y, label = label),
    size = 2,
    color = "#334155",
    alpha = 0.85,
    max.overlaps = 80,
    box.padding = 0.12,
    point.padding = 0.08,
    min.segment.length = 0.02,
    segment.size = 0.15
  ) +
  scale_color_manual(values = method_colors) +
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
  guides(color = guide_legend(override.aes = list(linewidth = 1.2, alpha = 1))) +
  labs(
    title = "BH-rejected promoter-compound network",
    subtitle = "Promoter size: number of rejected compounds; compound placement: shared promoter barycenter",
    color = "Evidence method"
  )

write.table(
  nodes,
  file.path(out_dir, "bh_rejected_network_nodes.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
write.table(
  edge_base[, c("promoter", "srn_code", "method", "x", "y", "xend", "yend")],
  file.path(out_dir, "bh_rejected_network_edges.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

ggsave(
  file.path(out_dir, "bh_rejected_promoter_compound_network.png"),
  p,
  width = 14,
  height = 12,
  dpi = 300,
  bg = "white"
)
ggsave(
  file.path(out_dir, "bh_rejected_promoter_compound_network.pdf"),
  p,
  width = 14,
  height = 12,
  bg = "white"
)

message("Wrote network outputs to: ", out_dir)
message("Promoters: ", nrow(promoter_nodes), "; compounds: ", nrow(compound_nodes), "; method-specific edges: ", nrow(edge_base))
