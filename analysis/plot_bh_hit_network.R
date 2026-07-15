#!/usr/bin/env Rscript

source(file.path("analysis", "_helpers.R"))

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(scales)
})

methods <- c("median_polish", "destress_standard", "destress_moderated")
input_file <- file.path(comparison_results_dir("hit_overlap"), "differential_pair_list.tsv")
out_dir <- comparison_results_dir("hit_network")

if (!file.exists(input_file)) {
  stop(
    "Missing differential pair list: ", input_file,
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

hits <- read_tsv_base(input_file)
for (method in methods) {
  hits[[method]] <- truthy(hits[[method]])
}
hits$compound_label <- if ("ProductName" %in% names(hits)) {
  ifelse(is.na(hits$ProductName) | hits$ProductName == "" | hits$ProductName == "NA",
         hits$compound, hits$ProductName)
} else {
  hits$compound
}

promoters <- sort(unique(hits$promoter))
compounds <- sort(unique(hits$compound))

incidence <- xtabs(~ promoter + compound, hits)
incidence <- incidence[promoters, compounds, drop = FALSE]
promoter_order <- if (nrow(incidence) > 2) {
  rownames(incidence)[stats::hclust(stats::dist(as.matrix(incidence > 0), method = "binary"), method = "average")$order]
} else {
  promoters
}

theta <- seq(0, 2 * pi, length.out = length(promoter_order) + 1)[-length(promoter_order) - 1]
theta <- theta + pi / 2
promoter_counts <- aggregate(compound ~ promoter, hits, function(x) length(unique(x)))
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

compound_summary <- do.call(rbind, lapply(split(hits, hits$compound), function(d) {
  ps <- unique(d$promoter)
  data.frame(
    node_id = d$compound[1],
    label = d$compound_label[1],
    type = "compound",
    n_promoters = length(ps),
    promoter_signature = paste(sort(ps), collapse = ";"),
    stringsAsFactors = FALSE
  )
}))

set.seed(17)
compound_coords <- do.call(rbind, lapply(seq_len(nrow(compound_summary)), function(i) {
  row <- compound_summary[i, ]
  ps <- unlist(strsplit(row$promoter_signature, ";", fixed = TRUE))
  pxy <- promoter_nodes[match(ps, promoter_nodes$node_id), c("x", "y", "angle"), drop = FALSE]
  center_x <- mean(pxy$x)
  center_y <- mean(pxy$y)
  uv <- unit_vector(center_x, center_y)
  if (row$n_promoters == 1) {
    angle <- pxy$angle[1] + stats::runif(1, -0.12, 0.12)
    radius <- stats::runif(1, 0.62, 0.88)
    x <- radius * cos(angle)
    y <- radius * sin(angle)
  } else {
    radius <- rescale(row$n_promoters, to = c(0.52, 0.18), from = range(compound_summary$n_promoters))
    if (!is.finite(radius)) {
      radius <- 0.42
    }
    jitter <- stats::runif(2, -0.025, 0.025)
    x <- radius * uv$x + jitter[1]
    y <- radius * uv$y + jitter[2]
  }
  data.frame(node_id = row$node_id, x = x, y = y)
}))

compound_nodes <- merge(compound_summary, compound_coords, by = "node_id", all.x = TRUE, sort = FALSE)
compound_nodes$degree <- compound_nodes$n_promoters
compound_nodes$node_size <- rescale(sqrt(compound_nodes$n_promoters), to = c(1.1, 4.5))

nodes <- rbind(
  promoter_nodes[, c("node_id", "label", "type", "x", "y", "degree", "node_size")],
  compound_nodes[, c("node_id", "label", "type", "x", "y", "degree", "node_size")]
)

method_long <- do.call(rbind, lapply(methods, function(method) {
  data.frame(
    promoter = hits$promoter,
    compound = hits$compound,
    method = method_label(method),
    significant = hits[[method]],
    stringsAsFactors = FALSE
  )
}))
method_long <- method_long[method_long$significant, , drop = FALSE]
method_levels <- vapply(methods, method_label, character(1))
method_long$method <- factor(method_long$method, levels = method_levels)
offsets <- stats::setNames(seq(-0.012, 0.012, length.out = length(methods)), method_levels)
method_long$method_offset <- offsets[as.character(method_long$method)]

edge_base <- merge(method_long, promoter_nodes[, c("node_id", "x", "y")],
                   by.x = "promoter", by.y = "node_id", all.x = TRUE, sort = FALSE)
names(edge_base)[names(edge_base) == "x"] <- "x_promoter"
names(edge_base)[names(edge_base) == "y"] <- "y_promoter"
edge_base <- merge(edge_base, compound_nodes[, c("node_id", "x", "y")],
                   by.x = "compound", by.y = "node_id", all.x = TRUE, sort = FALSE)
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
    compound_nodes$node_id %in% names(sort(table(hits$compound), decreasing = TRUE))[seq_len(min(18, length(unique(hits$compound))))],
]

method_colors <- stats::setNames(c("#64748b", "#2563eb", "#16a34a")[seq_along(methods)], method_levels)

p <- ggplot() +
  geom_curve(data = edge_base,
             aes(x = x, y = y, xend = xend, yend = yend, color = method),
             curvature = 0.08, linewidth = 0.28, alpha = 0.5) +
  geom_point(data = subset(nodes, type == "compound"),
             aes(x = x, y = y, size = node_size),
             shape = 21, stroke = 0.2, fill = "#f8fafc", color = "#475569", alpha = 0.85) +
  geom_point(data = subset(nodes, type == "promoter"),
             aes(x = x, y = y, size = node_size),
             shape = 21, stroke = 0.7, fill = "#fee2e2", color = "#991b1b") +
  geom_text_repel(data = promoter_nodes,
                  aes(x = x * 1.06, y = y * 1.06, label = label),
                  size = 3.1, segment.color = NA, max.overlaps = Inf,
                  box.padding = 0.15, point.padding = 0.1) +
  geom_text_repel(data = compound_label_nodes,
                  aes(x = x, y = y, label = label),
                  size = 2, color = "#334155", alpha = 0.85,
                  max.overlaps = 80, box.padding = 0.12, point.padding = 0.08,
                  min.segment.length = 0.02, segment.size = 0.15) +
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
    title = "Differential promoter-compound network",
    subtitle = "Package outputs; adjusted p < 0.05",
    color = "Evidence method"
  )

write.table(nodes, file.path(out_dir, "hit_network_nodes.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(edge_base[, c("promoter", "compound", "method", "x", "y", "xend", "yend")],
            file.path(out_dir, "hit_network_edges.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

ggsave(file.path(out_dir, "hit_network.png"), p, width = 14, height = 12, dpi = 300, bg = "white")
ggsave(file.path(out_dir, "hit_network.pdf"), p, width = 14, height = 12, bg = "white")

message("Wrote network outputs to: ", out_dir)
message("Promoters: ", nrow(promoter_nodes), "; compounds: ", nrow(compound_nodes), "; method-specific edges: ", nrow(edge_base))
