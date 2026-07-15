#!/usr/bin/env Rscript

source(file.path("analysis", "_helpers.R"))

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(grid)
  library(gridExtra)
})

out_dir <- comparison_results_dir("significance_summary")
adjustment <- comparison_adjustment()
methods <- c("destress_standard", "median_polish")
method_labels <- c(
  destress_standard = "DStressR standard model",
  median_polish = "Median-polish max-p model"
)
table_labels <- c(
  destress_standard = "DStressR standard\nmodel",
  median_polish = "Median-polish\nmax-p model"
)

libmap <- read_tsv_base(libmap_path())
libmap$libplate <- paste0("lp", libmap[["Library plate"]])
libmap$compound <- paste(libmap$libplate, libmap[["Well"]], sep = "_")
libmap$ProductName <- ifelse(
  is.na(libmap$ProductName) | libmap$ProductName == "NA" | libmap$ProductName == "",
  libmap[["Catalog Number"]],
  libmap$ProductName
)
compound_lookup <- unique(libmap[, c("compound", "ProductName", "Catalog Number", "Target"), drop = FALSE])

tabs <- lapply(methods, read_package_pair_results)
names(tabs) <- methods
pair_table <- Reduce(
  function(x, y) merge(x, y, by = c("promoter", "compound", "pair_id"), all = FALSE, sort = FALSE),
  tabs
)
pair_table <- merge(pair_table, compound_lookup, by = "compound", all.x = TRUE, sort = FALSE)
pair_table$compound_label <- ifelse(
  is.na(pair_table$ProductName) | pair_table$ProductName == "" | pair_table$ProductName == "NA",
  pair_table$compound,
  pair_table$ProductName
)

for (method in methods) {
  padj <- padj_column(method, adjustment)
  pair_table[[paste0(method, "_hit")]] <- is.finite(pair_table[[padj]]) & pair_table[[padj]] < 0.05
  pair_table[[paste0(method, "_direction")]] <- ifelse(
    pair_table[[paste0(method, "_effect")]] > 0,
    "Positive",
    "Negative"
  )
}

highlight_file <- file.path(
  comparison_results_dir("hit_network"),
  "moderated_network_no_pcmeA_unique_top25_compounds.tsv"
)
highlight_compounds <- character()
if (file.exists(highlight_file)) {
  highlight_compounds <- read_tsv_base(highlight_file)$node_id
}

clean_label <- function(x, max_chars = 34) {
  out <- iconv(x, from = "", to = "ASCII//TRANSLIT")
  out[is.na(out)] <- x[is.na(out)]
  out <- gsub("\u03b3", "gamma", out, fixed = TRUE)
  too_long <- nchar(out) > max_chars
  out[too_long] <- paste0(substr(out[too_long], 1, max_chars - 3), "...")
  out
}

label_pairs <- function(tab, method, max_labels = 32) {
  padj <- padj_column(method, adjustment)
  hit <- paste0(method, "_hit")
  d <- tab[tab[[hit]] & is.finite(tab[[padj]]), , drop = FALSE]
  if (nrow(d) == 0) {
    return(d[FALSE, , drop = FALSE])
  }

  priority <- d[d$compound %in% highlight_compounds, , drop = FALSE]
  if (nrow(priority) > 0) {
    priority <- do.call(rbind, lapply(split(priority, priority$compound), function(x) {
      x[order(x[[padj]], -abs(x[[paste0(method, "_effect")]])), , drop = FALSE][1, , drop = FALSE]
    }))
  }

  most_sig <- d[order(d[[padj]], -abs(d[[paste0(method, "_effect")]])), , drop = FALSE]
  labeled <- rbind(priority, most_sig)
  labeled <- labeled[!duplicated(labeled$pair_id), , drop = FALSE]
  labeled <- head(labeled, max_labels)
  labeled$.label <- paste0(labeled$promoter, ": ", clean_label(labeled$compound_label, 28))
  labeled
}

volcano_plot <- function(tab, method) {
  effect <- paste0(method, "_effect")
  pvalue <- paste0(method, "_pvalue")
  padj <- padj_column(method, adjustment)
  hit <- paste0(method, "_hit")
  direction <- paste0(method, "_direction")

  d <- tab
  d$.effect <- d[[effect]]
  d$.neglog10p <- safe_neglog10(d[[pvalue]])
  d$.class <- "Not significant"
  d$.class[d[[hit]] & d[[direction]] == "Positive"] <- "Positive hit"
  d$.class[d[[hit]] & d[[direction]] == "Negative"] <- "Negative hit"
  d$.class <- factor(d$.class, levels = c("Negative hit", "Not significant", "Positive hit"))

  labels <- label_pairs(d, method)
  labels$.effect <- labels[[effect]]
  labels$.neglog10p <- safe_neglog10(labels[[pvalue]])

  ggplot(d, aes(.effect, .neglog10p)) +
    geom_point(aes(color = .class), alpha = 0.28, size = 0.45) +
    geom_text_repel(
      data = labels,
      aes(label = .label),
      size = 2.1,
      max.overlaps = Inf,
      box.padding = 0.22,
      point.padding = 0.08,
      min.segment.length = 0.01,
      segment.size = 0.18,
      seed = 17
    ) +
    scale_color_manual(
      values = c("Negative hit" = "#2563eb", "Not significant" = "#94a3b8", "Positive hit" = "#b91c1c"),
      name = NULL
    ) +
    theme_bw(base_size = 9) +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    ) +
    labs(
      title = method_labels[[method]],
      subtitle = paste0("Global BH FDR < 0.05; labeled pairs prioritize highlighted compounds"),
      x = "Effect size",
      y = expression(-log[10](p))
    )
}

circle_df <- function(cx, cy, r, n = 240) {
  theta <- seq(0, 2 * pi, length.out = n)
  data.frame(x = cx + r * cos(theta), y = cy + r * sin(theta))
}

venn_counts <- function(stats_set, median_set) {
  c(
    stats_only = length(setdiff(stats_set, median_set)),
    overlap = length(intersect(stats_set, median_set)),
    median_only = length(setdiff(median_set, stats_set))
  )
}

venn_plot <- function(tab, direction_filter = NULL, title = "All hits") {
  sets <- lapply(methods, function(method) {
    hit <- paste0(method, "_hit")
    direction <- paste0(method, "_direction")
    keep <- tab[[hit]]
    if (!is.null(direction_filter)) {
      keep <- keep & tab[[direction]] == direction_filter
    }
    tab$pair_id[keep]
  })
  names(sets) <- methods
  counts <- venn_counts(sets$destress_standard, sets$median_polish)
  left <- circle_df(-0.55, 0, 1)
  right <- circle_df(0.55, 0, 1)

  ggplot() +
    geom_polygon(data = left, aes(x, y), fill = "#60a5fa", alpha = 0.35, color = "#1d4ed8", linewidth = 0.5) +
    geom_polygon(data = right, aes(x, y), fill = "#f87171", alpha = 0.35, color = "#b91c1c", linewidth = 0.5) +
    annotate("text", x = -1.0, y = 0, label = counts[["stats_only"]], size = 5, fontface = "bold") +
    annotate("text", x = 0, y = 0, label = counts[["overlap"]], size = 5, fontface = "bold") +
    annotate("text", x = 1.0, y = 0, label = counts[["median_only"]], size = 5, fontface = "bold") +
    annotate("text", x = -0.75, y = -1.18, label = "DStressR standard model", size = 3.1) +
    annotate("text", x = 0.75, y = -1.18, label = "Median-polish max-p model", size = 3.1) +
    coord_equal(xlim = c(-1.8, 1.8), ylim = c(-1.35, 1.2), clip = "off") +
    theme_void(base_size = 9) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    labs(title = title)
}

count_promoter_hits <- function(tab, method) {
  hit <- paste0(method, "_hit")
  direction <- paste0(method, "_direction")
  d <- tab[tab[[hit]], c("promoter", direction), drop = FALSE]
  names(d)[2] <- "direction"
  if (nrow(d) == 0) {
    return(data.frame())
  }
  out <- as.data.frame.matrix(table(d$promoter, d$direction), stringsAsFactors = FALSE)
  out$promoter <- rownames(out)
  if (!"Positive" %in% names(out)) out$Positive <- 0
  if (!"Negative" %in% names(out)) out$Negative <- 0
  out <- out[, c("promoter", "Positive", "Negative"), drop = FALSE]
  names(out)[2:3] <- paste(method, c("positive_hits", "negative_hits"), sep = "_")
  rownames(out) <- NULL
  out
}

promoters <- data.frame(promoter = sort(unique(pair_table$promoter)), stringsAsFactors = FALSE)
hit_counts <- Reduce(
  function(x, y) merge(x, y, by = "promoter", all = TRUE, sort = FALSE),
  c(list(promoters), lapply(methods, function(method) count_promoter_hits(pair_table, method)))
)
hit_counts[is.na(hit_counts)] <- 0
hit_counts <- hit_counts[order(
  -hit_counts$destress_standard_positive_hits - hit_counts$destress_standard_negative_hits,
  hit_counts$promoter
), , drop = FALSE]

write.table(
  hit_counts,
  file.path(out_dir, "promoter_positive_negative_hit_counts.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

plot_count_table <- function(hit_counts, method, stem) {
  color_cap <- as.numeric(Sys.getenv("DSTRESSR_HIT_TABLE_COLOR_CAP", unset = "20"))
  if (!is.finite(color_cap) || color_cap <= 0) {
    color_cap <- 20
  }
  positive_col <- paste0(method, "_positive_hits")
  negative_col <- paste0(method, "_negative_hits")
  d <- hit_counts[, c("promoter", positive_col, negative_col), drop = FALSE]
  names(d) <- c("promoter", "Up-regulated", "Down-regulated")
  d$total <- d[["Up-regulated"]] + d[["Down-regulated"]]
  d <- d[order(-d$total, -d[["Up-regulated"]], -d[["Down-regulated"]], d$promoter), , drop = FALSE]

  long <- rbind(
    data.frame(promoter = d$promoter, direction = "Down-regulated", count = d[["Down-regulated"]]),
    data.frame(promoter = d$promoter, direction = "Up-regulated", count = d[["Up-regulated"]])
  )
  long$promoter <- factor(long$promoter, levels = rev(d$promoter))
  long$direction <- factor(long$direction, levels = c("Down-regulated", "Up-regulated"))
  long$fill_score <- log1p(pmin(long$count, color_cap)) / log1p(color_cap)
  long$text_color <- ifelse(long$fill_score >= 0.68, "white", "#111827")

  p <- ggplot(long, aes(direction, promoter)) +
    geom_tile(aes(fill = fill_score), color = "#52525b", linewidth = 0.28, width = 0.98, height = 0.98) +
    geom_text(aes(label = count, color = text_color), size = 2.35) +
    scale_fill_gradient(
      low = "#f8fafc",
      high = "#b91c1c",
      limits = c(0, 1),
      guide = "none"
    ) +
    scale_color_identity() +
    coord_equal(clip = "off") +
    theme_void(base_size = 8) +
    theme(
      plot.background = element_rect(fill = "black", color = NA),
      panel.background = element_rect(fill = "black", color = NA),
      plot.title = element_text(color = "white", face = "bold", size = 8.2, hjust = 0.5, lineheight = 0.95),
      axis.text.y = element_text(color = "white", size = 6.8, hjust = 1),
      axis.text.x = element_text(color = "white", size = 6.6, angle = 45, hjust = 1, vjust = 1),
      axis.ticks = element_blank(),
      plot.margin = margin(6, 7, 14, 6)
    ) +
    labs(title = table_labels[[method]], x = NULL, y = NULL)

  ggsave(file.path(out_dir, paste0(stem, ".png")), p, width = 2.35, height = 6.2, dpi = 300, bg = "black")
  ggsave(file.path(out_dir, paste0(stem, ".pdf")), p, width = 2.35, height = 6.2, bg = "black")
  write.table(
    d[, c("promoter", "Down-regulated", "Up-regulated"), drop = FALSE],
    file.path(out_dir, paste0(stem, ".tsv")),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  p
}

method_promoter_order <- function(hit_counts, method) {
  positive_col <- paste0(method, "_positive_hits")
  negative_col <- paste0(method, "_negative_hits")
  d <- hit_counts[, c("promoter", positive_col, negative_col), drop = FALSE]
  names(d) <- c("promoter", "positive", "negative")
  d$total <- d$positive + d$negative
  d <- d[order(-d$total, -d$positive, -d$negative, d$promoter), , drop = FALSE]
  d$promoter
}

pvalue_histogram_plot <- function(tab, method, promoter_order) {
  pvalue <- paste0(method, "_pvalue")
  hit <- paste0(method, "_hit")
  direction <- paste0(method, "_direction")
  d <- data.frame(
    promoter = tab$promoter,
    pvalue = tab[[pvalue]],
    hit = tab[[hit]],
    direction = tab[[direction]],
    stringsAsFactors = FALSE
  )
  d <- d[is.finite(d$pvalue), , drop = FALSE]
  d$promoter <- factor(d$promoter, levels = promoter_order)
  hits <- d[d$hit, , drop = FALSE]
  hits$direction <- factor(hits$direction, levels = c("Negative", "Positive"))

  ggplot(d, aes(pvalue)) +
    geom_histogram(binwidth = 0.05, boundary = 0, fill = "#d1d5db", color = "white", linewidth = 0.08) +
    geom_histogram(
      data = hits,
      aes(fill = direction),
      binwidth = 0.05,
      boundary = 0,
      color = "white",
      linewidth = 0.08,
      alpha = 0.9
    ) +
    facet_wrap(~ promoter, ncol = 6, drop = FALSE) +
    scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), expand = c(0.01, 0.01)) +
    scale_fill_manual(values = c("Negative" = "#2563eb", "Positive" = "#b91c1c"), name = "Significant hit") +
    theme_bw(base_size = 7) +
    theme(
      strip.background = element_rect(fill = "#e5e7eb", color = "#52525b", linewidth = 0.25),
      strip.text = element_text(face = "bold", size = 5.2),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "#e5e7eb", linewidth = 0.18),
      axis.text = element_text(size = 5),
      axis.title = element_text(size = 7),
      legend.position = "bottom",
      legend.title = element_text(size = 6.5),
      legend.text = element_text(size = 6.2),
      plot.title = element_text(face = "bold", size = 10),
      plot.subtitle = element_text(size = 8),
      plot.margin = margin(5, 4, 5, 5)
    ) +
    labs(
      title = method_labels[[method]],
      subtitle = "Raw p-value distributions by promoter; significant hits overlaid by direction",
      x = "p-value",
      y = "Count"
    )
}

write_histogram_table_combined <- function(tab, hit_counts, method, table_plot, stem) {
  promoter_order <- method_promoter_order(hit_counts, method)
  hist_plot <- pvalue_histogram_plot(tab, method, promoter_order)
  compact_table <- table_plot + theme(plot.title = element_blank())
  combined <- arrangeGrob(hist_plot, compact_table, ncol = 2, widths = c(3.1, 1.15))
  png(file.path(out_dir, paste0(stem, ".png")), width = 14, height = 6.4, units = "in", res = 300)
  grid.draw(combined)
  dev.off()
  pdf(file.path(out_dir, paste0(stem, ".pdf")), width = 14, height = 6.4)
  grid.draw(combined)
  dev.off()
  invisible(combined)
}

volcano_stats <- volcano_plot(pair_table, "destress_standard")
volcano_median <- volcano_plot(pair_table, "median_polish")
venn_all <- venn_plot(pair_table, NULL, "All significant pairs")
venn_pos <- venn_plot(pair_table, "Positive", "Positive hits")
venn_neg <- venn_plot(pair_table, "Negative", "Negative hits")

volcano_grob <- arrangeGrob(volcano_stats, volcano_median, ncol = 2)
venn_grob <- arrangeGrob(venn_all, venn_pos, venn_neg, ncol = 3)
combined <- arrangeGrob(
  volcano_grob,
  venn_grob,
  ncol = 1,
  heights = c(2.2, 1)
)

png(file.path(out_dir, "volcano_venn_significance_summary.png"), width = 16, height = 11, units = "in", res = 300)
grid.draw(combined)
dev.off()
pdf(file.path(out_dir, "volcano_venn_significance_summary.pdf"), width = 16, height = 11)
grid.draw(combined)
dev.off()

dstressr_table <- plot_count_table(hit_counts, "destress_standard", "dstressr_standard_model_promoter_hit_count_table")
dstressr_table_compat <- plot_count_table(hit_counts, "destress_standard", "dstressr_model_promoter_hit_count_table")
median_table <- plot_count_table(hit_counts, "median_polish", "median_polish_model_promoter_hit_count_table")
write_histogram_table_combined(
  pair_table,
  hit_counts,
  "destress_standard",
  dstressr_table,
  "dstressr_standard_model_pvalue_histograms_with_hit_table"
)
write_histogram_table_combined(
  pair_table,
  hit_counts,
  "destress_standard",
  dstressr_table_compat,
  "dstressr_model_pvalue_histograms_with_hit_table"
)
write_histogram_table_combined(
  pair_table,
  hit_counts,
  "median_polish",
  median_table,
  "median_polish_model_pvalue_histograms_with_hit_table"
)

method_summary <- data.frame(
  method = unname(method_labels[methods]),
  total_hits = vapply(methods, function(method) sum(pair_table[[paste0(method, "_hit")]]), integer(1)),
  positive_hits = vapply(methods, function(method) {
    sum(pair_table[[paste0(method, "_hit")]] & pair_table[[paste0(method, "_direction")]] == "Positive")
  }, integer(1)),
  negative_hits = vapply(methods, function(method) {
    sum(pair_table[[paste0(method, "_hit")]] & pair_table[[paste0(method, "_direction")]] == "Negative")
  }, integer(1)),
  stringsAsFactors = FALSE
)
write.table(
  method_summary,
  file.path(out_dir, "method_hit_summary.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

message("Wrote significance summary to: ", out_dir)
print(method_summary)
