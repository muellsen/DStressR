#!/usr/bin/env Rscript

source(file.path("analysis", "_helpers.R"))

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package `ggplot2` is required for p-value histograms.", call. = FALSE)
}
if (!requireNamespace("gridExtra", quietly = TRUE)) {
  stop("Package `gridExtra` is required for p-value histograms.", call. = FALSE)
}

ggplot2 <- asNamespace("ggplot2")
gridExtra <- asNamespace("gridExtra")

out_dir <- comparison_results_dir("background_rank")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_pair_file <- function(path, method) {
  tab <- read_tsv_base(path)
  required <- c("promoter", "compound", "effect", "pvalue", "padj_global")
  missing <- setdiff(required, names(tab))
  if (length(missing) > 0) {
    stop("Missing columns in ", path, ": ", paste(missing, collapse = ", "), call. = FALSE)
  }
  tab <- tab[, required, drop = FALSE]
  tab$method <- method
  tab$hit <- is.finite(tab$padj_global) & tab$padj_global < 0.05
  tab$direction <- ifelse(tab$effect < 0, "Down-regulated", "Up-regulated")
  tab
}

rank0 <- read_pair_file(
  package_results_dir("destress_moderated_rank0_pair_results.tsv"),
  "DStressR default k=0"
)
rank1 <- read_pair_file(
  package_results_dir("destress_moderated_rank1_pair_results.tsv"),
  "DStressR default k=1"
)
median_polish <- read_pair_file(
  package_results_dir("median_polish_pair_results.tsv"),
  "Median polish"
)

all_pvalues <- rbind(rank0, rank1, median_polish)
all_pvalues$method <- factor(
  all_pvalues$method,
  levels = c("DStressR default k=0", "DStressR default k=1", "Median polish")
)

promoter_order <- names(sort(table(all_pvalues$promoter), decreasing = TRUE))
all_pvalues$promoter <- factor(all_pvalues$promoter, levels = promoter_order)

summary <- do.call(rbind, lapply(split(all_pvalues, all_pvalues$method), function(d) {
  data.frame(
    method = as.character(d$method[1]),
    n_pairs = nrow(d),
    raw_p_lt_0_01 = sum(d$pvalue < 0.01, na.rm = TRUE),
    raw_p_lt_0_05 = sum(d$pvalue < 0.05, na.rm = TRUE),
    global_bh_hits = sum(d$hit, na.rm = TRUE),
    positive_hits = sum(d$hit & d$effect > 0, na.rm = TRUE),
    negative_hits = sum(d$hit & d$effect < 0, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))
write.table(
  summary,
  file.path(out_dir, "rank1_pvalue_histogram_summary.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

plot_one <- function(tab, method) {
  d <- tab[tab$method == method, , drop = FALSE]
  ggplot2$ggplot(d, ggplot2$aes(pvalue)) +
    ggplot2$geom_histogram(
      bins = 24,
      boundary = 0,
      fill = "#e5e7eb",
      color = "white",
      linewidth = 0.12
    ) +
    ggplot2$geom_histogram(
      data = d[d$hit, , drop = FALSE],
      ggplot2$aes(fill = direction),
      bins = 24,
      boundary = 0,
      color = "white",
      linewidth = 0.12
    ) +
    ggplot2$scale_fill_manual(
      values = c("Down-regulated" = "#2563eb", "Up-regulated" = "#b91c1c"),
      drop = FALSE
    ) +
    ggplot2$scale_x_continuous(limits = c(0, 1), expand = c(0.01, 0.01)) +
    ggplot2$facet_wrap(~ promoter, ncol = 6, drop = FALSE) +
    ggplot2$labs(
      title = method,
      subtitle = "Raw p-value histograms by promoter; global-BH significant hits overlaid by direction",
      x = "Raw p-value",
      y = "Count",
      fill = NULL
    ) +
    ggplot2$theme_minimal(base_size = 8) +
    ggplot2$theme(
      plot.title = ggplot2$element_text(face = "bold"),
      legend.position = "bottom",
      strip.text = ggplot2$element_text(size = 6.4),
      panel.grid.minor = ggplot2$element_blank()
    )
}

plots <- lapply(levels(all_pvalues$method), function(method) plot_one(all_pvalues, method))
combined <- gridExtra$grid.arrange(grobs = plots, ncol = 1)

ggplot2$ggsave(
  file.path(out_dir, "rank1_pvalue_histograms.png"),
  combined,
  width = 11,
  height = 14.5,
  dpi = 300,
  bg = "white"
)
ggplot2$ggsave(
  file.path(out_dir, "rank1_pvalue_histograms.pdf"),
  combined,
  width = 11,
  height = 14.5,
  bg = "white"
)

message("Wrote rank-1 p-value histograms to: ", out_dir)
