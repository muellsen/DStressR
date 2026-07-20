#!/usr/bin/env Rscript

source(file.path("analysis", "_helpers.R"))

suppressPackageStartupMessages({
  library(ggplot2)
})

methods <- c("median_polish", "destress_standard", "destress_moderated")
out_dir <- comparison_results_dir("promoter_padj_histograms")

tabs <- lapply(methods, function(method) {
  tab <- read_package_pair_results(method)
  data.frame(
    promoter = tab$promoter,
    compound = tab$compound,
    method = method_label(method),
    pvalue = tab[[paste0(method, "_pvalue")]],
    padj_global = tab[[paste0(method, "_padj_global")]],
    padj_by_promoter = tab[[paste0(method, "_padj_by_promoter")]],
    stringsAsFactors = FALSE
  )
})
hist_data <- do.call(rbind, tabs)
hist_data$method <- factor(hist_data$method, levels = vapply(methods, method_label, character(1)))

plot_histograms <- function(label, x_column, x_label) {
  d <- hist_data[is.finite(hist_data[[x_column]]), , drop = FALSE]
  d$p <- d[[x_column]]
  d$promoter <- factor(d$promoter, levels = sort(unique(d$promoter)))

  summary <- do.call(rbind, lapply(split(d, list(d$method, d$promoter), drop = TRUE), function(x) {
    data.frame(
      method = as.character(x$method[1]),
      promoter = as.character(x$promoter[1]),
      n = nrow(x),
      p_lt_0.05 = sum(x$p < 0.05, na.rm = TRUE),
      median_p = stats::median(x$p, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))

  p <- ggplot(d, aes(p, fill = method)) +
    geom_histogram(binwidth = 0.05, boundary = 0, position = "identity",
                   alpha = 0.42, color = NA) +
    facet_wrap(~ promoter, ncol = 6) +
    scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
    scale_fill_manual(values = c(
      "Median-polish max-p model" = "#64748b",
      "DStressR ordinary model" = "#2563eb",
      "DStressR default (moderated)" = "#16a34a"
    )) +
    theme_bw(base_size = 9) +
    theme(
      legend.position = "bottom",
      strip.background = element_rect(fill = "#f8fafc", color = "#cbd5e1"),
      strip.text = element_text(face = "bold", size = 8),
      panel.grid.minor = element_blank()
    ) +
    labs(
      title = paste(x_label, "histograms by promoter"),
      x = x_label,
      y = "Promoter-compound pairs",
      fill = "Method"
    )

  prefix <- paste0("promoter_", label, "_pvalue_histograms")
  write.table(
    summary,
    file.path(out_dir, paste0(prefix, "_summary.tsv")),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  ggsave(file.path(out_dir, paste0(prefix, ".png")),
         p, width = 14, height = 10, dpi = 300, bg = "white")
  ggsave(file.path(out_dir, paste0(prefix, ".pdf")),
         p, width = 14, height = 10, bg = "white")
  invisible(summary)
}

raw_summary <- plot_histograms("raw", "pvalue", "Raw p-value")
global_summary <- plot_histograms("global_adjusted", "padj_global", "Global adjusted p-value")
by_promoter_summary <- plot_histograms("by_promoter_adjusted", "padj_by_promoter", "Within-promoter adjusted p-value")

message("Wrote promoter p-value histograms to: ", out_dir)
print(raw_summary[order(raw_summary$method, -raw_summary$p_lt_0.05), ])
print(global_summary[order(global_summary$method, -global_summary$p_lt_0.05), ])
print(by_promoter_summary[order(by_promoter_summary$method, -by_promoter_summary$p_lt_0.05), ])
