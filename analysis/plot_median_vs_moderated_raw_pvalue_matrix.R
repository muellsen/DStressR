#!/usr/bin/env Rscript

source(file.path("analysis", "_helpers.R"))

suppressPackageStartupMessages({
  library(ggplot2)
})

methods <- c("median_polish", "destress_moderated")
out_dir <- comparison_results_dir("promoter_pvalue_matrix")

tabs <- lapply(methods, function(method) {
  tab <- read_package_pair_results(method)
  data.frame(
    promoter = tab$promoter,
    compound = tab$compound,
    method = method_label(method),
    pvalue = tab[[paste0(method, "_pvalue")]],
    stringsAsFactors = FALSE
  )
})

plot_data <- do.call(rbind, tabs)
plot_data <- plot_data[is.finite(plot_data$pvalue), , drop = FALSE]
plot_data$method <- factor(plot_data$method, levels = vapply(methods, method_label, character(1)))
plot_data$promoter <- factor(plot_data$promoter, levels = sort(unique(plot_data$promoter)))

summary <- do.call(rbind, lapply(split(plot_data, list(plot_data$method, plot_data$promoter), drop = TRUE), function(x) {
  data.frame(
    method = as.character(x$method[1]),
    promoter = as.character(x$promoter[1]),
    n = nrow(x),
    p_lt_0.05 = sum(x$pvalue < 0.05, na.rm = TRUE),
    median_pvalue = stats::median(x$pvalue, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))

p <- ggplot(plot_data, aes(pvalue)) +
  geom_histogram(binwidth = 0.05, boundary = 0, fill = "#2563eb", color = "white", linewidth = 0.12) +
  facet_grid(promoter ~ method) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  theme_bw(base_size = 8) +
  theme(
    strip.background = element_rect(fill = "#f8fafc", color = "#cbd5e1"),
    strip.text.x = element_text(face = "bold", size = 9),
    strip.text.y = element_text(face = "bold", size = 7, angle = 0),
    panel.grid.minor = element_blank(),
    panel.spacing.y = unit(0.08, "lines"),
    panel.spacing.x = unit(0.35, "lines")
  ) +
  labs(
    title = "Raw p-value histograms by promoter",
    subtitle = "Median-polish max-p model vs DStressR moderated model",
    x = "Raw p-value",
    y = "Promoter-compound pairs"
  )

write.table(
  summary,
  file.path(out_dir, "median_polish_vs_destress_moderated_raw_pvalue_matrix_summary.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

ggsave(
  file.path(out_dir, "median_polish_vs_destress_moderated_raw_pvalue_matrix.png"),
  p,
  width = 8.5,
  height = 18,
  dpi = 300,
  bg = "white"
)
ggsave(
  file.path(out_dir, "median_polish_vs_destress_moderated_raw_pvalue_matrix.pdf"),
  p,
  width = 8.5,
  height = 18,
  bg = "white"
)

message("Wrote median-polish vs moderated p-value matrix to: ", out_dir)
print(summary[order(summary$promoter, summary$method), ])
