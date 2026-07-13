#!/usr/bin/env Rscript

source(file.path("analysis", "_helpers.R"))

suppressPackageStartupMessages({
  library(ggplot2)
})

if (requireNamespace("DStressR", quietly = TRUE)) {
  suppressPackageStartupMessages(library(DStressR))
} else {
  load_destress_package()
}

out_dir <- file.path(getwd(), "analysis", "outputs", "normalized_matrix", "effect_histograms")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

matrix_long_file <- file.path(
  getwd(),
  "analysis",
  "outputs",
  "normalized_matrix",
  "normalized_promoter_compound_matrix_long.tsv"
)
tab <- read.delim(matrix_long_file, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
tab$destress_eb_effect_centered <- as.numeric(tab$destress_eb_effect_centered)

effect_limit <- stats::quantile(abs(tab$destress_eb_effect_centered), 0.995, na.rm = TRUE)
if (!is.finite(effect_limit) || effect_limit <= 0) {
  effect_limit <- max(abs(tab$destress_eb_effect_centered), na.rm = TRUE)
}

p_all <- plot_effect_histogram(
  tab,
  value = "destress_eb_effect_centered",
  by = "all",
  bins = 120,
  xlim = c(-effect_limit, effect_limit),
  title = "Distribution of DStressR EB effects",
  subtitle = "All promoter-compound entries; x-axis clipped at the 99.5th percentile of absolute effects",
  xlab = "Centered DStressR EB effect"
)

p_promoter <- plot_effect_histogram(
  tab,
  value = "destress_eb_effect_centered",
  promoter = "promoter",
  by = "promoter",
  bins = 70,
  xlim = c(-effect_limit, effect_limit),
  title = "DStressR EB effect distributions per promoter",
  subtitle = "Each panel shows compounds for one promoter; x-axis shared and clipped at the pooled 99.5th percentile",
  xlab = "Centered DStressR EB effect"
)

ggsave(
  file.path(out_dir, "effect_histogram_all_entries.png"),
  p_all,
  width = 8.5,
  height = 5,
  dpi = 300,
  bg = "white"
)
ggsave(
  file.path(out_dir, "effect_histogram_all_entries.pdf"),
  p_all,
  width = 8.5,
  height = 5,
  bg = "white"
)
ggsave(
  file.path(out_dir, "effect_histograms_by_promoter.png"),
  p_promoter,
  width = 12,
  height = 10,
  dpi = 300,
  bg = "white"
)
ggsave(
  file.path(out_dir, "effect_histograms_by_promoter.pdf"),
  p_promoter,
  width = 12,
  height = 10,
  bg = "white"
)

summary_df <- data.frame(
  scope = c("all", sort(unique(tab$promoter))),
  n = c(
    sum(is.finite(tab$destress_eb_effect_centered)),
    as.integer(stats::aggregate(
      destress_eb_effect_centered ~ promoter,
      tab,
      function(x) sum(is.finite(x))
    )$destress_eb_effect_centered)
  ),
  mean = c(
    mean(tab$destress_eb_effect_centered, na.rm = TRUE),
    stats::aggregate(destress_eb_effect_centered ~ promoter, tab, mean, na.rm = TRUE)$destress_eb_effect_centered
  ),
  median = c(
    stats::median(tab$destress_eb_effect_centered, na.rm = TRUE),
    stats::aggregate(destress_eb_effect_centered ~ promoter, tab, stats::median, na.rm = TRUE)$destress_eb_effect_centered
  ),
  sd = c(
    stats::sd(tab$destress_eb_effect_centered, na.rm = TRUE),
    stats::aggregate(destress_eb_effect_centered ~ promoter, tab, stats::sd, na.rm = TRUE)$destress_eb_effect_centered
  ),
  mad = c(
    stats::mad(tab$destress_eb_effect_centered, na.rm = TRUE),
    stats::aggregate(destress_eb_effect_centered ~ promoter, tab, stats::mad, na.rm = TRUE)$destress_eb_effect_centered
  )
)
write.table(
  summary_df,
  file.path(out_dir, "effect_histogram_summary.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

message("Wrote effect histograms to: ", out_dir)
print(utils::head(summary_df, 10), row.names = FALSE)
