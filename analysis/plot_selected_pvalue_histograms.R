#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

out_dir <- file.path(getwd(), "analysis", "outputs", "pvalue_histograms")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

comparison_file <- file.path(
  getwd(),
  "analysis",
  "outputs",
  "workflow_vs_destress_replicate_pvalues.tsv"
)

if (!file.exists(comparison_file)) {
  stop(
    "Missing comparison table: ", comparison_file,
    "\nRun analysis/compare_workflow_pvalues.R first.",
    call. = FALSE
  )
}

comparison <- read.delim(comparison_file, sep = "\t", check.names = FALSE)
comparison <- comparison[
  is.finite(comparison$pvalue) &
    is.finite(comparison$destress_pvalue),
]
if ("ProductName" %in% names(comparison)) {
  comparison <- comparison[!(comparison$ProductName %in% c("DMSO", "DMSO noisy")), ]
}

group_stats <- do.call(
  rbind,
  lapply(split(comparison, comparison$promoter_libplate_replicate), function(d) {
    data.frame(
      promoter_libplate_replicate = d$promoter_libplate_replicate[1],
      promoter = d$promoter[1],
      libplate = d$libplate[1],
      replicate = d$replicate[1],
      n = nrow(d),
      workflow_p_lt_005 = mean(d$pvalue < 0.05, na.rm = TRUE),
      destress_p_lt_005 = mean(d$destress_pvalue < 0.05, na.rm = TRUE),
      ks_distance = as.numeric(
        suppressWarnings(stats::ks.test(d$pvalue, d$destress_pvalue)$statistic)
      )
    )
  })
)
group_stats$delta_small_p <- group_stats$destress_p_lt_005 - group_stats$workflow_p_lt_005
group_stats$abs_delta_small_p <- abs(group_stats$delta_small_p)
group_stats$n <- as.numeric(group_stats$n)
group_stats$abs_delta_small_p <- as.numeric(group_stats$abs_delta_small_p)
group_stats <- group_stats[order(-group_stats$abs_delta_small_p, -group_stats$n), ]

selected <- head(group_stats$promoter_libplate_replicate, 10)

plot_df <- comparison[comparison$promoter_libplate_replicate %in% selected, ]
plot_df$example <- factor(
  plot_df$promoter_libplate_replicate,
  levels = selected
)

long_df <- rbind(
  data.frame(
    example = plot_df$example,
    promoter_libplate_replicate = plot_df$promoter_libplate_replicate,
    pvalue = plot_df$pvalue,
    method = "Median-polish"
  ),
  data.frame(
    example = plot_df$example,
    promoter_libplate_replicate = plot_df$promoter_libplate_replicate,
    pvalue = plot_df$destress_pvalue,
    method = "DStressR"
  )
)
long_df$method <- factor(long_df$method, levels = c("Median-polish", "DStressR"))

selected_stats <- group_stats[group_stats$promoter_libplate_replicate %in% selected, ]
selected_stats <- selected_stats[match(selected, selected_stats$promoter_libplate_replicate), ]
write.table(
  selected_stats,
  file.path(out_dir, "selected_pvalue_histogram_examples.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

p <- ggplot(long_df, aes(pvalue)) +
  geom_histogram(
    breaks = seq(0, 1, by = 0.025),
    fill = "#334155",
    color = "white",
    linewidth = 0.15
  ) +
  geom_vline(xintercept = 0.05, color = "#b91c1c", linewidth = 0.35) +
  facet_grid(example ~ method, scales = "free_y") +
  coord_cartesian(xlim = c(0, 1)) +
  theme_bw(base_size = 9) +
  theme(
    strip.text.y = element_text(angle = 0, hjust = 0),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "P-value histograms for selected promoter-library-replicate examples",
    subtitle = "Left: original median-polish workflow; right: DStressR compound-global-adjusted residual test",
    x = "p-value",
    y = "Count"
  )

ggsave(
  file.path(out_dir, "selected_pvalue_histograms_medianpolish_vs_destress.png"),
  p,
  width = 9,
  height = 12,
  dpi = 300
)
ggsave(
  file.path(out_dir, "selected_pvalue_histograms_medianpolish_vs_destress.pdf"),
  p,
  width = 9,
  height = 12
)

message("Wrote selected p-value histograms to: ", out_dir)
print(selected_stats)
