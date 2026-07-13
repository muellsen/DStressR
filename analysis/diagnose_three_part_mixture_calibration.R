#!/usr/bin/env Rscript

source(file.path("analysis", "_helpers.R"))

suppressPackageStartupMessages({
  library(ggplot2)
})

out_dir <- file.path(getwd(), "analysis", "outputs", "three_part_mixture", "calibration")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

mixture_file <- file.path(
  getwd(),
  "analysis",
  "outputs",
  "three_part_mixture",
  "promoter_three_part_mixture_parameters.tsv"
)
dmso_summary_file <- file.path(
  getwd(),
  "analysis",
  "outputs",
  "normalized_matrix",
  "effect_histograms",
  "effect_histogram_compounds_vs_dmso_summary.tsv"
)

mix <- read.delim(mixture_file, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
summary <- read.delim(dmso_summary_file, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
dmso <- summary[summary$source == "DMSO controls", c("promoter", "n", "median", "sd", "mad")]
names(dmso) <- c("promoter", "dmso_n", "dmso_median", "dmso_sd", "dmso_mad")
compound <- summary[summary$source == "Compounds", c("promoter", "n", "median", "sd", "mad")]
names(compound) <- c("promoter", "compound_n", "compound_median", "compound_sd", "compound_mad")

diagnostic <- merge(mix, dmso, by = "promoter", all.x = TRUE, sort = FALSE)
diagnostic <- merge(diagnostic, compound, by = "promoter", all.x = TRUE, sort = FALSE)
diagnostic$nonnull_prior <- diagnostic$prior_repressed + diagnostic$prior_activated
diagnostic$null_shift_vs_dmso <- diagnostic$location_null - diagnostic$dmso_median
diagnostic$abs_null_shift_vs_dmso <- abs(diagnostic$null_shift_vs_dmso)
diagnostic$null_shift_in_dmso_mad <- diagnostic$abs_null_shift_vs_dmso / pmax(diagnostic$dmso_mad, .Machine$double.eps)
diagnostic$null_scale_vs_dmso_mad <- diagnostic$scale_null / pmax(diagnostic$dmso_mad, .Machine$double.eps)
diagnostic$flag_null_drift <- diagnostic$null_shift_in_dmso_mad > 1
diagnostic$flag_large_nonnull_prior <- diagnostic$nonnull_prior > 0.25
diagnostic$flag_unreliable_lfdr <- diagnostic$flag_null_drift | diagnostic$flag_large_nonnull_prior

write.table(
  diagnostic,
  file.path(out_dir, "three_part_mixture_null_calibration_diagnostic.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

p_shift <- ggplot(
  diagnostic,
  aes(dmso_median, location_null)
) +
  geom_abline(slope = 1, intercept = 0, linetype = "longdash", color = "#303030", linewidth = 0.35) +
  geom_point(aes(size = nonnull_prior, color = flag_unreliable_lfdr), alpha = 0.8) +
  geom_text(
    data = diagnostic[diagnostic$flag_unreliable_lfdr, , drop = FALSE],
    aes(label = promoter),
    size = 2.6,
    vjust = -0.7,
    check_overlap = TRUE
  ) +
  scale_color_manual(values = c("FALSE" = "#4E79A7", "TRUE" = "#B2182B")) +
  theme_light(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom", plot.title.position = "plot") +
  labs(
    title = "Three-part mixture null is not always anchored to DMSO",
    subtitle = "Point size is fitted non-null prior mass; red flags null drift or implausibly high non-null mass",
    x = "DMSO median effect",
    y = "Fitted mixture null location",
    color = "Flagged",
    size = "Fitted non-null prior"
  )
ggsave(file.path(out_dir, "three_part_null_location_vs_dmso.png"), p_shift, width = 7.2, height = 5.8, dpi = 300, bg = "white")
ggsave(file.path(out_dir, "three_part_null_location_vs_dmso.pdf"), p_shift, width = 7.2, height = 5.8, bg = "white")

p_prior <- ggplot(diagnostic, aes(reorder(promoter, nonnull_prior), nonnull_prior)) +
  geom_col(aes(fill = flag_unreliable_lfdr), width = 0.8) +
  geom_hline(yintercept = 0.25, linetype = "longdash", color = "#303030", linewidth = 0.35) +
  coord_flip() +
  scale_fill_manual(values = c("FALSE" = "#4E79A7", "TRUE" = "#B2182B"), guide = "none") +
  theme_light(base_size = 10) +
  theme(panel.grid.minor = element_blank(), plot.title.position = "plot") +
  labs(
    title = "Fitted three-part non-null mass by promoter",
    subtitle = "Large side-component mass suggests histogram partitioning rather than calibrated hit calling",
    x = "Promoter",
    y = "Prior(repressed) + Prior(activated)"
  )
ggsave(file.path(out_dir, "three_part_nonnull_prior_by_promoter.png"), p_prior, width = 7.2, height = 6.8, dpi = 300, bg = "white")
ggsave(file.path(out_dir, "three_part_nonnull_prior_by_promoter.pdf"), p_prior, width = 7.2, height = 6.8, bg = "white")

summary_counts <- data.frame(
  quantity = c(
    "promoters",
    "flag_null_drift_gt_1_dmso_mad",
    "flag_nonnull_prior_gt_0.25",
    "flag_unreliable_lfdr",
    "median_nonnull_prior",
    "max_nonnull_prior"
  ),
  value = c(
    nrow(diagnostic),
    sum(diagnostic$flag_null_drift, na.rm = TRUE),
    sum(diagnostic$flag_large_nonnull_prior, na.rm = TRUE),
    sum(diagnostic$flag_unreliable_lfdr, na.rm = TRUE),
    stats::median(diagnostic$nonnull_prior, na.rm = TRUE),
    max(diagnostic$nonnull_prior, na.rm = TRUE)
  )
)
write.table(
  summary_counts,
  file.path(out_dir, "three_part_mixture_calibration_summary.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

message("Wrote three-part mixture calibration diagnostics to: ", out_dir)
print(summary_counts, row.names = FALSE)
