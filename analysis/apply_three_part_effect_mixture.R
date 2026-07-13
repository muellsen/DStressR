#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

if (requireNamespace("DStressR", quietly = TRUE)) {
  suppressPackageStartupMessages(library(DStressR))
} else {
  suppressPackageStartupMessages(devtools::load_all(".", quiet = TRUE))
}

out_dir <- file.path(getwd(), "analysis", "outputs", "three_part_mixture")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

input_file <- file.path(
  getwd(),
  "analysis",
  "outputs",
  "normalized_matrix",
  "normalized_promoter_compound_matrix_long.tsv"
)
tab <- read.delim(input_file, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)

mix <- fit_effect_mixture(
  tab,
  value = "destress_eb_effect_centered",
  promoter = "promoter",
  df = 4,
  max_iter = 2000
)
mix_summary <- attr(mix, "mixture_summary")
mix$three_part_hit <- is.finite(mix$local_fdr) &
  mix$local_fdr <= 0.20 &
  mix$posterior_class %in% c("repressed", "activated")
mix$three_part_q_hit <- is.finite(mix$local_fdr_qvalue_by_promoter) &
  mix$local_fdr_qvalue_by_promoter <= 0.20 &
  mix$posterior_class %in% c("repressed", "activated")
mix$three_part_q05_hit <- is.finite(mix$local_fdr_qvalue_by_promoter) &
  mix$local_fdr_qvalue_by_promoter <= 0.05 &
  mix$posterior_class %in% c("repressed", "activated")
mix$three_part_bh_hit <- is.finite(mix$empirical_null_padj_by_promoter) &
  mix$empirical_null_padj_by_promoter < 0.05

write.table(
  mix,
  file.path(out_dir, "promoter_compound_three_part_mixture_results.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
write.table(
  mix_summary,
  file.path(out_dir, "promoter_three_part_mixture_parameters.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
write.table(
  mix[mix$three_part_q05_hit, , drop = FALSE],
  file.path(out_dir, "promoter_compound_three_part_local_fdr_q05_hits.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
write.table(
  mix[mix$three_part_q_hit, , drop = FALSE],
  file.path(out_dir, "promoter_compound_three_part_local_fdr_q20_hits.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

class_counts <- as.data.frame(
  xtabs(~ promoter + posterior_class, data = mix),
  stringsAsFactors = FALSE
)
names(class_counts) <- c("promoter", "posterior_class", "n")
write.table(
  class_counts,
  file.path(out_dir, "promoter_three_part_class_counts.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

summary_counts <- data.frame(
  quantity = c(
    "rows",
    "finite_local_fdr",
    "local_fdr_le_0.20_nonnull",
    "local_fdr_q_le_0.05_nonnull",
    "local_fdr_q_le_0.20_nonnull",
    "empirical_null_bh_lt_0.05",
    "eb_bh_lt_0.05",
    "both_local_fdr_q_and_eb_bh",
    "local_fdr_q_only",
    "eb_bh_only"
  ),
  value = c(
    nrow(mix),
    sum(is.finite(mix$local_fdr)),
    sum(mix$three_part_hit, na.rm = TRUE),
    sum(mix$three_part_q05_hit, na.rm = TRUE),
    sum(mix$three_part_q_hit, na.rm = TRUE),
    sum(mix$three_part_bh_hit, na.rm = TRUE),
    sum(mix$estimated_alpha_eb_padj_by_promoter < 0.05, na.rm = TRUE),
    sum(mix$three_part_q_hit & mix$estimated_alpha_eb_padj_by_promoter < 0.05, na.rm = TRUE),
    sum(mix$three_part_q_hit & !(mix$estimated_alpha_eb_padj_by_promoter < 0.05), na.rm = TRUE),
    sum(!mix$three_part_q_hit & mix$estimated_alpha_eb_padj_by_promoter < 0.05, na.rm = TRUE)
  )
)
write.table(
  summary_counts,
  file.path(out_dir, "three_part_mixture_summary_counts.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

selected_promoters <- mix_summary[order(-(mix_summary$prior_repressed + mix_summary$prior_activated)), "promoter"]
selected_promoters <- utils::head(selected_promoters, 12)
plot_df <- mix[mix$promoter %in% selected_promoters, , drop = FALSE]
density_grid <- do.call(
  rbind,
  lapply(selected_promoters, function(prom) {
    d <- plot_df[plot_df$promoter == prom, , drop = FALSE]
    s <- mix_summary[mix_summary$promoter == prom, , drop = FALSE]
    lim <- stats::quantile(abs(d$destress_eb_effect_centered), 0.995, na.rm = TRUE)
    if (!is.finite(lim) || lim <= 0) lim <- max(abs(d$destress_eb_effect_centered), na.rm = TRUE)
    x <- seq(-lim, lim, length.out = 250)
    dens_t <- function(x, location, scale, df) stats::dt((x - location) / scale, df = df) / scale
    data.frame(
      promoter = prom,
      effect = rep(x, 4),
      density = c(
        s$prior_repressed * dens_t(x, s$location_repressed, s$scale_repressed, s$df),
        s$prior_null * dens_t(x, s$location_null, s$scale_null, s$df),
        s$prior_activated * dens_t(x, s$location_activated, s$scale_activated, s$df),
        s$prior_repressed * dens_t(x, s$location_repressed, s$scale_repressed, s$df) +
          s$prior_null * dens_t(x, s$location_null, s$scale_null, s$df) +
          s$prior_activated * dens_t(x, s$location_activated, s$scale_activated, s$df)
      ),
      component = rep(c("repressed", "null", "activated", "mixture"), each = length(x)),
      stringsAsFactors = FALSE
    )
  })
)

p_fit <- ggplot(plot_df, aes(destress_eb_effect_centered)) +
  geom_histogram(aes(y = after_stat(density)), bins = 55, fill = "#D8DEE9", color = "white", linewidth = 0.15) +
  geom_line(
    data = density_grid,
    aes(effect, density, color = component, linewidth = component),
    inherit.aes = FALSE
  ) +
  scale_color_manual(values = c(
    repressed = "#2166AC",
    null = "#303030",
    activated = "#B2182B",
    mixture = "#009E73"
  )) +
  scale_linewidth_manual(values = c(repressed = 0.45, null = 0.55, activated = 0.45, mixture = 0.75), guide = "none") +
  facet_wrap(vars(promoter), scales = "free_y", ncol = 4) +
  theme_light(base_size = 9) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    plot.title.position = "plot"
  ) +
  labs(
    title = "Three-part Student-t mixture fits",
    subtitle = "Selected promoters with largest estimated non-null mixture mass",
    x = "Centered DStressR EB effect",
    y = "Density",
    color = "Component"
  )
ggsave(file.path(out_dir, "three_part_mixture_selected_promoter_fits.png"), p_fit, width = 12, height = 8, dpi = 300, bg = "white")
ggsave(file.path(out_dir, "three_part_mixture_selected_promoter_fits.pdf"), p_fit, width = 12, height = 8, bg = "white")

p_lfdr <- ggplot(mix, aes(destress_eb_effect_centered, local_fdr)) +
  geom_point(aes(color = posterior_class), alpha = 0.35, size = 0.6) +
  geom_hline(yintercept = 0.20, linetype = "longdash", color = "#303030", linewidth = 0.35) +
  scale_color_manual(values = c(activated = "#B2182B", null = "#9AA0A6", repressed = "#2166AC")) +
  theme_light(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom") +
  labs(
    title = "Three-part mixture local FDR",
    subtitle = "Dashed line marks local FDR = 0.20",
    x = "Centered DStressR EB effect",
    y = "Posterior null probability / local FDR",
    color = "Posterior class"
  )
ggsave(file.path(out_dir, "three_part_mixture_local_fdr_vs_effect.png"), p_lfdr, width = 8, height = 5.2, dpi = 300, bg = "white")
ggsave(file.path(out_dir, "three_part_mixture_local_fdr_vs_effect.pdf"), p_lfdr, width = 8, height = 5.2, bg = "white")

message("Wrote three-part mixture outputs to: ", out_dir)
print(summary_counts, row.names = FALSE)
