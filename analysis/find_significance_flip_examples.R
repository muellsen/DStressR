#!/usr/bin/env Rscript

source(file.path("analysis", "_helpers.R"))

suppressPackageStartupMessages({
  library(ggplot2)
})

read_tsv_base <- function(path) {
  read.delim(path, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
}

out_dir <- file.path(getwd(), "analysis", "outputs", "growth_exponent", "educational_examples")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fixed <- read_tsv_base(file.path(
  getwd(),
  "analysis",
  "outputs",
  "eb_moderated_variance",
  "workflow_vs_destress_eb_promoter_compound_pvalues.tsv"
))
adjusted <- read_tsv_base(file.path(
  getwd(),
  "analysis",
  "outputs",
  "growth_exponent",
  "workflow_vs_destress_eb_estimated_growth_alpha_promoter_compound_pvalues.tsv"
))

key <- paste(adjusted$promoter, adjusted$srn_code, sep = "\r")
fixed_key <- paste(fixed$promoter, fixed$srn_code, sep = "\r")
fixed_idx <- match(key, fixed_key)

comparison <- data.frame(
  promoter = adjusted$promoter,
  compound = adjusted$srn_code,
  replicate = adjusted$replicate,
  adjusted_alpha = adjusted$growth_alpha,
  median_polish_effect = adjusted$log2FC.polished,
  median_polish_p = adjusted$pvalue,
  fixed_effect = fixed$destress_eb_effect_centered[fixed_idx],
  fixed_p = fixed$destress_eb_pvalue[fixed_idx],
  fixed_bh = fixed$destress_eb_padj_by_promoter[fixed_idx],
  adjusted_effect = adjusted$destress_eb_effect_centered,
  adjusted_p = adjusted$destress_eb_pvalue,
  adjusted_bh = adjusted$estimated_alpha_eb_padj_by_promoter,
  stringsAsFactors = FALSE
)
comparison$fixed_sig <- comparison$fixed_bh < 0.05
comparison$adjusted_sig <- comparison$adjusted_bh < 0.05
comparison$effect_delta <- comparison$adjusted_effect - comparison$fixed_effect
comparison$abs_effect_delta <- abs(comparison$effect_delta)
comparison$change_type <- ifelse(
  comparison$fixed_sig & !comparison$adjusted_sig,
  "Fixed alpha only",
  ifelse(!comparison$fixed_sig & comparison$adjusted_sig, "Adjusted alpha only", "No flip")
)

flips <- comparison[comparison$change_type != "No flip" & is.finite(comparison$abs_effect_delta), ]
fixed_only <- flips[flips$change_type == "Fixed alpha only", ]
adjusted_only <- flips[flips$change_type == "Adjusted alpha only", ]
fixed_only <- fixed_only[order(-fixed_only$abs_effect_delta), ]
adjusted_only <- adjusted_only[order(-adjusted_only$abs_effect_delta), ]
examples <- rbind(head(fixed_only, 6), head(adjusted_only, 6))

examples <- examples[, c(
  "change_type",
  "promoter",
  "compound",
  "replicate",
  "adjusted_alpha",
  "median_polish_effect",
  "median_polish_p",
  "fixed_effect",
  "fixed_p",
  "fixed_bh",
  "adjusted_effect",
  "adjusted_p",
  "adjusted_bh",
  "effect_delta"
)]

write.table(
  examples,
  file.path(out_dir, "significance_flip_examples.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

plot_df <- rbind(
  data.frame(
    change_type = examples$change_type,
    promoter = examples$promoter,
    compound = examples$compound,
    method = "DStressR EB, alpha = 1",
    effect = examples$fixed_effect,
    pvalue = examples$fixed_p,
    bh = examples$fixed_bh,
    stringsAsFactors = FALSE
  ),
  data.frame(
    change_type = examples$change_type,
    promoter = examples$promoter,
    compound = examples$compound,
    method = "DStressR EB, adjusted alpha_g",
    effect = examples$adjusted_effect,
    pvalue = examples$adjusted_p,
    bh = examples$adjusted_bh,
    stringsAsFactors = FALSE
  )
)
plot_df$pair <- paste(plot_df$promoter, plot_df$compound, sep = " / ")
plot_df$pair <- factor(plot_df$pair, levels = unique(plot_df$pair))
plot_df$method <- factor(
  plot_df$method,
  levels = c("DStressR EB, alpha = 1", "DStressR EB, adjusted alpha_g")
)
plot_df$neglog10_bh <- -log10(pmax(plot_df$bh, .Machine$double.xmin))

p_effect <- ggplot(plot_df, aes(method, effect, fill = method)) +
  geom_col(width = 0.7, color = "white", linewidth = 0.2) +
  geom_hline(yintercept = 0, color = "#334155", linewidth = 0.25) +
  facet_grid(change_type ~ pair, scales = "free_y") +
  scale_fill_manual(values = c("#2563eb", "#dc2626"), guide = "none") +
  theme_bw(base_size = 8) +
  theme(
    axis.text.x = element_text(angle = 22, hjust = 1),
    panel.grid.minor = element_blank(),
    strip.text.x = element_text(size = 6)
  ) +
  labs(
    title = "Examples where adjusted growth normalization changes effect size and BH decision",
    x = NULL,
    y = "EB-centered effect"
  )

p_bh <- ggplot(plot_df, aes(method, neglog10_bh, fill = method)) +
  geom_col(width = 0.7, color = "white", linewidth = 0.2) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, color = "#334155", linewidth = 0.25) +
  facet_grid(change_type ~ pair, scales = "free_y") +
  scale_fill_manual(values = c("#2563eb", "#dc2626"), guide = "none") +
  theme_bw(base_size = 8) +
  theme(
    axis.text.x = element_text(angle = 22, hjust = 1),
    panel.grid.minor = element_blank(),
    strip.text.x = element_text(size = 6)
  ) +
  labs(
    title = "The same examples on the BH-adjusted p-value scale",
    subtitle = "Dashed line marks BH = 0.05",
    x = NULL,
    y = expression(-log[10](BH))
  )

ggsave(
  file.path(out_dir, "significance_flip_examples_effects.png"),
  p_effect,
  width = 13,
  height = 5,
  dpi = 300,
  bg = "white"
)
ggsave(
  file.path(out_dir, "significance_flip_examples_effects.pdf"),
  p_effect,
  width = 13,
  height = 5,
  bg = "white"
)
ggsave(
  file.path(out_dir, "significance_flip_examples_bh.png"),
  p_bh,
  width = 13,
  height = 5,
  dpi = 300,
  bg = "white"
)
ggsave(
  file.path(out_dir, "significance_flip_examples_bh.pdf"),
  p_bh,
  width = 13,
  height = 5,
  bg = "white"
)

message("Wrote significance-flip examples to: ", out_dir)
message("Fixed alpha only BH hits: ", sum(comparison$fixed_sig & !comparison$adjusted_sig, na.rm = TRUE))
message("Adjusted alpha only BH hits: ", sum(!comparison$fixed_sig & comparison$adjusted_sig, na.rm = TRUE))
print(examples, digits = 4, row.names = FALSE)
