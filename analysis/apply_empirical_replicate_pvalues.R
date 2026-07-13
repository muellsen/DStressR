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

out_dir <- file.path(getwd(), "analysis", "outputs", "empirical_replicate_pvalues")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

replicate_file <- file.path(
  getwd(),
  "analysis",
  "outputs",
  "growth_exponent",
  "workflow_vs_destress_eb_estimated_growth_alpha_replicate_pvalues.tsv"
)
compound_file <- file.path(
  getwd(),
  "analysis",
  "outputs",
  "growth_exponent",
  "workflow_vs_destress_eb_estimated_growth_alpha_promoter_compound_pvalues.tsv"
)
libmap_file <- libmap_path()

replicate <- read.delim(replicate_file, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
compound_level <- read.delim(compound_file, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
libmap <- read.delim(libmap_file, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
libmap$libplate <- paste0("lp", libmap[["Library plate"]])
libmap$srn_code <- paste(libmap$libplate, libmap[["Well"]], sep = "_")
libmap$compound_label <- ifelse(
  is.na(libmap$ProductName) | libmap$ProductName == "NA" | !nzchar(libmap$ProductName),
  libmap[["Catalog Number"]],
  libmap$ProductName
)
dmso_srn_codes <- libmap$srn_code[libmap$compound_label %in% c("DMSO", "DMSO noisy")]

emp <- empirical_replicate_pvalues(
  replicate,
  value = "destress_eb_effect_centered",
  promoter = "promoter",
  compound = "srn_code",
  control = dmso_srn_codes,
  replicate = "replicate",
  strata = "libplate",
  min_replicates = 2,
  min_null = 5,
  permutation = TRUE,
  B = 1000,
  seed = 1,
  alternative = "two.sided"
)

emp <- merge(
  emp,
  libmap[, c("srn_code", "compound_label"), drop = FALSE],
  by = "srn_code",
  all.x = TRUE,
  sort = FALSE
)
emp <- merge(
  emp,
  compound_level[, c(
    "promoter",
    "srn_code",
    "estimated_alpha_eb_padj_by_promoter",
    "estimated_alpha_eb_neglog10padj",
    "destress_eb_pvalue"
  )],
  by = c("promoter", "srn_code"),
  all.x = TRUE,
  sort = FALSE
)
emp$empirical_neglog10p <- -log10(pmax(emp$empirical_pvalue, .Machine$double.xmin))
emp$empirical_neglog10padj <- -log10(pmax(emp$empirical_padj_by_promoter, .Machine$double.xmin))
emp$permutation_neglog10p <- -log10(pmax(emp$permutation_pvalue, .Machine$double.xmin))
emp$permutation_neglog10padj <- -log10(pmax(emp$permutation_padj_by_promoter, .Machine$double.xmin))
emp$empirical_bh_hit <- is.finite(emp$empirical_padj_by_promoter) & emp$empirical_padj_by_promoter < 0.05
emp$permutation_bh_hit <- is.finite(emp$permutation_padj_by_promoter) & emp$permutation_padj_by_promoter < 0.05
emp$eb_bh_hit <- is.finite(emp$estimated_alpha_eb_padj_by_promoter) &
  emp$estimated_alpha_eb_padj_by_promoter < 0.05

write.table(
  emp,
  file.path(out_dir, "promoter_compound_empirical_replicate_pvalues.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
write.table(
  emp[emp$empirical_bh_hit, , drop = FALSE],
  file.path(out_dir, "promoter_compound_empirical_replicate_bh_hits.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
write.table(
  emp[emp$permutation_bh_hit, , drop = FALSE],
  file.path(out_dir, "promoter_compound_permutation_replicate_bh_hits.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

summary_counts <- data.frame(
  quantity = c(
    "rows",
    "finite_empirical_p",
    "median_null_n",
    "min_empirical_p",
    "min_permutation_p",
    "empirical_raw_p_lt_0.05",
    "permutation_raw_p_lt_0.05",
    "empirical_bh_lt_0.05",
    "permutation_bh_lt_0.05",
    "eb_bh_lt_0.05",
    "both_permutation_and_eb_bh",
    "empirical_bh_only",
    "permutation_bh_only",
    "eb_bh_only"
  ),
  value = c(
    nrow(emp),
    sum(is.finite(emp$empirical_pvalue)),
    stats::median(emp$null_n, na.rm = TRUE),
    min(emp$empirical_pvalue, na.rm = TRUE),
    min(emp$permutation_pvalue, na.rm = TRUE),
    sum(emp$empirical_pvalue < 0.05, na.rm = TRUE),
    sum(emp$permutation_pvalue < 0.05, na.rm = TRUE),
    sum(emp$empirical_bh_hit, na.rm = TRUE),
    sum(emp$permutation_bh_hit, na.rm = TRUE),
    sum(emp$eb_bh_hit, na.rm = TRUE),
    sum(emp$permutation_bh_hit & emp$eb_bh_hit, na.rm = TRUE),
    sum(emp$empirical_bh_hit & !emp$eb_bh_hit, na.rm = TRUE),
    sum(emp$permutation_bh_hit & !emp$eb_bh_hit, na.rm = TRUE),
    sum(!emp$permutation_bh_hit & emp$eb_bh_hit, na.rm = TRUE)
  )
)
write.table(
  summary_counts,
  file.path(out_dir, "empirical_replicate_pvalue_summary.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

p_hist <- ggplot(emp, aes(empirical_pvalue)) +
  geom_histogram(bins = 40, fill = "#4E79A7", color = "white", linewidth = 0.15) +
  theme_light(base_size = 10) +
  theme(panel.grid.minor = element_blank(), plot.title.position = "plot") +
  labs(
    title = "Empirical replicate-averaged p-values",
    subtitle = "Compound/promoter replicate averages compared with matched promoter-library-plate DMSO averages",
    x = "Empirical p-value",
    y = "Count"
  )
ggsave(file.path(out_dir, "empirical_replicate_pvalue_histogram.png"), p_hist, width = 7, height = 4.5, dpi = 300, bg = "white")
ggsave(file.path(out_dir, "empirical_replicate_pvalue_histogram.pdf"), p_hist, width = 7, height = 4.5, bg = "white")

p_perm_hist <- ggplot(emp, aes(permutation_pvalue)) +
  geom_histogram(bins = 50, fill = "#009E73", color = "white", linewidth = 0.15) +
  theme_light(base_size = 10) +
  theme(panel.grid.minor = element_blank(), plot.title.position = "plot") +
  labs(
    title = "Permutation replicate-averaged p-values",
    subtitle = "B = 1000 DMSO replicate-sized draws within matched promoter-library-plate nulls",
    x = "Permutation p-value",
    y = "Count"
  )
ggsave(file.path(out_dir, "permutation_replicate_pvalue_histogram.png"), p_perm_hist, width = 7, height = 4.5, dpi = 300, bg = "white")
ggsave(file.path(out_dir, "permutation_replicate_pvalue_histogram.pdf"), p_perm_hist, width = 7, height = 4.5, bg = "white")

p_exact_perm <- ggplot(emp, aes(empirical_pvalue, permutation_pvalue)) +
  geom_point(aes(color = permutation_bh_hit), alpha = 0.45, size = 0.65) +
  geom_abline(slope = 1, intercept = 0, linetype = "longdash", color = "#303030", linewidth = 0.35) +
  scale_color_manual(values = c("FALSE" = "#9AA0A6", "TRUE" = "#009E73"), guide = "none") +
  theme_light(base_size = 10) +
  theme(panel.grid.minor = element_blank(), plot.title.position = "plot") +
  labs(
    title = "Exact matched empirical vs B=1000 permutation p-values",
    subtitle = "Permutation p-values avoid the 1/(n_null+1) exact empirical grid limit",
    x = "Exact matched empirical p-value",
    y = "Permutation p-value"
  )
ggsave(file.path(out_dir, "exact_empirical_vs_permutation_pvalues.png"), p_exact_perm, width = 6.5, height = 5.5, dpi = 300, bg = "white")
ggsave(file.path(out_dir, "exact_empirical_vs_permutation_pvalues.pdf"), p_exact_perm, width = 6.5, height = 5.5, bg = "white")

p_scatter <- ggplot(emp, aes(estimated_alpha_eb_neglog10padj, permutation_neglog10padj)) +
  geom_point(aes(color = permutation_bh_hit | eb_bh_hit), alpha = 0.45, size = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "longdash", color = "#303030", linewidth = 0.35) +
  scale_color_manual(values = c("FALSE" = "#9AA0A6", "TRUE" = "#B2182B"), guide = "none") +
  theme_light(base_size = 10) +
  theme(panel.grid.minor = element_blank(), plot.title.position = "plot") +
  labs(
    title = "EB BH p-values vs permutation replicate-averaged p-values",
    x = "DStressR EB -log10 adjusted p-value",
    y = "Permutation replicate -log10 adjusted p-value"
  )
ggsave(file.path(out_dir, "eb_vs_permutation_replicate_adjusted_pvalues.png"), p_scatter, width = 6.5, height = 5.5, dpi = 300, bg = "white")
ggsave(file.path(out_dir, "eb_vs_permutation_replicate_adjusted_pvalues.pdf"), p_scatter, width = 6.5, height = 5.5, bg = "white")

message("Wrote empirical replicate p-values to: ", out_dir)
print(summary_counts, row.names = FALSE)
