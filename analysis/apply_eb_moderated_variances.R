#!/usr/bin/env Rscript

source(file.path("analysis", "_helpers.R"))

suppressPackageStartupMessages({
  library(ggplot2)
})

root <- analysis_data_root()
out_dir <- file.path(getwd(), "analysis", "outputs", "eb_moderated_variance")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

comparison_file <- file.path(getwd(), "analysis", "outputs", "workflow_vs_destress_replicate_pvalues.tsv")
libmap_file <- file.path(root, "00-import", "Campylobacter", "LibMap.txt")

read_tsv_base <- function(path) {
  read.delim(path, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
}

safe_log10p <- function(p) {
  -log10(pmax(p, .Machine$double.xmin))
}

squeeze_variances <- function(s2, df) {
  keep <- is.finite(s2) & s2 > 0 & is.finite(df) & df > 0
  if (sum(keep) < 3) {
    stop("Need at least three finite variance estimates for empirical-Bayes moderation.", call. = FALSE)
  }
  z <- log(s2[keep])
  observed_var <- stats::var(z)
  target <- max(observed_var - mean(trigamma(df[keep] / 2)), 1e-8)
  objective <- function(df0) (trigamma(df0 / 2) - target)^2
  df0 <- stats::optimize(objective, interval = c(0.1, 1e5))$minimum
  s20 <- exp(mean(z - digamma(df[keep] / 2) + log(df[keep] / 2)))
  list(df0 = df0, s20 = s20)
}

libmap <- read_tsv_base(libmap_file)
libmap$libplate <- paste0("lp", libmap[["Library plate"]])
libmap$srn_code <- paste(libmap$libplate, libmap[["Well"]], sep = "_")
libmap$ProductName <- ifelse(
  is.na(libmap$ProductName) | libmap$ProductName == "NA" | libmap$ProductName == "",
  libmap[["Catalog Number"]],
  libmap$ProductName
)
dmso_srn_codes <- libmap$srn_code[libmap$ProductName == "DMSO"]
dmso_noisy_srn_codes <- libmap$srn_code[libmap$ProductName == "DMSO noisy"]

comparison <- read_tsv_base(comparison_file)
comparison <- comparison[is.finite(comparison$destress_specific_effect), ]

dmso <- comparison[comparison$srn_code %in% dmso_srn_codes, ]

group_mean <- aggregate(
  destress_specific_effect ~ promoter_libplate_replicate,
  dmso,
  mean,
  na.rm = TRUE
)
names(group_mean)[2] <- "eb_dmso_group_mean"

dmso <- merge(dmso, group_mean, by = "promoter_libplate_replicate", all.x = TRUE, sort = FALSE)
dmso$dmso_centered_specific <- dmso$destress_specific_effect - dmso$eb_dmso_group_mean

promoter_variance <- do.call(
  rbind,
  lapply(split(dmso, dmso$promoter), function(d) {
    df <- nrow(d) - length(unique(d$promoter_libplate_replicate))
    s2 <- sum(d$dmso_centered_specific^2, na.rm = TRUE) / df
    data.frame(
      promoter = d$promoter[1],
      dmso_n = nrow(d),
      dmso_groups = length(unique(d$promoter_libplate_replicate)),
      variance_df = df,
      promoter_s2 = s2
    )
  })
)
rownames(promoter_variance) <- NULL

prior <- squeeze_variances(promoter_variance$promoter_s2, promoter_variance$variance_df)
promoter_variance$eb_prior_df <- prior$df0
promoter_variance$eb_prior_s2 <- prior$s20
promoter_variance$moderated_s2 <- (
  prior$df0 * prior$s20 +
    promoter_variance$variance_df * promoter_variance$promoter_s2
) / (prior$df0 + promoter_variance$variance_df)
promoter_variance$moderated_df <- prior$df0 + promoter_variance$variance_df
promoter_variance$variance_shrinkage_ratio <- promoter_variance$moderated_s2 / promoter_variance$promoter_s2

write.table(
  promoter_variance,
  file.path(out_dir, "promoter_empirical_bayes_variances.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

comparison <- merge(
  comparison,
  group_mean,
  by = "promoter_libplate_replicate",
  all.x = TRUE,
  sort = FALSE
)
comparison <- merge(
  comparison,
  promoter_variance[, c("promoter", "promoter_s2", "variance_df", "moderated_s2", "moderated_df")],
  by = "promoter",
  all.x = TRUE,
  sort = FALSE
)

comparison$destress_eb_effect_centered <- comparison$destress_specific_effect -
  comparison$eb_dmso_group_mean
comparison$destress_eb_t <- comparison$destress_eb_effect_centered / sqrt(comparison$moderated_s2)
comparison$destress_eb_pvalue <- 2 * stats::pt(
  abs(comparison$destress_eb_t),
  df = comparison$moderated_df,
  lower.tail = FALSE
)
comparison$destress_eb_neglog10p <- safe_log10p(comparison$destress_eb_pvalue)

write.table(
  comparison,
  file.path(out_dir, "workflow_vs_destress_eb_replicate_pvalues.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

compound_level <- comparison[!(comparison$srn_code %in% c(dmso_srn_codes, dmso_noisy_srn_codes)), ]
compound_level <- compound_level[order(compound_level$promoter, compound_level$srn_code, -compound_level$pvalue), ]
compound_level <- compound_level[!duplicated(paste(compound_level$promoter, compound_level$srn_code)), ]

compound_level$workflow_padj_by_promoter <- ave(
  compound_level$pvalue,
  compound_level$promoter,
  FUN = function(x) p.adjust(x, method = "BH")
)
compound_level$destress_gaussian_padj_by_promoter <- ave(
  compound_level$destress_pvalue,
  compound_level$promoter,
  FUN = function(x) p.adjust(x, method = "BH")
)
compound_level$destress_eb_padj_by_promoter <- ave(
  compound_level$destress_eb_pvalue,
  compound_level$promoter,
  FUN = function(x) p.adjust(x, method = "BH")
)

compound_level$workflow_neglog10p <- safe_log10p(compound_level$pvalue)
compound_level$destress_gaussian_neglog10p <- safe_log10p(compound_level$destress_pvalue)
compound_level$destress_eb_neglog10p <- safe_log10p(compound_level$destress_eb_pvalue)
compound_level$workflow_neglog10padj <- safe_log10p(compound_level$workflow_padj_by_promoter)
compound_level$destress_gaussian_neglog10padj <- safe_log10p(compound_level$destress_gaussian_padj_by_promoter)
compound_level$destress_eb_neglog10padj <- safe_log10p(compound_level$destress_eb_padj_by_promoter)

write.table(
  compound_level,
  file.path(out_dir, "workflow_vs_destress_eb_promoter_compound_pvalues.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

summary_table <- data.frame(
  metric = c(
    "promoters_with_variance_estimates",
    "eb_prior_df",
    "eb_prior_sd",
    "median_raw_promoter_sd",
    "median_moderated_promoter_sd",
    "median_variance_shrinkage_ratio",
    "promoter_compound_rows",
    "workflow_raw_p_lt_0.05",
    "destress_gaussian_raw_p_lt_0.05",
    "destress_eb_raw_p_lt_0.05",
    "workflow_bh_lt_0.05",
    "destress_gaussian_bh_lt_0.05",
    "destress_eb_bh_lt_0.05",
    "workflow_vs_eb_spearman_neglog10p",
    "gaussian_vs_eb_spearman_neglog10p"
  ),
  value = c(
    nrow(promoter_variance),
    prior$df0,
    sqrt(prior$s20),
    stats::median(sqrt(promoter_variance$promoter_s2), na.rm = TRUE),
    stats::median(sqrt(promoter_variance$moderated_s2), na.rm = TRUE),
    stats::median(promoter_variance$variance_shrinkage_ratio, na.rm = TRUE),
    nrow(compound_level),
    sum(compound_level$pvalue < 0.05, na.rm = TRUE),
    sum(compound_level$destress_pvalue < 0.05, na.rm = TRUE),
    sum(compound_level$destress_eb_pvalue < 0.05, na.rm = TRUE),
    sum(compound_level$workflow_padj_by_promoter < 0.05, na.rm = TRUE),
    sum(compound_level$destress_gaussian_padj_by_promoter < 0.05, na.rm = TRUE),
    sum(compound_level$destress_eb_padj_by_promoter < 0.05, na.rm = TRUE),
    stats::cor(compound_level$workflow_neglog10p, compound_level$destress_eb_neglog10p, method = "spearman"),
    stats::cor(compound_level$destress_gaussian_neglog10p, compound_level$destress_eb_neglog10p, method = "spearman")
  )
)

write.table(
  summary_table,
  file.path(out_dir, "eb_moderated_variance_summary.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

p_var <- ggplot(promoter_variance, aes(sqrt(promoter_s2), sqrt(moderated_s2))) +
  geom_point(color = "#1f2937", size = 2) +
  geom_abline(slope = 1, intercept = 0, color = "#b91c1c", linewidth = 0.5) +
  theme_bw(base_size = 10) +
  labs(
    title = "Promoter-specific DMSO SD before and after EB moderation",
    subtitle = paste0("Prior df = ", round(prior$df0, 2), "; prior SD = ", round(sqrt(prior$s20), 3)),
    x = "Raw promoter DMSO residual SD",
    y = "Moderated promoter DMSO residual SD"
  )
ggsave(file.path(out_dir, "promoter_sd_raw_vs_moderated.png"), p_var, width = 6, height = 5, dpi = 300)
ggsave(file.path(out_dir, "promoter_sd_raw_vs_moderated.pdf"), p_var, width = 6, height = 5)

p_scatter <- ggplot(compound_level, aes(destress_gaussian_neglog10p, destress_eb_neglog10p)) +
  geom_point(alpha = 0.12, size = 0.45, color = "#1f2937") +
  geom_abline(slope = 1, intercept = 0, color = "#b91c1c", linewidth = 0.5) +
  coord_equal() +
  theme_bw(base_size = 10) +
  labs(
    title = "DStressR Gaussian vs EB-moderated promoter variance p-values",
    x = "DStressR per-promoter-library-replicate Gaussian -log10(p)",
    y = "DStressR EB moderated promoter-variance -log10(p)"
  )
ggsave(file.path(out_dir, "scatter_destress_gaussian_vs_eb_raw_pvalues.png"), p_scatter, width = 6.5, height = 6, dpi = 300)
ggsave(file.path(out_dir, "scatter_destress_gaussian_vs_eb_raw_pvalues.pdf"), p_scatter, width = 6.5, height = 6)

p_adjusted <- ggplot(compound_level, aes(destress_gaussian_neglog10padj, destress_eb_neglog10padj)) +
  geom_point(alpha = 0.15, size = 0.5, color = "#1f2937") +
  geom_abline(slope = 1, intercept = 0, color = "#b91c1c", linewidth = 0.5) +
  coord_equal() +
  theme_bw(base_size = 10) +
  labs(
    title = "DStressR Gaussian vs EB-moderated promoter variance adjusted p-values",
    x = "DStressR Gaussian -log10(BH adjusted p)",
    y = "DStressR EB moderated -log10(BH adjusted p)"
  )
ggsave(file.path(out_dir, "scatter_destress_gaussian_vs_eb_adjusted_pvalues.png"), p_adjusted, width = 6.5, height = 6, dpi = 300)
ggsave(file.path(out_dir, "scatter_destress_gaussian_vs_eb_adjusted_pvalues.pdf"), p_adjusted, width = 6.5, height = 6)

p_workflow <- ggplot(compound_level, aes(workflow_neglog10p, destress_eb_neglog10p)) +
  geom_point(alpha = 0.12, size = 0.45, color = "#1f2937") +
  geom_abline(slope = 1, intercept = 0, color = "#b91c1c", linewidth = 0.5) +
  coord_equal() +
  theme_bw(base_size = 10) +
  labs(
    title = "Workflow median-polish vs DStressR EB-moderated p-values",
    x = "Workflow median-polish -log10(p)",
    y = "DStressR EB moderated -log10(p)"
  )
ggsave(file.path(out_dir, "scatter_workflow_vs_destress_eb_raw_pvalues.png"), p_workflow, width = 6.5, height = 6, dpi = 300)
ggsave(file.path(out_dir, "scatter_workflow_vs_destress_eb_raw_pvalues.pdf"), p_workflow, width = 6.5, height = 6)

print(summary_table)
message("Wrote EB moderated variance outputs to: ", out_dir)
