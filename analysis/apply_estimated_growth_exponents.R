#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
})

source(file.path(getwd(), "R", "growth.R"))

root <- "/Users/cmueller/Documents/GitHub/campylobacter_stressregnet/workflow/data"
out_dir <- file.path(getwd(), "analysis", "outputs", "growth_exponent")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

expression_file <- file.path(root, "02-lux_expression", "expression_values.tsv.gz")
old_pvalue_file <- file.path(root, "03-hit_determination", "expression_df.pvalues.tsv.gz")
libmap_file <- file.path(root, "00-import", "Campylobacter", "LibMap.txt")
fixed_eb_file <- file.path(
  getwd(),
  "analysis",
  "outputs",
  "eb_moderated_variance",
  "workflow_vs_destress_eb_promoter_compound_pvalues.tsv"
)

read_tsv_base <- function(path) {
  read.delim(path, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
}

safe_log10p <- function(p) {
  -log10(pmax(p, .Machine$double.xmin))
}

squeeze_variances <- function(s2, df) {
  keep <- is.finite(s2) & s2 > 0 & is.finite(df) & df > 0
  z <- log(s2[keep])
  observed_var <- stats::var(z)
  target <- max(observed_var - mean(trigamma(df[keep] / 2)), 1e-8)
  df0 <- stats::optimize(function(x) (trigamma(x / 2) - target)^2, interval = c(0.1, 1e5))$minimum
  s20 <- exp(mean(z - digamma(df[keep] / 2) + log(df[keep] / 2)))
  list(df0 = df0, s20 = s20)
}

run_destress_eb <- function(expr, dmso_srn_codes, old_pvalue_file) {
  dmso_expr <- expr[expr$srn_code %in% dmso_srn_codes, ]
  dmso_mean <- aggregate(
    normalized.lux ~ promoter_libplate_replicate,
    dmso_expr,
    mean
  )
  names(dmso_mean)[2] <- "dmso_mean_normalized_lux"

  modeled <- merge(expr, dmso_mean, by = "promoter_libplate_replicate", all.x = TRUE, sort = FALSE)
  modeled$log2FC <- modeled$normalized.lux - modeled$dmso_mean_normalized_lux

  compound_global <- aggregate(log2FC ~ srn_code, modeled, mean, na.rm = TRUE)
  names(compound_global)[2] <- "destress_global_effect"
  modeled <- merge(modeled, compound_global, by = "srn_code", all.x = TRUE, sort = FALSE)
  modeled$destress_specific_effect <- modeled$log2FC - modeled$destress_global_effect

  dmso_specific <- modeled[modeled$srn_code %in% dmso_srn_codes, ]
  group_mean <- aggregate(
    destress_specific_effect ~ promoter_libplate_replicate,
    dmso_specific,
    mean,
    na.rm = TRUE
  )
  names(group_mean)[2] <- "eb_dmso_group_mean"
  dmso_specific <- merge(dmso_specific, group_mean, by = "promoter_libplate_replicate", all.x = TRUE, sort = FALSE)
  dmso_specific$dmso_centered_specific <- dmso_specific$destress_specific_effect -
    dmso_specific$eb_dmso_group_mean

  promoter_variance <- do.call(
    rbind,
    lapply(split(dmso_specific, dmso_specific$promoter), function(d) {
      df <- nrow(d) - length(unique(d$promoter_libplate_replicate))
      s2 <- sum(d$dmso_centered_specific^2, na.rm = TRUE) / df
      data.frame(
        promoter = d$promoter[1],
        dmso_n = nrow(d),
        dmso_groups = length(unique(d$promoter_libplate_replicate)),
        variance_df = df,
        promoter_s2 = s2,
        stringsAsFactors = FALSE
      )
    })
  )
  rownames(promoter_variance) <- NULL
  prior <- squeeze_variances(promoter_variance$promoter_s2, promoter_variance$variance_df)
  promoter_variance$eb_prior_df <- prior$df0
  promoter_variance$eb_prior_s2 <- prior$s20
  promoter_variance$moderated_s2 <- (
    prior$df0 * prior$s20 + promoter_variance$variance_df * promoter_variance$promoter_s2
  ) / (prior$df0 + promoter_variance$variance_df)
  promoter_variance$moderated_df <- prior$df0 + promoter_variance$variance_df

  modeled <- merge(modeled, group_mean, by = "promoter_libplate_replicate", all.x = TRUE, sort = FALSE)
  modeled <- merge(
    modeled,
    promoter_variance[, c("promoter", "promoter_s2", "variance_df", "moderated_s2", "moderated_df")],
    by = "promoter",
    all.x = TRUE,
    sort = FALSE
  )
  modeled$destress_eb_effect_centered <- modeled$destress_specific_effect -
    modeled$eb_dmso_group_mean
  modeled$destress_eb_t <- modeled$destress_eb_effect_centered / sqrt(modeled$moderated_s2)
  modeled$destress_eb_pvalue <- 2 * stats::pt(
    abs(modeled$destress_eb_t),
    df = modeled$moderated_df,
    lower.tail = FALSE
  )

  old <- read_tsv_base(old_pvalue_file)
  comparison <- merge(
    old,
    modeled[, c(
      "promoter_libplate_replicate",
      "promoter",
      "libplate",
      "replicate",
      "srn_code",
      "normalized.lux",
      "growth_alpha",
      "log2FC",
      "destress_global_effect",
      "destress_specific_effect",
      "destress_eb_effect_centered",
      "destress_eb_t",
      "destress_eb_pvalue"
    )],
    by = c("promoter_libplate_replicate", "promoter", "libplate", "replicate", "srn_code"),
    all = FALSE,
    sort = FALSE
  )

  list(replicate = comparison, promoter_variance = promoter_variance, prior = prior)
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

expr <- read_tsv_base(expression_file)
expr <- merge(
  expr,
  libmap[, c("srn_code", "ProductName", "Catalog Number")],
  by = "srn_code",
  all.x = TRUE,
  sort = FALSE
)
expr <- expr[
  !(expr[["Catalog Number"]] %in% "DMSO noisy") &
    !(expr$promoter %in% c("PCJnc20", "PCjas704")) &
    is.finite(expr$LUX.AUC_16) &
    is.finite(expr$od_16h.measured) &
    expr$LUX.AUC_16 > 0 &
    expr$od_16h.measured > 0,
]

alpha_est <- estimate_growth_exponents(
  expr,
  promoter = "promoter",
  compound = "srn_code",
  lux = "LUX.AUC_16",
  growth = "od_16h.measured",
  covariates = c("libplate", "replicate"),
  controls = dmso_srn_codes,
  min_control_n = 20,
  shrink = TRUE,
  alpha_bounds = c(-2, 3)
)
alpha_est$alpha_raw_lower <- alpha_est$alpha_raw - 1.96 * alpha_est$alpha_raw_se
alpha_est$alpha_raw_upper <- alpha_est$alpha_raw + 1.96 * alpha_est$alpha_raw_se
alpha_est$alpha_shrunk_lower <- alpha_est$alpha_shrunk - 1.96 * alpha_est$alpha_shrunk_se
alpha_est$alpha_shrunk_upper <- alpha_est$alpha_shrunk + 1.96 * alpha_est$alpha_shrunk_se
alpha_est <- alpha_est[order(-alpha_est$alpha_shrunk), ]
alpha_est$promoter <- factor(alpha_est$promoter, levels = alpha_est$promoter)

write.table(
  alpha_est,
  file.path(out_dir, "growth_exponent_estimates.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

expr <- expr[order(expr$promoter, expr$srn_code), ]
expr$original_replicate <- expr$replicate
expr$experiment_id_original <- expr$experiment_id
expr$replicate <- rep(c("r1", "r2"), length.out = nrow(expr))
expr$promoter_libplate_replicate <- paste(expr$promoter, expr$libplate, expr$replicate, sep = "_")
expr$growth_alpha <- alpha_est$alpha_shrunk[match(expr$promoter, as.character(alpha_est$promoter))]
expr$normalized.lux <- log2(expr$LUX.AUC_16) - expr$growth_alpha * log2(expr$od_16h.measured)
expr$fixed_alpha_response <- log2(expr$LUX.AUC_16) - log2(expr$od_16h.measured)
expr$response_delta_vs_fixed_alpha <- expr$normalized.lux - expr$fixed_alpha_response

write.table(
  expr[, c(
    "experiment_id",
    "promoter",
    "srn_code",
    "libplate",
    "replicate",
    "original_replicate",
    "LUX.AUC_16",
    "od_16h.measured",
    "growth_alpha",
    "fixed_alpha_response",
    "normalized.lux",
    "response_delta_vs_fixed_alpha"
  )],
  file.path(out_dir, "expression_values_estimated_growth_alpha.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

fit <- run_destress_eb(expr, dmso_srn_codes, old_pvalue_file)
replicate_level <- fit$replicate
replicate_level <- replicate_level[is.finite(replicate_level$destress_eb_pvalue), ]
replicate_level$estimated_alpha_neglog10p <- safe_log10p(replicate_level$destress_eb_pvalue)

write.table(
  replicate_level,
  file.path(out_dir, "workflow_vs_destress_eb_estimated_growth_alpha_replicate_pvalues.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

compound_level <- replicate_level[!(replicate_level$srn_code %in% c(dmso_srn_codes, dmso_noisy_srn_codes)), ]
compound_level <- compound_level[order(compound_level$promoter, compound_level$srn_code, -compound_level$pvalue), ]
compound_level <- compound_level[!duplicated(paste(compound_level$promoter, compound_level$srn_code)), ]
compound_level$estimated_alpha_eb_padj_by_promoter <- ave(
  compound_level$destress_eb_pvalue,
  compound_level$promoter,
  FUN = function(x) p.adjust(x, method = "BH")
)
compound_level$estimated_alpha_eb_neglog10padj <- safe_log10p(compound_level$estimated_alpha_eb_padj_by_promoter)

if (file.exists(fixed_eb_file)) {
  fixed_eb <- read_tsv_base(fixed_eb_file)
  fixed_cols <- fixed_eb[, c("promoter", "srn_code", "destress_eb_pvalue", "destress_eb_padj_by_promoter")]
  names(fixed_cols) <- c("promoter", "srn_code", "fixed_alpha_eb_pvalue", "fixed_alpha_eb_padj_by_promoter")
  compound_level <- merge(compound_level, fixed_cols, by = c("promoter", "srn_code"), all.x = TRUE, sort = FALSE)
  compound_level$fixed_alpha_eb_neglog10p <- safe_log10p(compound_level$fixed_alpha_eb_pvalue)
  compound_level$fixed_alpha_eb_neglog10padj <- safe_log10p(compound_level$fixed_alpha_eb_padj_by_promoter)
}

write.table(
  compound_level,
  file.path(out_dir, "workflow_vs_destress_eb_estimated_growth_alpha_promoter_compound_pvalues.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

p_alpha <- ggplot(alpha_est, aes(promoter, alpha_shrunk)) +
  geom_hline(yintercept = 1, color = "#b91c1c", linewidth = 0.45, linetype = "dashed") +
  geom_hline(yintercept = unique(alpha_est$alpha_global)[1], color = "#334155", linewidth = 0.45, linetype = "dotted") +
  geom_errorbar(aes(ymin = alpha_shrunk_lower, ymax = alpha_shrunk_upper), width = 0.2, color = "#16a34a", alpha = 0.75) +
  geom_point(aes(y = alpha_raw), color = "#64748b", alpha = 0.55, size = 1.8) +
  geom_point(color = "#16a34a", size = 2.2) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1), legend.position = "none") +
  labs(
    title = "Estimated promoter-specific growth exponents",
    subtitle = "Adjusted for libplate and replicate; green: EB-shrunken alpha_g with 95% intervals; gray: raw adjusted slope; red dashed: alpha = 1",
    x = "Promoter, ordered high-to-low by shrunken alpha_g",
    y = expression(alpha[g])
  )

ggsave(file.path(out_dir, "growth_exponent_alpha_estimates.png"), p_alpha, width = 11, height = 5.5, dpi = 300, bg = "white")
ggsave(file.path(out_dir, "growth_exponent_alpha_estimates.pdf"), p_alpha, width = 11, height = 5.5, bg = "white")

alpha_plot_df <- alpha_est
alpha_plot_df$significant_vs_one <- alpha_plot_df$alpha_shrunk_lower > 1 | alpha_plot_df$alpha_shrunk_upper < 1
p_alpha_scatter <- ggplot(alpha_plot_df, aes(alpha_raw, alpha_shrunk, label = promoter, color = significant_vs_one)) +
  geom_abline(slope = 1, intercept = 0, color = "#525252", linewidth = 0.4) +
  geom_hline(yintercept = 1, color = "#b91c1c", linetype = "dashed", linewidth = 0.35) +
  geom_vline(xintercept = 1, color = "#b91c1c", linetype = "dashed", linewidth = 0.35) +
  geom_point(size = 2.2) +
  geom_text_repel(size = 2.5, max.overlaps = 30) +
  scale_color_manual(values = c("FALSE" = "#64748b", "TRUE" = "#16a34a")) +
  theme_bw(base_size = 10) +
  labs(
    title = "Raw vs shrunken growth exponents",
    x = "Raw DMSO slope",
    y = "Shrunken alpha_g",
    color = "95% interval excludes 1"
  )
ggsave(file.path(out_dir, "growth_exponent_raw_vs_shrunken.png"), p_alpha_scatter, width = 7, height = 6, dpi = 300, bg = "white")
ggsave(file.path(out_dir, "growth_exponent_raw_vs_shrunken.pdf"), p_alpha_scatter, width = 7, height = 6, bg = "white")

response_sample <- expr
if (nrow(response_sample) > 50000) {
  set.seed(31)
  response_sample <- response_sample[sample.int(nrow(response_sample), 50000), ]
}
p_response <- ggplot(response_sample, aes(fixed_alpha_response, normalized.lux)) +
  geom_point(alpha = 0.12, size = 0.35, color = "#1f2937") +
  geom_abline(slope = 1, intercept = 0, color = "#b91c1c", linewidth = 0.45) +
  theme_bw(base_size = 10) +
  coord_equal() +
  labs(
    title = "Response change from estimating growth exponents",
    x = expression(log[2](LUX.AUC / OD)~"(fixed alpha = 1)"),
    y = expression(log[2](LUX.AUC) - hat(alpha)[g]~log[2](OD))
  )
ggsave(file.path(out_dir, "response_fixed_alpha_vs_estimated_alpha.png"), p_response, width = 6.5, height = 6, dpi = 300, bg = "white")
ggsave(file.path(out_dir, "response_fixed_alpha_vs_estimated_alpha.pdf"), p_response, width = 6.5, height = 6, bg = "white")

if ("fixed_alpha_eb_pvalue" %in% names(compound_level)) {
  p_pvalue <- ggplot(compound_level, aes(fixed_alpha_eb_neglog10p, estimated_alpha_neglog10p)) +
    geom_point(alpha = 0.12, size = 0.45, color = "#1f2937") +
    geom_abline(slope = 1, intercept = 0, color = "#b91c1c", linewidth = 0.45) +
    coord_equal() +
    theme_bw(base_size = 10) +
    labs(
      title = "EB p-values: fixed alpha = 1 vs estimated alpha_g",
      x = "Fixed-alpha DStressR EB -log10(p)",
      y = "Estimated-alpha DStressR EB -log10(p)"
    )
  ggsave(file.path(out_dir, "pvalues_fixed_alpha_vs_estimated_alpha.png"), p_pvalue, width = 6.5, height = 6, dpi = 300, bg = "white")
  ggsave(file.path(out_dir, "pvalues_fixed_alpha_vs_estimated_alpha.pdf"), p_pvalue, width = 6.5, height = 6, bg = "white")

  p_padj <- ggplot(compound_level, aes(fixed_alpha_eb_neglog10padj, estimated_alpha_eb_neglog10padj)) +
    geom_point(alpha = 0.15, size = 0.5, color = "#1f2937") +
    geom_abline(slope = 1, intercept = 0, color = "#b91c1c", linewidth = 0.45) +
    coord_equal() +
    theme_bw(base_size = 10) +
    labs(
      title = "EB adjusted p-values: fixed alpha = 1 vs estimated alpha_g",
      x = "Fixed-alpha DStressR EB -log10(BH adjusted p)",
      y = "Estimated-alpha DStressR EB -log10(BH adjusted p)"
    )
  ggsave(file.path(out_dir, "padj_fixed_alpha_vs_estimated_alpha.png"), p_padj, width = 6.5, height = 6, dpi = 300, bg = "white")
  ggsave(file.path(out_dir, "padj_fixed_alpha_vs_estimated_alpha.pdf"), p_padj, width = 6.5, height = 6, bg = "white")
}

summary_table <- data.frame(
  metric = c(
    "promoters_with_alpha_estimates",
    "global_alpha",
    "alpha_prior_sd",
    "median_alpha_raw",
    "median_alpha_shrunk",
    "promoters_alpha_ci_excludes_one",
    "estimated_alpha_eb_bh_lt_0.05",
    "fixed_alpha_eb_bh_lt_0.05",
    "both_bh_lt_0.05",
    "estimated_alpha_only_bh_lt_0.05",
    "fixed_alpha_only_bh_lt_0.05",
    "spearman_fixed_vs_estimated_neglog10p"
  ),
  value = c(
    nrow(alpha_est),
    unique(alpha_est$alpha_global)[1],
    unique(alpha_est$alpha_prior_sd)[1],
    stats::median(alpha_est$alpha_raw, na.rm = TRUE),
    stats::median(alpha_est$alpha_shrunk, na.rm = TRUE),
    sum(alpha_est$alpha_shrunk_lower > 1 | alpha_est$alpha_shrunk_upper < 1, na.rm = TRUE),
    sum(compound_level$estimated_alpha_eb_padj_by_promoter < 0.05, na.rm = TRUE),
    if ("fixed_alpha_eb_padj_by_promoter" %in% names(compound_level)) sum(compound_level$fixed_alpha_eb_padj_by_promoter < 0.05, na.rm = TRUE) else NA_real_,
    if ("fixed_alpha_eb_padj_by_promoter" %in% names(compound_level)) sum(compound_level$estimated_alpha_eb_padj_by_promoter < 0.05 & compound_level$fixed_alpha_eb_padj_by_promoter < 0.05, na.rm = TRUE) else NA_real_,
    if ("fixed_alpha_eb_padj_by_promoter" %in% names(compound_level)) sum(compound_level$estimated_alpha_eb_padj_by_promoter < 0.05 & !(compound_level$fixed_alpha_eb_padj_by_promoter < 0.05), na.rm = TRUE) else NA_real_,
    if ("fixed_alpha_eb_padj_by_promoter" %in% names(compound_level)) sum(!(compound_level$estimated_alpha_eb_padj_by_promoter < 0.05) & compound_level$fixed_alpha_eb_padj_by_promoter < 0.05, na.rm = TRUE) else NA_real_,
    if ("fixed_alpha_eb_pvalue" %in% names(compound_level)) stats::cor(compound_level$fixed_alpha_eb_neglog10p, compound_level$estimated_alpha_neglog10p, method = "spearman", use = "complete.obs") else NA_real_
  )
)

write.table(
  summary_table,
  file.path(out_dir, "estimated_growth_alpha_summary.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

print(summary_table)
message("Wrote estimated growth-alpha outputs to: ", out_dir)
