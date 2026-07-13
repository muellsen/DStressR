#!/usr/bin/env Rscript

source(file.path("analysis", "_helpers.R"))

suppressPackageStartupMessages({
  library(ggplot2)
})

root <- analysis_data_root()
out_dir <- file.path(getwd(), "analysis", "outputs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

expression_file <- file.path(root, "02-lux_expression", "expression_values.tsv.gz")
old_pvalue_file <- file.path(root, "03-hit_determination", "expression_df.pvalues.tsv.gz")
libmap_file <- file.path(root, "00-import", "Campylobacter", "LibMap.txt")

read_tsv_base <- function(path) {
  read.delim(path, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
}

safe_log10p <- function(p) {
  -log10(pmax(p, .Machine$double.xmin))
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
expr$normalized.lux <- expr$log2.auc.16hmeasured.normed
expr <- merge(
  expr,
  libmap[, c("srn_code", "ProductName", "Catalog Number")],
  by = "srn_code",
  all.x = TRUE,
  sort = FALSE
)

# Reproduce the workflow's filtering and two-replicate recoding before the old
# median-polish p-values were generated.
expr <- expr[
  expr[["Catalog Number"]] != "DMSO noisy" &
    !(expr$promoter %in% c("PCJnc20", "PCjas704")),
]
expr <- expr[order(expr$promoter, expr$srn_code), ]
expr$original_replicate <- expr$replicate
expr$experiment_id_original <- expr$experiment_id
expr$replicate <- rep(c("r1", "r2"), length.out = nrow(expr))
expr$promoter_libplate_replicate <- paste(expr$promoter, expr$libplate, expr$replicate, sep = "_")

dmso_expr <- expr[expr$srn_code %in% dmso_srn_codes, ]
dmso_mean <- aggregate(
  normalized.lux ~ promoter_libplate_replicate,
  dmso_expr,
  mean
)
names(dmso_mean)[2] <- "dmso_mean_normalized_lux"

modeled <- merge(expr, dmso_mean, by = "promoter_libplate_replicate", all.x = TRUE, sort = FALSE)
modeled$log2FC <- modeled$normalized.lux - modeled$dmso_mean_normalized_lux

# Scalable DStressR-style decomposition:
# log2FC[g,c,r] = compound_global[c] + promoter_specific[g,c,r] + residual.
# This is the same estimand as the fixed-effect framework: separate the
# compound-wide perturbation from promoter-specific stress.
compound_global <- aggregate(log2FC ~ srn_code, modeled, mean, na.rm = TRUE)
names(compound_global)[2] <- "destress_global_effect"
modeled <- merge(modeled, compound_global, by = "srn_code", all.x = TRUE, sort = FALSE)
modeled$destress_specific_effect <- modeled$log2FC - modeled$destress_global_effect

dmso_specific <- modeled[modeled$srn_code %in% dmso_srn_codes, ]
dmso_params <- aggregate(
  destress_specific_effect ~ promoter_libplate_replicate,
  dmso_specific,
  function(x) c(mean = mean(x, na.rm = TRUE), sd = stats::sd(x, na.rm = TRUE), n = sum(!is.na(x)))
)
dmso_params <- do.call(data.frame, dmso_params)
names(dmso_params) <- c(
  "promoter_libplate_replicate",
  "destress_dmso_mean",
  "destress_dmso_sd",
  "destress_dmso_n"
)

modeled <- merge(modeled, dmso_params, by = "promoter_libplate_replicate", all.x = TRUE, sort = FALSE)
modeled$destress_z <- (
  modeled$destress_specific_effect - modeled$destress_dmso_mean
) / modeled$destress_dmso_sd
modeled$destress_pvalue <- 2 * pnorm(abs(modeled$destress_z), lower.tail = FALSE)

old <- read_tsv_base(old_pvalue_file)
comparison <- merge(
  old,
  modeled[, c(
    "promoter_libplate_replicate",
    "promoter",
    "libplate",
    "replicate",
    "srn_code",
    "log2FC",
    "destress_global_effect",
    "destress_specific_effect",
    "destress_z",
    "destress_pvalue"
  )],
  by = c("promoter_libplate_replicate", "promoter", "libplate", "replicate", "srn_code"),
  suffixes = c("_workflow", "_destress"),
  all = FALSE,
  sort = FALSE
)

comparison <- comparison[is.finite(comparison$pvalue) & is.finite(comparison$destress_pvalue), ]
comparison$workflow_neglog10p <- safe_log10p(comparison$pvalue)
comparison$destress_neglog10p <- safe_log10p(comparison$destress_pvalue)
comparison$delta_neglog10p <- comparison$destress_neglog10p - comparison$workflow_neglog10p

write.table(
  comparison,
  file.path(out_dir, "workflow_vs_destress_replicate_pvalues.tsv"),
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
compound_level$destress_padj_by_promoter <- ave(
  compound_level$destress_pvalue,
  compound_level$promoter,
  FUN = function(x) p.adjust(x, method = "BH")
)
compound_level$workflow_neglog10padj <- safe_log10p(compound_level$workflow_padj_by_promoter)
compound_level$destress_neglog10padj <- safe_log10p(compound_level$destress_padj_by_promoter)

write.table(
  compound_level,
  file.path(out_dir, "workflow_vs_destress_promoter_compound_pvalues.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

cor_raw <- stats::cor(
  comparison$workflow_neglog10p,
  comparison$destress_neglog10p,
  method = "spearman"
)
cor_compound <- stats::cor(
  compound_level$workflow_neglog10p,
  compound_level$destress_neglog10p,
  method = "spearman"
)

scatter_raw <- ggplot(comparison, aes(workflow_neglog10p, destress_neglog10p)) +
  geom_point(alpha = 0.08, size = 0.35, color = "#1f2937") +
  geom_abline(slope = 1, intercept = 0, color = "#b91c1c", linewidth = 0.5) +
  coord_equal() +
  theme_bw(base_size = 10) +
  labs(
    title = "Replicate-level p-value comparison",
    subtitle = paste0("Spearman rho = ", round(cor_raw, 3), "; n = ", nrow(comparison)),
    x = "Workflow median-polish -log10(p)",
    y = "DStressR-style global-adjusted -log10(p)"
  )

ggsave(
  file.path(out_dir, "scatter_replicate_raw_pvalues.png"),
  scatter_raw,
  width = 6.5,
  height = 6,
  dpi = 300
)
ggsave(
  file.path(out_dir, "scatter_replicate_raw_pvalues.pdf"),
  scatter_raw,
  width = 6.5,
  height = 6
)

scatter_compound <- ggplot(compound_level, aes(workflow_neglog10p, destress_neglog10p)) +
  geom_point(alpha = 0.12, size = 0.45, color = "#1f2937") +
  geom_abline(slope = 1, intercept = 0, color = "#b91c1c", linewidth = 0.5) +
  coord_equal() +
  theme_bw(base_size = 10) +
  labs(
    title = "Promoter-compound p-value comparison",
    subtitle = paste0("Old workflow max-p replicate selection; Spearman rho = ", round(cor_compound, 3)),
    x = "Workflow median-polish -log10(p)",
    y = "DStressR-style global-adjusted -log10(p)"
  )

ggsave(
  file.path(out_dir, "scatter_promoter_compound_raw_pvalues.png"),
  scatter_compound,
  width = 6.5,
  height = 6,
  dpi = 300
)
ggsave(
  file.path(out_dir, "scatter_promoter_compound_raw_pvalues.pdf"),
  scatter_compound,
  width = 6.5,
  height = 6
)

scatter_adjusted <- ggplot(compound_level, aes(workflow_neglog10padj, destress_neglog10padj)) +
  geom_point(alpha = 0.15, size = 0.5, color = "#1f2937") +
  geom_abline(slope = 1, intercept = 0, color = "#b91c1c", linewidth = 0.5) +
  coord_equal() +
  theme_bw(base_size = 10) +
  labs(
    title = "Promoter-compound adjusted p-value comparison",
    x = "Workflow median-polish -log10(BH adjusted p)",
    y = "DStressR-style global-adjusted -log10(BH adjusted p)"
  )

ggsave(
  file.path(out_dir, "scatter_promoter_compound_adjusted_pvalues.png"),
  scatter_adjusted,
  width = 6.5,
  height = 6,
  dpi = 300
)
ggsave(
  file.path(out_dir, "scatter_promoter_compound_adjusted_pvalues.pdf"),
  scatter_adjusted,
  width = 6.5,
  height = 6
)

summary_table <- data.frame(
  metric = c(
    "replicate_rows_compared",
    "promoter_compound_rows_compared",
    "replicate_spearman_neglog10p",
    "promoter_compound_spearman_neglog10p",
    "workflow_raw_p_lt_0.05_replicate",
    "destress_raw_p_lt_0.05_replicate",
    "workflow_bh_lt_0.05_promoter_compound",
    "destress_bh_lt_0.05_promoter_compound"
  ),
  value = c(
    nrow(comparison),
    nrow(compound_level),
    cor_raw,
    cor_compound,
    sum(comparison$pvalue < 0.05, na.rm = TRUE),
    sum(comparison$destress_pvalue < 0.05, na.rm = TRUE),
    sum(compound_level$workflow_padj_by_promoter < 0.05, na.rm = TRUE),
    sum(compound_level$destress_padj_by_promoter < 0.05, na.rm = TRUE)
  )
)

write.table(
  summary_table,
  file.path(out_dir, "workflow_vs_destress_summary.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

print(summary_table)
message("Wrote comparison outputs to: ", out_dir)
