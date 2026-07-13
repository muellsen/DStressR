#!/usr/bin/env Rscript

source(file.path("analysis", "_helpers.R"))

suppressPackageStartupMessages({
  library(ggplot2)
})

root <- analysis_data_root()
out_dir <- file.path(getwd(), "analysis", "outputs", "student_t")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

comparison_file <- file.path(getwd(), "analysis", "outputs", "workflow_vs_destress_replicate_pvalues.tsv")
libmap_file <- file.path(root, "00-import", "Campylobacter", "LibMap.txt")

read_tsv_base <- function(path) {
  read.delim(path, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
}

safe_log10p <- function(p) {
  -log10(pmax(p, .Machine$double.xmin))
}

fit_t_tail <- function(z) {
  z <- z[is.finite(z)]
  nll <- function(par) {
    df <- 2 + exp(par[1])
    scale <- exp(par[2])
    -sum(stats::dt(z / scale, df = df, log = TRUE) - log(scale))
  }
  fit <- stats::optim(
    par = c(log(5 - 2), log(1)),
    fn = nll,
    method = "Nelder-Mead",
    control = list(maxit = 5000)
  )
  c(df = 2 + exp(fit$par[1]), scale = exp(fit$par[2]), nll = fit$value, n = length(z))
}

fit_t_df_scale1 <- function(z) {
  z <- z[is.finite(z)]
  nll <- function(par) {
    df <- 2 + exp(par[1])
    -sum(stats::dt(z, df = df, log = TRUE))
  }
  fit <- stats::optim(
    par = log(5 - 2),
    fn = nll,
    method = "Brent",
    lower = log(2.01 - 2),
    upper = log(1e6 - 2)
  )
  c(df = 2 + exp(fit$par[1]), scale = 1, nll = fit$value, n = length(z))
}

if (!file.exists(comparison_file)) {
  stop("Missing comparison table. Run analysis/compare_workflow_pvalues.R first.", call. = FALSE)
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
comparison <- comparison[is.finite(comparison$destress_z), ]

dmso_z <- comparison$destress_z[comparison$srn_code %in% dmso_srn_codes]
t_fit <- fit_t_df_scale1(dmso_z)
t_fit_free_scale <- fit_t_tail(dmso_z)
all_residual_fit <- fit_t_tail(comparison$destress_z)

comparison$destress_t_statistic <- comparison$destress_z / unname(t_fit["scale"])
comparison$destress_t_pvalue <- 2 * stats::pt(
  abs(comparison$destress_t_statistic),
  df = unname(t_fit["df"]),
  lower.tail = FALSE
)
comparison$destress_t_all_residual_statistic <- comparison$destress_z / unname(all_residual_fit["scale"])
comparison$destress_t_all_residual_pvalue <- 2 * stats::pt(
  abs(comparison$destress_t_all_residual_statistic),
  df = unname(all_residual_fit["df"]),
  lower.tail = FALSE
)
comparison$destress_t_neglog10p <- safe_log10p(comparison$destress_t_pvalue)
comparison$destress_t_all_residual_neglog10p <- safe_log10p(comparison$destress_t_all_residual_pvalue)

write.table(
  comparison,
  file.path(out_dir, "workflow_vs_destress_student_t_replicate_pvalues.tsv"),
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
compound_level$destress_t_padj_by_promoter <- ave(
  compound_level$destress_t_pvalue,
  compound_level$promoter,
  FUN = function(x) p.adjust(x, method = "BH")
)
compound_level$destress_t_all_residual_padj_by_promoter <- ave(
  compound_level$destress_t_all_residual_pvalue,
  compound_level$promoter,
  FUN = function(x) p.adjust(x, method = "BH")
)
compound_level$workflow_neglog10padj <- safe_log10p(compound_level$workflow_padj_by_promoter)
compound_level$destress_gaussian_neglog10padj <- safe_log10p(compound_level$destress_gaussian_padj_by_promoter)
compound_level$destress_t_neglog10padj <- safe_log10p(compound_level$destress_t_padj_by_promoter)
compound_level$destress_t_all_residual_neglog10padj <- safe_log10p(
  compound_level$destress_t_all_residual_padj_by_promoter
)

write.table(
  compound_level,
  file.path(out_dir, "workflow_vs_destress_student_t_promoter_compound_pvalues.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

rho_workflow_t <- stats::cor(
  safe_log10p(compound_level$pvalue),
  compound_level$destress_t_neglog10p,
  method = "spearman"
)
rho_gaussian_t <- stats::cor(
  safe_log10p(compound_level$destress_pvalue),
  compound_level$destress_t_neglog10p,
  method = "spearman"
)

summary_table <- data.frame(
  metric = c(
    "dmso_z_rows_used_for_t_fit",
    "student_t_df",
    "student_t_scale_multiplier",
    "student_t_free_scale_df",
    "student_t_free_scale_multiplier",
    "student_t_all_residual_df_sensitivity",
    "student_t_all_residual_scale_sensitivity",
    "promoter_compound_rows",
    "workflow_raw_p_lt_0.05",
    "destress_gaussian_raw_p_lt_0.05",
    "destress_student_t_raw_p_lt_0.05",
    "destress_student_t_all_residual_raw_p_lt_0.05",
    "workflow_bh_lt_0.05",
    "destress_gaussian_bh_lt_0.05",
    "destress_student_t_bh_lt_0.05",
    "destress_student_t_all_residual_bh_lt_0.05",
    "workflow_vs_student_t_spearman_neglog10p",
    "gaussian_vs_student_t_spearman_neglog10p"
  ),
  value = c(
    unname(t_fit["n"]),
    unname(t_fit["df"]),
    unname(t_fit["scale"]),
    unname(t_fit_free_scale["df"]),
    unname(t_fit_free_scale["scale"]),
    unname(all_residual_fit["df"]),
    unname(all_residual_fit["scale"]),
    nrow(compound_level),
    sum(compound_level$pvalue < 0.05, na.rm = TRUE),
    sum(compound_level$destress_pvalue < 0.05, na.rm = TRUE),
    sum(compound_level$destress_t_pvalue < 0.05, na.rm = TRUE),
    sum(compound_level$destress_t_all_residual_pvalue < 0.05, na.rm = TRUE),
    sum(compound_level$workflow_padj_by_promoter < 0.05, na.rm = TRUE),
    sum(compound_level$destress_gaussian_padj_by_promoter < 0.05, na.rm = TRUE),
    sum(compound_level$destress_t_padj_by_promoter < 0.05, na.rm = TRUE),
    sum(compound_level$destress_t_all_residual_padj_by_promoter < 0.05, na.rm = TRUE),
    rho_workflow_t,
    rho_gaussian_t
  )
)

write.table(
  summary_table,
  file.path(out_dir, "student_t_pvalue_summary.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

scatter_workflow_t <- ggplot(
  compound_level,
  aes(safe_log10p(pvalue), destress_t_neglog10p)
) +
  geom_point(alpha = 0.12, size = 0.45, color = "#1f2937") +
  geom_abline(slope = 1, intercept = 0, color = "#b91c1c", linewidth = 0.5) +
  coord_equal() +
  theme_bw(base_size = 10) +
  labs(
    title = "Workflow vs DStressR Student-t p-values",
    subtitle = paste0(
      "Student-t df = ", round(unname(t_fit["df"]), 2),
      ", scale = ", round(unname(t_fit["scale"]), 3),
      "; Spearman rho = ", round(rho_workflow_t, 3)
    ),
    x = "Workflow median-polish -log10(p)",
    y = "DStressR Student-t -log10(p)"
  )

ggsave(
  file.path(out_dir, "scatter_workflow_vs_destress_student_t_raw_pvalues.png"),
  scatter_workflow_t,
  width = 6.5,
  height = 6,
  dpi = 300
)
ggsave(
  file.path(out_dir, "scatter_workflow_vs_destress_student_t_raw_pvalues.pdf"),
  scatter_workflow_t,
  width = 6.5,
  height = 6
)

scatter_gaussian_t <- ggplot(
  compound_level,
  aes(safe_log10p(destress_pvalue), destress_t_neglog10p)
) +
  geom_point(alpha = 0.12, size = 0.45, color = "#1f2937") +
  geom_abline(slope = 1, intercept = 0, color = "#b91c1c", linewidth = 0.5) +
  coord_equal() +
  theme_bw(base_size = 10) +
  labs(
    title = "DStressR Gaussian vs Student-t p-values",
    subtitle = paste0("Spearman rho = ", round(rho_gaussian_t, 3)),
    x = "DStressR Gaussian -log10(p)",
    y = "DStressR Student-t -log10(p)"
  )

ggsave(
  file.path(out_dir, "scatter_destress_gaussian_vs_student_t_raw_pvalues.png"),
  scatter_gaussian_t,
  width = 6.5,
  height = 6,
  dpi = 300
)
ggsave(
  file.path(out_dir, "scatter_destress_gaussian_vs_student_t_raw_pvalues.pdf"),
  scatter_gaussian_t,
  width = 6.5,
  height = 6
)

scatter_adjusted <- ggplot(
  compound_level,
  aes(destress_gaussian_neglog10padj, destress_t_neglog10padj)
) +
  geom_point(alpha = 0.15, size = 0.5, color = "#1f2937") +
  geom_abline(slope = 1, intercept = 0, color = "#b91c1c", linewidth = 0.5) +
  coord_equal() +
  theme_bw(base_size = 10) +
  labs(
    title = "DStressR Gaussian vs Student-t adjusted p-values",
    x = "DStressR Gaussian -log10(BH adjusted p)",
    y = "DStressR Student-t -log10(BH adjusted p)"
  )

ggsave(
  file.path(out_dir, "scatter_destress_gaussian_vs_student_t_adjusted_pvalues.png"),
  scatter_adjusted,
  width = 6.5,
  height = 6,
  dpi = 300
)
ggsave(
  file.path(out_dir, "scatter_destress_gaussian_vs_student_t_adjusted_pvalues.pdf"),
  scatter_adjusted,
  width = 6.5,
  height = 6
)

p_hist <- ggplot(data.frame(dmso_z = dmso_z), aes(dmso_z)) +
  geom_histogram(aes(y = after_stat(density)), bins = 90, fill = "#334155", color = "white") +
  stat_function(fun = dnorm, color = "#2563eb", linewidth = 0.7) +
  stat_function(
    fun = function(x) stats::dt(x / unname(t_fit["scale"]), df = unname(t_fit["df"])) / unname(t_fit["scale"]),
    color = "#b91c1c",
    linewidth = 0.7
  ) +
  coord_cartesian(xlim = c(-6, 6)) +
  theme_bw(base_size = 10) +
  labs(
    title = "DMSO standardized residual tails",
    subtitle = "Blue: standard normal; red: fitted Student-t",
    x = "DStressR DMSO z residual",
    y = "Density"
  )

ggsave(
  file.path(out_dir, "dmso_z_student_t_fit.png"),
  p_hist,
  width = 7,
  height = 5,
  dpi = 300
)
ggsave(
  file.path(out_dir, "dmso_z_student_t_fit.pdf"),
  p_hist,
  width = 7,
  height = 5
)

print(summary_table)
message("Wrote Student-t outputs to: ", out_dir)
