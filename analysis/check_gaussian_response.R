#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

root <- "/Users/cmueller/Documents/GitHub/campylobacter_stressregnet/workflow/data"
out_dir <- file.path(getwd(), "analysis", "outputs", "gaussian_diagnostics")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

expression_file <- file.path(root, "02-lux_expression", "expression_values.tsv.gz")
libmap_file <- file.path(root, "00-import", "Campylobacter", "LibMap.txt")

read_tsv_base <- function(path) {
  read.delim(path, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
}

skewness <- function(x) {
  x <- x[is.finite(x)]
  m <- mean(x)
  s <- stats::sd(x)
  mean(((x - m) / s)^3)
}

excess_kurtosis <- function(x) {
  x <- x[is.finite(x)]
  m <- mean(x)
  s <- stats::sd(x)
  mean(((x - m) / s)^4) - 3
}

sample_for_plot <- function(x, n = 50000, seed = 11) {
  if (nrow(x) <= n) {
    return(x)
  }
  set.seed(seed)
  x[sample.int(nrow(x), n), , drop = FALSE]
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

expr <- read_tsv_base(expression_file)
expr <- merge(
  expr,
  libmap[, c("srn_code", "ProductName", "Catalog Number")],
  by = "srn_code",
  all.x = TRUE,
  sort = FALSE
)
expr <- expr[
  expr[["Catalog Number"]] != "DMSO noisy" &
    !(expr$promoter %in% c("PCJnc20", "PCjas704")) &
    is.finite(expr$LUX.AUC_16) &
    is.finite(expr$od_16h.measured) &
    expr$LUX.AUC_16 > 0 &
    expr$od_16h.measured > 0,
]

expr$Y_raw_lux <- expr$LUX.AUC_16
expr$Y_log_lux <- log2(expr$LUX.AUC_16)
expr$Y_log_lux_over_growth <- log2(expr$LUX.AUC_16 / expr$od_16h.measured)
expr$Y_workflow <- expr$log2.auc.16hmeasured.normed

expr <- expr[order(expr$promoter, expr$srn_code), ]
expr$workflow_replicate <- rep(c("r1", "r2"), length.out = nrow(expr))
expr$promoter_libplate_replicate <- paste(
  expr$promoter,
  expr$libplate,
  expr$workflow_replicate,
  sep = "_"
)

dmso_expr <- expr[expr$srn_code %in% dmso_srn_codes, ]
dmso_mean <- aggregate(
  Y_workflow ~ promoter_libplate_replicate,
  dmso_expr,
  mean
)
names(dmso_mean)[2] <- "dmso_mean_y"

diag_df <- merge(expr, dmso_mean, by = "promoter_libplate_replicate", all.x = TRUE, sort = FALSE)
diag_df$Y_dmso_centered <- diag_df$Y_workflow - diag_df$dmso_mean_y

compound_global <- aggregate(Y_dmso_centered ~ srn_code, diag_df, mean, na.rm = TRUE)
names(compound_global)[2] <- "compound_global"
diag_df <- merge(diag_df, compound_global, by = "srn_code", all.x = TRUE, sort = FALSE)
diag_df$Y_specific_residual <- diag_df$Y_dmso_centered - diag_df$compound_global

response_long <- rbind(
  data.frame(response = "raw_lux_auc", value = diag_df$Y_raw_lux),
  data.frame(response = "log2_lux_auc", value = diag_df$Y_log_lux),
  data.frame(response = "log2_lux_over_od", value = diag_df$Y_log_lux_over_growth),
  data.frame(response = "dmso_centered_log2_lux_over_od", value = diag_df$Y_dmso_centered),
  data.frame(response = "compound_global_adjusted_residual", value = diag_df$Y_specific_residual)
)
response_long <- response_long[is.finite(response_long$value), ]

summary_stats <- do.call(
  rbind,
  lapply(split(response_long$value, response_long$response), function(x) {
    data.frame(
      n = length(x),
      mean = mean(x),
      sd = stats::sd(x),
      median = stats::median(x),
      mad = stats::mad(x),
      skewness = skewness(x),
      excess_kurtosis = excess_kurtosis(x),
      q01 = unname(stats::quantile(x, 0.01)),
      q99 = unname(stats::quantile(x, 0.99))
    )
  })
)
summary_stats$response <- rownames(summary_stats)
summary_stats <- summary_stats[, c("response", setdiff(names(summary_stats), "response"))]
rownames(summary_stats) <- NULL

write.table(
  summary_stats,
  file.path(out_dir, "response_distribution_summary.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

p_hist <- ggplot(sample_for_plot(response_long), aes(value)) +
  geom_histogram(bins = 80, color = "white", fill = "#334155") +
  facet_wrap(~response, scales = "free", ncol = 1) +
  theme_bw(base_size = 10) +
  labs(
    title = "Response distributions",
    x = "Value",
    y = "Count"
  )
ggsave(file.path(out_dir, "response_histograms.png"), p_hist, width = 7, height = 10, dpi = 300)
ggsave(file.path(out_dir, "response_histograms.pdf"), p_hist, width = 7, height = 10)

qq_responses <- response_long[
  response_long$response %in% c(
    "log2_lux_over_od",
    "dmso_centered_log2_lux_over_od",
    "compound_global_adjusted_residual"
  ),
]
p_qq <- ggplot(sample_for_plot(qq_responses, n = 25000), aes(sample = value)) +
  stat_qq(alpha = 0.18, size = 0.35, color = "#1f2937") +
  stat_qq_line(color = "#b91c1c", linewidth = 0.5) +
  facet_wrap(~response, scales = "free", ncol = 3) +
  theme_bw(base_size = 10) +
  labs(
    title = "Normal Q-Q diagnostics for transformed responses",
    x = "Theoretical normal quantile",
    y = "Observed quantile"
  )
ggsave(file.path(out_dir, "qq_transformed_responses.png"), p_qq, width = 12, height = 4, dpi = 300)
ggsave(file.path(out_dir, "qq_transformed_responses.pdf"), p_qq, width = 12, height = 4)

mean_var <- aggregate(
  cbind(mean_lux = LUX.AUC_16, var_lux = LUX.AUC_16, mean_log = Y_log_lux, var_log = Y_log_lux) ~ promoter,
  diag_df,
  function(x) c(mean = mean(x), var = stats::var(x))
)

# Base aggregate with matrix columns is awkward; compute explicitly.
mv <- do.call(
  rbind,
  lapply(split(diag_df, diag_df$promoter), function(d) {
    data.frame(
      promoter = d$promoter[1],
      mean_lux = mean(d$LUX.AUC_16, na.rm = TRUE),
      var_lux = stats::var(d$LUX.AUC_16, na.rm = TRUE),
      mean_log2_lux_over_od = mean(d$Y_workflow, na.rm = TRUE),
      var_log2_lux_over_od = stats::var(d$Y_workflow, na.rm = TRUE)
    )
  })
)

write.table(
  mv,
  file.path(out_dir, "mean_variance_by_promoter.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

p_mv_raw <- ggplot(mv, aes(mean_lux, var_lux, label = promoter)) +
  geom_point(color = "#1f2937") +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw(base_size = 10) +
  labs(
    title = "Raw LUX AUC mean-variance relation by promoter",
    x = "Mean raw LUX AUC",
    y = "Variance raw LUX AUC"
  )
ggsave(file.path(out_dir, "mean_variance_raw_lux.png"), p_mv_raw, width = 6, height = 5, dpi = 300)

p_mv_log <- ggplot(mv, aes(mean_log2_lux_over_od, var_log2_lux_over_od, label = promoter)) +
  geom_point(color = "#1f2937") +
  theme_bw(base_size = 10) +
  labs(
    title = "Log2(LUX / OD) mean-variance relation by promoter",
    x = "Mean log2(LUX / OD)",
    y = "Variance log2(LUX / OD)"
  )
ggsave(file.path(out_dir, "mean_variance_log_response.png"), p_mv_log, width = 6, height = 5, dpi = 300)

message("Wrote Gaussian diagnostics to: ", out_dir)
print(summary_stats)
