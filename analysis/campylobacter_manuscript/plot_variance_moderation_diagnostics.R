source(file.path("analysis", "_helpers.R"))
load_destress_package()

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package `ggplot2` is required for Campylobacter variance diagnostics.", call. = FALSE)
}
if (!requireNamespace("scales", quietly = TRUE)) {
  stop("Package `scales` is required for Campylobacter variance diagnostics.", call. = FALSE)
}

ggplot2 <- asNamespace("ggplot2")

out_dir <- analysis_output_dir("campylobacter_variance_diagnostics")

data_root <- analysis_data_root()
expr_path <- file.path(data_root, "03-hit_determination", "expression_table.tsv.gz")
if (!file.exists(expr_path)) {
  stop("Campylobacter expression table not found: ", expr_path, call. = FALSE)
}

expr <- read_tsv_base(expr_path)
expr <- expr[!(expr$promoter %in% c("PCJnc20", "PCjas704")), , drop = FALSE]
expr$libplate <- sub("_.*$", "", expr$srn_code)
expr$compound_model <- ifelse(expr$ProductName == "DMSO", "DMSO", expr$srn_code)

assay <- prepare_assay(
  expr,
  promoter = "promoter",
  compound = "compound_model",
  control = "DMSO",
  lux = "lux_auc_until16h",
  growth = "od_at_16h",
  growth_exponent = "estimate",
  plate = "libplate",
  replicate = "replicate"
)

fit_pair_results <- function(empirical_bayes) {
  fit <- fit_destress(
    assay,
    technical = c("libplate", "replicate"),
    empirical_bayes = empirical_bayes,
    interaction = FALSE,
    adjustment = "global",
    background_rank = 0
  )
  list(fit = fit, results = results(fit))
}

raw_fit <- fit_pair_results(empirical_bayes = FALSE)
moderated_fit <- fit_pair_results(empirical_bayes = TRUE)
raw <- raw_fit$results
moderated <- moderated_fit$results
residual_df <- min(raw_fit$fit$promoter_effects$residual_df, na.rm = TRUE)
prior <- estimate_eb_prior_variance(raw$specific_se^2, residual_df)

diagnostics <- data.frame(
  workflow = "DStressR Campylobacter workflow",
  promoter = raw$promoter,
  compound = raw$compound,
  response_contrast = raw$total_effect,
  effect = raw$specific_effect,
  raw_variance = raw$specific_se^2,
  moderated_variance = moderated$specific_se^2,
  raw_pvalue = raw$specific_pvalue,
  moderated_pvalue = moderated$specific_pvalue,
  residual_df = residual_df,
  diagnostic_prior_variance = prior$var,
  diagnostic_prior_df = prior$df,
  stringsAsFactors = FALSE
)
diagnostics$moderated_to_raw_ratio <- diagnostics$moderated_variance / diagnostics$raw_variance
diagnostics <- diagnostics[order(diagnostics$raw_variance, decreasing = TRUE), , drop = FALSE]
diagnostics$rank_by_raw_variance <- seq_len(nrow(diagnostics))

utils::write.table(
  diagnostics,
  file.path(out_dir, "campylobacter_raw_vs_moderated_variance.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

variance_plot_data <- rbind(
  data.frame(
    rank_by_raw_variance = diagnostics$rank_by_raw_variance,
    variance_type = "Raw variance",
    variance = diagnostics$raw_variance,
    stringsAsFactors = FALSE
  ),
  data.frame(
    rank_by_raw_variance = diagnostics$rank_by_raw_variance,
    variance_type = "Moderated variance",
    variance = diagnostics$moderated_variance,
    stringsAsFactors = FALSE
  )
)
variance_plot_data$variance_type <- factor(
  variance_plot_data$variance_type,
  levels = c("Raw variance", "Moderated variance")
)

sorted_variance_plot <- ggplot2$ggplot(
  variance_plot_data,
  ggplot2$aes(rank_by_raw_variance, variance, color = variance_type)
) +
  ggplot2$geom_line(linewidth = 0.55, alpha = 0.95) +
  ggplot2$geom_point(size = 0.8, alpha = 0.42) +
  ggplot2$scale_y_log10(labels = scales::label_number()) +
  ggplot2$scale_color_manual(
    values = c("Raw variance" = "#64748b", "Moderated variance" = "#dc2626")
  ) +
  ggplot2$theme_light(base_size = 11) +
  ggplot2$theme(
    panel.grid.minor = ggplot2$element_blank(),
    plot.title = ggplot2$element_text(face = "bold"),
    legend.position = "top"
  ) +
  ggplot2$labs(
    title = "Raw and moderated specific-effect variances in the Campylobacter screen",
    subtitle = "Promoter-compound pairs are sorted from high to low raw variance",
    x = "Promoter-compound pair index, sorted by raw variance",
    y = "Specific-effect variance (log10 scale)",
    color = NULL
  )

ggplot2$ggsave(
  file.path(out_dir, "campylobacter_raw_vs_moderated_variance.png"),
  sorted_variance_plot,
  width = 8.2,
  height = 4.8,
  dpi = 320
)
ggplot2$ggsave(
  file.path(out_dir, "campylobacter_raw_vs_moderated_variance.pdf"),
  sorted_variance_plot,
  width = 8.2,
  height = 4.8
)

contrast_variance_plot <- ggplot2$ggplot(
  diagnostics,
  ggplot2$aes(response_contrast, raw_variance)
) +
  ggplot2$geom_point(
    ggplot2$aes(y = moderated_variance),
    color = "#2563eb",
    size = 1.0,
    alpha = 0.52
  ) +
  ggplot2$geom_point(color = "#111827", shape = 1, size = 1.25, alpha = 0.55, stroke = 0.32) +
  ggplot2$geom_smooth(
    color = "#dc2626",
    method = "loess",
    formula = y ~ x,
    se = FALSE,
    linewidth = 0.9
  ) +
  ggplot2$scale_y_log10(labels = scales::label_number()) +
  ggplot2$theme_light(base_size = 11) +
  ggplot2$theme(
    panel.grid.minor = ggplot2$element_blank(),
    plot.title = ggplot2$element_text(face = "bold")
  ) +
  ggplot2$labs(
    title = "Response-contrast variance diagnostic for the Campylobacter screen",
    subtitle = "Black open circles: raw pairwise variance; red: contrast-variance trend; blue: moderated variance",
    x = "Estimated promoter-compound response contrast",
    y = "Specific-effect variance (log10 scale)"
  )

ggplot2$ggsave(
  file.path(out_dir, "campylobacter_response_contrast_vs_variance.png"),
  contrast_variance_plot,
  width = 8.2,
  height = 4.8,
  dpi = 320
)
ggplot2$ggsave(
  file.path(out_dir, "campylobacter_response_contrast_vs_variance.pdf"),
  contrast_variance_plot,
  width = 8.2,
  height = 4.8
)

ratio_plot <- ggplot2$ggplot(
  diagnostics,
  ggplot2$aes(rank_by_raw_variance, moderated_to_raw_ratio)
) +
  ggplot2$geom_hline(yintercept = 1, color = "#111827", linewidth = 0.45) +
  ggplot2$geom_point(color = "#2563eb", size = 1.0, alpha = 0.62) +
  ggplot2$theme_light(base_size = 11) +
  ggplot2$theme(
    panel.grid.minor = ggplot2$element_blank(),
    plot.title = ggplot2$element_text(face = "bold")
  ) +
  ggplot2$labs(
    title = "Magnitude of variance moderation in the Campylobacter screen",
    subtitle = "Promoter-compound pairs are sorted by raw variance; the horizontal line marks no moderation",
    x = "Promoter-compound pair index, sorted by raw variance",
    y = "Moderated / raw variance"
  )

ggplot2$ggsave(
  file.path(out_dir, "campylobacter_moderated_to_raw_variance_ratio.png"),
  ratio_plot,
  width = 8.2,
  height = 4.8,
  dpi = 320
)
ggplot2$ggsave(
  file.path(out_dir, "campylobacter_moderated_to_raw_variance_ratio.pdf"),
  ratio_plot,
  width = 8.2,
  height = 4.8
)

summary_table <- data.frame(
  raw_variance_min = min(diagnostics$raw_variance, na.rm = TRUE),
  raw_variance_median = stats::median(diagnostics$raw_variance, na.rm = TRUE),
  raw_variance_max = max(diagnostics$raw_variance, na.rm = TRUE),
  moderated_variance_min = min(diagnostics$moderated_variance, na.rm = TRUE),
  moderated_variance_median = stats::median(diagnostics$moderated_variance, na.rm = TRUE),
  moderated_variance_max = max(diagnostics$moderated_variance, na.rm = TRUE),
  moderated_to_raw_ratio_min = min(diagnostics$moderated_to_raw_ratio, na.rm = TRUE),
  moderated_to_raw_ratio_median = stats::median(diagnostics$moderated_to_raw_ratio, na.rm = TRUE),
  moderated_to_raw_ratio_max = max(diagnostics$moderated_to_raw_ratio, na.rm = TRUE),
  residual_df = unique(diagnostics$residual_df),
  diagnostic_prior_variance = unique(diagnostics$diagnostic_prior_variance),
  diagnostic_prior_df = unique(diagnostics$diagnostic_prior_df)
)
utils::write.table(
  summary_table,
  file.path(out_dir, "campylobacter_raw_vs_moderated_variance_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Wrote Campylobacter variance moderation diagnostics to: ", out_dir)
