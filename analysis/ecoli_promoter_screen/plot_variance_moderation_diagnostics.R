source(file.path("analysis", "_helpers.R"))
load_destress_package()

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package `ggplot2` is required for E. coli variance diagnostics.", call. = FALSE)
}
if (!requireNamespace("scales", quietly = TRUE)) {
  stop("Package `scales` is required for E. coli variance diagnostics.", call. = FALSE)
}

ggplot2 <- asNamespace("ggplot2")

load(analysis_path("data", "binsfeld_reporter_data.rda"))

out_dir <- analysis_output_dir("binsfeld_variance_diagnostics")

wt_auc <- binsfeld_reporter_auc[
  binsfeld_reporter_auc$strain == "WT" &
    binsfeld_reporter_auc$removed == "No",
]

prepare_ecoli_assay <- function(use_ev_control) {
  assay_data <- if (isTRUE(use_ev_control)) {
    wt_auc
  } else {
    wt_auc[wt_auc$promoter != "EVC", , drop = FALSE]
  }
  args <- list(
    data = assay_data,
    promoter = "promoter",
    compound = "compound",
    control = "Water",
    lux = "lux_auc",
    growth = "od_auc",
    growth_exponent = "estimate",
    batch = "dose_level",
    replicate = "replicate",
    growth_covariates = "replicate",
    numeric_covariates = "dose_level"
  )
  if (isTRUE(use_ev_control)) {
    args$background_promoter <- "EVC"
    args$background_method <- "huber"
    args$background_by <- c("compound", "dose_level", "replicate")
  }
  do.call(prepare_assay, args)
}

fit_pair_results <- function(assay, empirical_bayes) {
  fit <- fit_destress(
    assay,
    technical = c("replicate", "dose_level"),
    empirical_bayes = empirical_bayes,
    adjustment = "by_promoter",
    interaction = FALSE
  )
  list(fit = fit, results = results(fit))
}

variance_table <- function(label, use_ev_control) {
  assay <- prepare_ecoli_assay(use_ev_control)
  raw_fit <- fit_pair_results(assay, empirical_bayes = FALSE)
  moderated_fit <- fit_pair_results(assay, empirical_bayes = TRUE)
  raw <- raw_fit$results
  moderated <- moderated_fit$results
  residual_df <- min(raw_fit$fit$promoter_effects$residual_df, na.rm = TRUE)
  prior <- estimate_eb_prior_variance(raw$specific_se^2, residual_df)

  tab <- data.frame(
    workflow = label,
    promoter = raw$promoter,
    compound = raw$compound,
    mean_response = raw$total_effect,
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
  tab$moderated_to_raw_ratio <- tab$moderated_variance / tab$raw_variance
  tab[order(tab$raw_variance, decreasing = TRUE), ]
}

diagnostics <- rbind(
  variance_table("DStressR without EV", use_ev_control = FALSE),
  variance_table("DStressR with EV", use_ev_control = TRUE)
)

diagnostics$rank_by_raw_variance <- ave(
  diagnostics$raw_variance,
  diagnostics$workflow,
  FUN = function(x) rank(-x, ties.method = "first")
)

utils::write.table(
  diagnostics,
  file.path(out_dir, "ecoli_raw_vs_moderated_variance.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

plot_data <- rbind(
  data.frame(
    workflow = diagnostics$workflow,
    rank_by_raw_variance = diagnostics$rank_by_raw_variance,
    variance_type = "Raw variance",
    variance = diagnostics$raw_variance,
    stringsAsFactors = FALSE
  ),
  data.frame(
    workflow = diagnostics$workflow,
    rank_by_raw_variance = diagnostics$rank_by_raw_variance,
    variance_type = "Moderated variance",
    variance = diagnostics$moderated_variance,
    stringsAsFactors = FALSE
  )
)

plot_data$workflow <- factor(
  plot_data$workflow,
  levels = c("DStressR without EV", "DStressR with EV")
)
plot_data$variance_type <- factor(
  plot_data$variance_type,
  levels = c("Raw variance", "Moderated variance")
)

variance_plot <- ggplot2$ggplot(
  plot_data,
  ggplot2$aes(rank_by_raw_variance, variance, color = variance_type)
) +
  ggplot2$geom_line(linewidth = 0.55, alpha = 0.95) +
  ggplot2$geom_point(size = 0.9, alpha = 0.55) +
  ggplot2$facet_wrap(ggplot2$vars(workflow), ncol = 1, scales = "free_y") +
  ggplot2$scale_y_log10(labels = scales::label_number()) +
  ggplot2$scale_color_manual(
    values = c("Raw variance" = "#64748b", "Moderated variance" = "#dc2626")
  ) +
  ggplot2$theme_light(base_size = 11) +
  ggplot2$theme(
    panel.grid.minor = ggplot2$element_blank(),
    strip.text = ggplot2$element_text(face = "bold"),
    plot.title = ggplot2$element_text(face = "bold"),
    legend.position = "top"
  ) +
  ggplot2$labs(
    title = "Raw and moderated specific-effect variances in the E. coli screen",
    subtitle = "Promoter-compound pairs are sorted from high to low raw variance within each workflow",
    x = "Promoter-compound pair index, sorted by raw variance",
    y = "Specific-effect variance (log10 scale)",
    color = NULL
  )

ggplot2$ggsave(
  file.path(out_dir, "ecoli_raw_vs_moderated_variance.png"),
  variance_plot,
  width = 8.2,
  height = 6.8,
  dpi = 320
)
ggplot2$ggsave(
  file.path(out_dir, "ecoli_raw_vs_moderated_variance.pdf"),
  variance_plot,
  width = 8.2,
  height = 6.8
)

diagnostics$workflow <- factor(
  diagnostics$workflow,
  levels = c("DStressR without EV", "DStressR with EV")
)

mean_variance_plot <- ggplot2$ggplot(
  diagnostics,
  ggplot2$aes(mean_response, raw_variance)
) +
  ggplot2$geom_point(
    ggplot2$aes(y = moderated_variance),
    color = "#2563eb",
    size = 1.15,
    alpha = 0.62
  ) +
  ggplot2$geom_point(color = "#111827", shape = 1, size = 1.45, alpha = 0.62, stroke = 0.35) +
  ggplot2$geom_smooth(
    color = "#dc2626",
    method = "loess",
    formula = y ~ x,
    se = FALSE,
    linewidth = 0.9
  ) +
  ggplot2$facet_wrap(ggplot2$vars(workflow), ncol = 1) +
  ggplot2$scale_y_log10(labels = scales::label_number()) +
  ggplot2$theme_light(base_size = 11) +
  ggplot2$theme(
    panel.grid.minor = ggplot2$element_blank(),
    strip.text = ggplot2$element_text(face = "bold"),
    plot.title = ggplot2$element_text(face = "bold"),
    legend.position = "none"
  ) +
  ggplot2$labs(
    title = "Response-contrast variance diagnostic for the E. coli screen",
    subtitle = "Black open circles: raw pairwise variance; red: contrast-variance trend; blue: moderated variance",
    x = "Estimated promoter-compound response contrast",
    y = "Specific-effect variance (log10 scale)"
  )

ggplot2$ggsave(
  file.path(out_dir, "ecoli_mean_response_vs_variance.png"),
  mean_variance_plot,
  width = 8.2,
  height = 6.8,
  dpi = 320
)
ggplot2$ggsave(
  file.path(out_dir, "ecoli_mean_response_vs_variance.pdf"),
  mean_variance_plot,
  width = 8.2,
  height = 6.8
)

ratio_plot <- ggplot2$ggplot(
  diagnostics,
  ggplot2$aes(rank_by_raw_variance, moderated_to_raw_ratio)
) +
  ggplot2$geom_hline(yintercept = 1, color = "#111827", linewidth = 0.45) +
  ggplot2$geom_point(color = "#2563eb", size = 1.25, alpha = 0.72) +
  ggplot2$facet_wrap(ggplot2$vars(workflow), ncol = 1) +
  ggplot2$theme_light(base_size = 11) +
  ggplot2$theme(
    panel.grid.minor = ggplot2$element_blank(),
    strip.text = ggplot2$element_text(face = "bold"),
    plot.title = ggplot2$element_text(face = "bold")
  ) +
  ggplot2$labs(
    title = "Magnitude of variance moderation in the E. coli screen",
    subtitle = "Promoter-compound pairs are sorted by raw variance; the horizontal line marks no moderation",
    x = "Promoter-compound pair index, sorted by raw variance",
    y = "Moderated / raw variance"
  )

ggplot2$ggsave(
  file.path(out_dir, "ecoli_moderated_to_raw_variance_ratio.png"),
  ratio_plot,
  width = 8.2,
  height = 5.8,
  dpi = 320
)
ggplot2$ggsave(
  file.path(out_dir, "ecoli_moderated_to_raw_variance_ratio.pdf"),
  ratio_plot,
  width = 8.2,
  height = 5.8
)

summary_table <- aggregate(
  cbind(raw_variance, moderated_variance, moderated_to_raw_ratio) ~ workflow,
  diagnostics,
  function(x) c(
    min = min(x, na.rm = TRUE),
    median = stats::median(x, na.rm = TRUE),
    max = max(x, na.rm = TRUE)
  )
)
prior_summary <- do.call(rbind, lapply(split(diagnostics, diagnostics$workflow), function(d) {
  data.frame(
    workflow = unique(d$workflow),
    residual_df = unique(d$residual_df),
    diagnostic_prior_variance = unique(d$diagnostic_prior_variance),
    diagnostic_prior_df = unique(d$diagnostic_prior_df),
    stringsAsFactors = FALSE
  )
}))
summary_table <- merge(summary_table, prior_summary, by = "workflow", all.x = TRUE, sort = FALSE)
utils::write.table(
  summary_table,
  file.path(out_dir, "ecoli_raw_vs_moderated_variance_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Wrote E. coli variance moderation diagnostics to: ", out_dir)
