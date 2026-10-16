# Compare a preserved submission baseline with revised inference on identical fits.
# Run from the repository root after creating the baseline source checkpoint.
source(file.path("analysis", "_helpers.R"))
load_binsfeld_paper_data()
out_dir <- analysis_output_dir("inference_audit_2026-10-16")
baseline_dir <- file.path(out_dir, "baseline_R")
if (!dir.exists(baseline_dir)) stop("The prechange baseline_R checkpoint is required.")

load_sources <- function(path) {
  env <- new.env(parent = .GlobalEnv)
  for (f in list.files(path, pattern = "[.]R$", full.names = TRUE)) sys.source(f, env)
  env
}
old <- load_sources(baseline_dir)
candidate <- load_sources(analysis_path("R"))
focal <- load_sources(analysis_path("R"))
focal$centered_contrast_df <- function(variance, df) df
write_tsv <- function(x, name) utils::write.table(
  x, file.path(out_dir, name), sep = "\t", quote = FALSE, row.names = FALSE
)

wt <- binsfeld_reporter_auc[
  binsfeld_reporter_auc$strain == "WT" & binsfeld_reporter_auc$removed == "No", ]
fits_path <- file.path(out_dir, "baseline_fits.rds")
if (file.exists(fits_path)) {
  fits <- readRDS(fits_path)
} else {
  fits <- lapply(c("without_EVC", "with_EVC"), function(workflow) {
    d <- if (workflow == "without_EVC") wt[wt$promoter != "EVC", ] else wt
    args <- list(data = d, reporter = "promoter", perturbation = "compound",
                 control = "Water", lux = "lux_auc", growth = "od_auc",
                 growth_exponent = "estimate", batch = "dose_level",
                 replicate = "replicate", growth_covariates = "replicate",
                 numeric_covariates = "dose_level")
    if (workflow == "with_EVC") {
      args$background_reporter <- "EVC"
      args$background_method <- "huber"
      args$background_by <- c("compound", "dose_level", "replicate")
    }
    assay <- do.call(old$prepare_assay, args)
    old$fit_destress(assay, technical = c("replicate", "dose_level"),
                     empirical_bayes = TRUE, adjustment = "by_reporter")
  })
  names(fits) <- c("without_EVC", "with_EVC")
  saveRDS(fits, fits_path)
}

tables <- list()
summary_rows <- list()
changes <- list()
for (workflow in names(fits)) {
  fit <- fits[[workflow]]
  ordinary <- fit
  ordinary$empirical_bayes <- FALSE
  runs <- list(
    submission_baseline = old$results(fit),
    ordinary_focal = old$results(ordinary),
    ordinary_WS = candidate$results(ordinary),
    single_focal = focal$results(fit),
    single_WS = candidate$results(fit)
  )
  baseline <- runs$submission_baseline
  baseline$key <- paste(baseline$reporter, baseline$perturbation, sep = "::")
  for (method in names(runs)) {
    x <- runs[[method]]
    x$key <- paste(x$reporter, x$perturbation, sep = "::")
    x <- x[match(baseline$key, x$key), ]
    stopifnot(!anyNA(x$key), identical(x$specific_effect, baseline$specific_effect))
    if (!"residual_df" %in% names(x)) {
      raw_key <- paste(fit$reporter_effects$reporter, fit$reporter_effects$perturbation, sep = "::")
      x$residual_df <- fit$reporter_effects$residual_df[match(x$key, raw_key)]
    }
    x$hit <- x$specific_padj_by_reporter < 0.05
    baseline_hit <- baseline$specific_padj_by_reporter < 0.05
    x$workflow <- workflow
    x$method <- method
    x$baseline_pvalue <- baseline$specific_pvalue
    x$baseline_padj <- baseline$specific_padj_by_reporter
    x$baseline_hit <- baseline_hit
    x$change <- ifelse(x$hit & !baseline_hit, "Gained",
                       ifelse(!x$hit & baseline_hit, "Lost", "Unchanged"))
    tables[[paste(workflow, method)]] <- x
    if (method == "single_WS") changes[[workflow]] <- x[x$change != "Unchanged", ]
    raw_df <- if ("specific_df_raw" %in% names(x)) x$specific_df_raw else x$residual_df
    fit_df <- if ("specific_df_for_test" %in% names(x)) x$specific_df_for_test else NA_real_
    prior <- if (grepl("single", method)) candidate$eb_moderate_se(sqrt(x$specific_var), raw_df) else NULL
    summary_rows[[paste(workflow, method)]] <- data.frame(
      workflow = workflow, method = method, n_pairs = nrow(x), hits = sum(x$hit),
      gained = sum(x$hit & !baseline_hit), lost = sum(!x$hit & baseline_hit),
      max_abs_effect_difference = max(abs(x$specific_effect - baseline$specific_effect)),
      median_abs_log10_p_change = median(abs(log10(pmax(x$specific_pvalue, 1e-300)) -
                                              log10(pmax(baseline$specific_pvalue, 1e-300)))),
      raw_df_min = min(raw_df), raw_df_median = median(raw_df), raw_df_max = max(raw_df),
      prior_df = if (is.null(prior)) NA_real_ else prior$prior_df,
      prior_variance = if (is.null(prior)) NA_real_ else prior$prior_var,
      test_df_median = median(fit_df)
    )
    write_tsv(x, paste0(workflow, "_", method, ".tsv"))
  }
}
summary_table <- do.call(rbind, summary_rows)
write_tsv(summary_table, "inference_comparison_summary.tsv")
all_changes <- do.call(rbind, changes)
write_tsv(all_changes, "candidate_changed_pairs.tsv")
saveRDS(tables, file.path(out_dir, "all_result_tables.rds"))

# Verify that the baseline really recovers the submission tables before drawing conclusions.
archive <- file.path(out_dir, "prechange_submission_snapshot.tar.gz")
published_path <- file.path("analysis", "outputs", "binsfeld_evc_calibrated",
                            "evc_huber_comparison_to_binsfeld_and_default.tsv")
verification_dir <- tempfile("dstressr-frozen-check-")
dir.create(verification_dir)
utils::untar(archive, files = published_path, exdir = verification_dir)
published <- read.delim(file.path(verification_dir, published_path))
checks <- lapply(names(fits), function(workflow) {
  x <- tables[[paste(workflow, "submission_baseline")]]
  ix <- match(x$key, paste(published$promoter, published$compound, sep = "::"))
  prefix <- if (workflow == "without_EVC") "modeled" else "evc_huber"
  p <- published[[paste0(prefix, "_pvalue")]][ix]
  q <- published[[paste0(prefix, "_padj_by_reporter")]][ix]
  data.frame(workflow = workflow, matched = sum(!is.na(ix)),
             max_abs_p_difference = max(abs(x$specific_pvalue - p), na.rm = TRUE),
             max_abs_q_difference = max(abs(x$specific_padj_by_reporter - q), na.rm = TRUE),
             hit_mismatches = sum(x$hit != (q < 0.05), na.rm = TRUE))
})
checks <- do.call(rbind, checks)
write_tsv(checks, "baseline_reproduction_checks.tsv")
stopifnot(all(checks$matched == 658L), all(checks$hit_mismatches == 0L),
          all(is.finite(checks$max_abs_p_difference)))
if (any(checks$max_abs_p_difference >= 1e-10)) {
  warning("The prechange code and archived submission p-values differ; see baseline_reproduction_checks.tsv.")
}

# Identify the historical no-EVC convention without modifying the frozen sources.
historical <- load_sources(baseline_dir)
baseline_moderator <- historical$eb_moderate_se
historical$eb_moderate_se <- function(se, df) baseline_moderator(se, min(df))
historical_no_evc <- historical$results(fits$without_EVC)
ix <- match(paste(historical_no_evc$reporter, historical_no_evc$perturbation, sep = "::"),
            paste(published$promoter, published$compound, sep = "::"))
historical_check <- data.frame(
  convention = "Minimum residual df used throughout moderation and testing",
  max_abs_p_difference = max(abs(historical_no_evc$specific_pvalue - published$modeled_pvalue[ix])),
  max_abs_q_difference = max(abs(historical_no_evc$specific_padj_by_reporter - published$modeled_padj_by_reporter[ix]))
)
write_tsv(historical_check, "historical_no_evc_reproduction.tsv")
stopifnot(historical_check$max_abs_p_difference < 1e-10)

frozen_files <- c(
  "paper/dstressr_package_manuscript/manuscript.qmd",
  "paper/dstressr_package_manuscript/manuscript.pdf",
  "analysis/outputs/binsfeld_evc_calibrated/evc_huber_comparison_to_binsfeld_and_default.tsv"
)
utils::untar(archive, files = frozen_files, exdir = verification_dir)
preservation <- data.frame(path = frozen_files,
  archived_md5 = unname(tools::md5sum(file.path(verification_dir, frozen_files))),
  current_md5 = unname(tools::md5sum(analysis_path(frozen_files))))
write_tsv(preservation, "submission_promotion_changes.tsv")
stopifnot(!anyNA(preservation$archived_md5))

if (requireNamespace("ggplot2", quietly = TRUE)) {
  selected <- tables[grepl("submission_baseline|single_WS", names(tables))]
  fields <- c("workflow", "method", "specific_pvalue", "baseline_pvalue", "hit", "change")
  plot_data <- do.call(rbind, lapply(selected, function(x) x[, fields]))
  plot_data$method <- factor(plot_data$method, levels = c("submission_baseline", "single_WS"),
                             labels = c("Prechange code", "Single moderation + effective df"))
  plot_data$workflow <- factor(plot_data$workflow, levels = names(fits),
                               labels = c("Without EVC", "With EVC"))
  plot_theme <- ggplot2::theme_bw(base_size = 11) + ggplot2::theme(
    panel.grid = ggplot2::element_blank(), strip.background = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(face = "bold"), legend.position = "bottom")
  hist_plot <- ggplot2::ggplot(plot_data, ggplot2::aes(specific_pvalue)) +
    ggplot2::geom_histogram(binwidth = 0.05, boundary = 0, fill = "#BFC3C7", colour = "white") +
    ggplot2::geom_histogram(data = plot_data[plot_data$hit, ], binwidth = 0.05,
                            boundary = 0, fill = "#0072B2", colour = "white") +
    ggplot2::facet_grid(workflow ~ method) +
    ggplot2::labs(x = "Raw p-value", y = "Number of reporter-compound pairs") + plot_theme
  scatter_data <- plot_data[plot_data$method == "Single moderation + effective df", ]
  scatter_plot <- ggplot2::ggplot(scatter_data, ggplot2::aes(
    -log10(pmax(baseline_pvalue, 1e-300)), -log10(pmax(specific_pvalue, 1e-300)), colour = change)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    ggplot2::geom_point(size = 1.3, alpha = 0.75) + ggplot2::facet_wrap(~workflow) +
    ggplot2::scale_colour_manual(values = c(Unchanged = "#9B9FA3", Gained = "#0072B2", Lost = "#D55E00")) +
    ggplot2::labs(x = expression(-log[10](p[baseline])), y = expression(-log[10](p[candidate])),
                  colour = NULL) + plot_theme
  ggplot2::ggsave(file.path(out_dir, "pvalue_histogram_comparison.png"), hist_plot,
                  width = 9, height = 6, dpi = 250)
  ggplot2::ggsave(file.path(out_dir, "pvalue_comparison.png"), scatter_plot,
                  width = 9, height = 4.3, dpi = 250)
}
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
print(checks)
print(summary_table, row.names = FALSE)
message("Candidate inference outputs are isolated in: ", out_dir)
