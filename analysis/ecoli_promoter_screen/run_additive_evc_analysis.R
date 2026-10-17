# Run from the repository root. Existing paper fits and figures are read only.
source(file.path("analysis", "_helpers.R"))
load_destress_package()
load_binsfeld_paper_data()
set.seed(1)

out_dir <- analysis_output_dir("binsfeld_additive_evc")
write_tsv <- function(x, name) {
  utils::write.table(x, file.path(out_dir, name), sep = "\t",
                     quote = FALSE, row.names = FALSE)
}
wt_auc <- binsfeld_reporter_auc[
  binsfeld_reporter_auc$strain == "WT" & binsfeld_reporter_auc$removed == "No", ]

assay <- prepare_assay(
  wt_auc, reporter = "promoter", perturbation = "compound", control = "Water",
  lux = "lux_auc", growth = "od_auc", growth_exponent = "estimate",
  batch = "dose_level", replicate = "replicate",
  growth_covariates = "replicate", numeric_covariates = "dose_level",
  background_reporter = "EVC", background_method = "subtract",
  background_by = c("compound", "dose_level", "replicate")
)
is_tested <- as.character(assay$.reporter) != "EVC"
matched <- is_tested & is.finite(assay$.background_response)
stopifnot(isTRUE(all.equal(
  assay$.response[matched],
  assay$.response_uncalibrated[matched] - assay$.background_response[matched]
)))
fit <- fit_destress(
  assay, technical = c("replicate", "dose_level"), empirical_bayes = TRUE,
  adjustment = "by_reporter", interaction = FALSE
)
saveRDS(fit, file.path(out_dir, "additive_evc_fit.rds"))
params <- model_parameters(fit)
write_tsv(params$growth_exponents, "growth_exponents.tsv")
write_tsv(params$background_calibration, "background_calibration.tsv")

fit_paths <- c(
  without_evc = analysis_path("analysis", "outputs", "binsfeld_workflow_comparison",
                              "destress_estimated_exponent_fit.rds"),
  huber_evc = analysis_path("analysis", "outputs", "binsfeld_evc_calibrated",
                            "evc_huber_fit.rds")
)
fits <- lapply(fit_paths, readRDS)
fits$additive_evc <- fit
stopifnot(isTRUE(all.equal(params$growth_exponents,
                          model_parameters(fits$huber_evc)$growth_exponents)))
stopifnot(nrow(fit$assay_data) == nrow(fits$huber_evc$assay_data))
tables <- lapply(fits, results)
pair_key <- function(d) paste(d$reporter, d$perturbation, sep = "\r")
keys <- pair_key(tables$additive_evc)
stopifnot(!anyDuplicated(keys))
tables <- lapply(tables, function(d) {
  stopifnot(!anyDuplicated(pair_key(d)), setequal(keys, pair_key(d)))
  d <- d[match(keys, pair_key(d)), , drop = FALSE]
  d$total_hit <- !is.na(d$total_padj_by_reporter) & d$total_padj_by_reporter < 0.05
  d$relative_hit <- !is.na(d$relative_padj_by_reporter) & d$relative_padj_by_reporter < 0.05
  d
})
write_tsv(tables$additive_evc, "additive_evc_full_results.tsv")
write_tsv(tables$additive_evc[tables$additive_evc$total_hit, ],
          "additive_evc_significant_total_pairs.tsv")
write_tsv(tables$additive_evc[tables$additive_evc$relative_hit, ],
          "additive_evc_significant_specific_pairs.tsv")

summary <- do.call(rbind, lapply(names(tables), function(method) {
  d <- tables[[method]]
  model_counts <- unique(fits[[method]]$reporter_effects[, c("reporter", "n_observations")])
  stopifnot(!anyDuplicated(model_counts$reporter))
  data.frame(method = method, input_observations = nrow(fits[[method]]$assay_data),
             modeled_observations = sum(model_counts$n_observations),
             pairs = nrow(d), total_hits = sum(d$total_hit),
             relative_hits = sum(d$relative_hit))
}))
write_tsv(summary, "hit_counts.tsv")
comparison <- tables$additive_evc[, c("reporter", "perturbation")]
for (method in names(tables)) {
  columns <- c("total_effect", "total_se", "total_pvalue", "total_padj_by_reporter",
               "total_hit", "relative_effect", "relative_padj_by_reporter", "relative_hit")
  for (column in columns) comparison[[paste(method, column, sep = "_")]] <- tables[[method]][[column]]
}
write_tsv(comparison, "three_method_pair_comparison.tsv")
overlaps <- do.call(rbind, lapply(c("without_evc", "huber_evc"), function(method) {
  a <- tables$additive_evc$total_hit
  b <- tables[[method]]$total_hit
  data.frame(comparator = method, shared_total_hits = sum(a & b),
             additive_only = sum(a & !b), comparator_only = sum(!a & b),
             total_effect_correlation = cor(tables$additive_evc$total_effect,
                                            tables[[method]]$total_effect))
}))
write_tsv(overlaps, "total_hit_overlaps.tsv")
write_tsv(data.frame(
  tested_input_rows = sum(is_tested), matched_rows = sum(matched),
  unmatched_rows = sum(is_tested & !is.finite(assay$.background_response)),
  subtraction_max_error = max(abs(assay$.response[matched] -
    (assay$.response_uncalibrated[matched] - assay$.background_response[matched])))
), "matching_checks.tsv")
saveRDS(list(results = tables, summary = summary, overlaps = overlaps,
             parameters = params, comparison = comparison,
             provenance = list(comparator_fits = fit_paths, input = "binsfeld_reporter_auc",
                               background_by = c("compound", "dose_level", "replicate"))),
        file.path(out_dir, "additive_evc_comparison.rds"))

if (!requireNamespace("ggplot2", quietly = TRUE) ||
    !requireNamespace("gridExtra", quietly = TRUE)) {
  stop("Install ggplot2 and gridExtra to draw comparison figures.")
}
labels <- c(without_evc = "Without EVC", huber_evc = "Huber EVC calibration",
            additive_evc = "Direct EVC subtraction")
plot_data <- do.call(rbind, lapply(names(tables), function(method) {
  d <- tables[[method]]
  d$method <- factor(labels[[method]], levels = unname(labels))
  d
}))
histogram <- ggplot2::ggplot(plot_data, ggplot2::aes(x = total_pvalue)) +
  ggplot2::geom_histogram(breaks = seq(0, 1, by = 0.05), fill = "#CBD5DC",
                          color = "white", linewidth = 0.25) +
  ggplot2::geom_histogram(data = plot_data[plot_data$total_hit, ],
                          breaks = seq(0, 1, by = 0.05), fill = "#0072B2",
                          color = "white", linewidth = 0.25) +
  ggplot2::facet_wrap(~method, nrow = 1) +
  ggplot2::labs(x = "P-value", y = "Number of pairs") +
  ggplot2::theme_light(base_size = 12) +
  ggplot2::theme(strip.background = ggplot2::element_rect(fill = "#354650"),
                 strip.text = ggplot2::element_text(color = "white"))
ggplot2::ggsave(file.path(out_dir, "total_pvalue_histograms.png"), histogram,
                width = 13.5, height = 3.4, dpi = 220, bg = "white")
xlim <- range(plot_data$total_effect, finite = TRUE) * 1.12
ylim <- c(0, max(-log10(pmax(plot_data$total_pvalue, .Machine$double.xmin))) * 1.12)
volcanos <- lapply(names(tables), function(method) {
  plot_volcano(tables[[method]], effect = "total_effect",
    padj = "total_padj_by_reporter", pvalue = "total_pvalue", top_n = 5,
    color_by = "significance", title = NULL, xlab = expression(hat(Delta)[aj]^tot),
    call_labels = c("Not significant", "Significant")) +
    ggplot2::coord_cartesian(xlim = xlim, ylim = ylim, clip = "off") +
    ggplot2::labs(subtitle = labels[[method]], color = NULL) +
    ggplot2::theme(plot.subtitle = ggplot2::element_text(size = 12, face = "bold"))
})
composite <- gridExtra::arrangeGrob(
  histogram, gridExtra::arrangeGrob(grobs = volcanos, ncol = 3),
  ncol = 1, heights = c(0.9, 1.35)
)
ggplot2::ggsave(file.path(out_dir, "total_effect_comparison.png"), composite,
                width = 13.5, height = 8.5, dpi = 220, bg = "white")
capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))
print(summary, row.names = FALSE)
print(overlaps, row.names = FALSE)
message("Saved additive EVC analysis to ", out_dir)
