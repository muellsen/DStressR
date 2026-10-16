#!/usr/bin/env Rscript

# Reproducible subset DStressR workflow for the French et al. PFIbox data.
#
# The primary response is the Zaslaver/French promoter-activity summary
#
#   y_i^{PA}(T) = |T|^{-1} int_T PA_i(t) dt,
#
# approximated by the mean of the deposited promoter.activity time series over a
# user-defined time window T. Reporter and perturbation limits are configurable:
# the defaults are 250 reporters and 8 non-control perturbations, i.e. about 2,000
# reporter-perturbation pairs before missing-value filtering.

source(file.path("analysis", "_helpers.R"))
load_destress_package()

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package `ggplot2` is required for the PFIbox PA subset workflow.", call. = FALSE)
}
if (!requireNamespace("gridExtra", quietly = TRUE)) {
  stop("Package `gridExtra` is required for the PFIbox PA subset workflow.", call. = FALSE)
}

source_file <- analysis_path("data-raw", "french_pfibox", "source", "gfpData.Rdata")
out_dir <- analysis_output_dir("french_pfibox", "promoter_activity_analysis")
processed_dir <- analysis_path("data-raw", "french_pfibox", "processed")
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

reporter_limit <- suppressWarnings(as.integer(Sys.getenv("PFIBOX_REPORTER_LIMIT", "250")))
perturbation_limit <- suppressWarnings(as.integer(Sys.getenv("PFIBOX_PERTURBATION_LIMIT", "8")))
t_start_h <- suppressWarnings(as.numeric(Sys.getenv("PFIBOX_T_START_H", "2")))
t_end_h <- suppressWarnings(as.numeric(Sys.getenv("PFIBOX_T_END_H", "18")))
design <- Sys.getenv("PFIBOX_DESIGN", "dose_replicate")
design <- match.arg(
  design,
  c("none", "replicate", "dose", "dose_replicate", "dose_replicate_drug")
)
perturbation_target <- Sys.getenv("PFIBOX_PERTURBATION_TARGET", "treatment")
perturbation_target <- match.arg(perturbation_target, c("treatment", "drug"))
model_path <- Sys.getenv("PFIBOX_MODEL_PATH", "scalable")
model_path <- match.arg(model_path, c("scalable", "interaction"))
empirical_bayes <- tolower(Sys.getenv("PFIBOX_EMPIRICAL_BAYES", "true"))
empirical_bayes <- empirical_bayes %in% c("1", "true", "yes", "y")

if (!is.finite(reporter_limit) || reporter_limit < 1) {
  stop("PFIBOX_REPORTER_LIMIT must be a positive integer.", call. = FALSE)
}
if (!is.finite(perturbation_limit) || perturbation_limit < 1) {
  stop("PFIBOX_PERTURBATION_LIMIT must be a positive integer.", call. = FALSE)
}
if (!is.finite(t_start_h) || !is.finite(t_end_h) || t_start_h >= t_end_h) {
  stop("PFIBOX_T_START_H and PFIBOX_T_END_H must define an increasing time window.", call. = FALSE)
}
response_cache <- Sys.getenv("PFIBOX_RESPONSE_TABLE", "")
if (!nzchar(response_cache) && !file.exists(source_file)) {
  stop("Missing PFIbox object: ", source_file, call. = FALSE)
}

run_tag <- paste0(
  "reporters", reporter_limit,
  "_perturbations", perturbation_limit,
  "_T", gsub("[.]", "p", t_start_h), "-", gsub("[.]", "p", t_end_h), "h",
  "_", perturbation_target,
  "_", design,
  "_", model_path,
  if (isTRUE(empirical_bayes)) "" else "_rawse"
)

parse_treatment <- function(x) {
  m <- regexec("^(.*)-([0-9]+)$", x)
  r <- regmatches(x, m)[[1]]
  if (length(r) == 3) {
    data.frame(
      treatment = x,
      drug = r[2],
      concentration_code = r[3],
      perturbation = if (identical(r[3], "0")) "Control" else x,
      is_control = identical(r[3], "0"),
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      treatment = x,
      drug = x,
      concentration_code = NA_character_,
      perturbation = x,
      is_control = FALSE,
      stringsAsFactors = FALSE
    )
  }
}

if (nzchar(response_cache)) {
  message("Reading previously summarized PA responses: ", response_cache)
  response_tab <- utils::read.delim(response_cache, check.names = FALSE)
  required <- c("reporter", "perturbation", "drug_perturbation", "y_pa",
                "dose_code_centered", "drug_factor", "replicate",
                "pa_window_start_h", "pa_window_end_h", "pa_n_timepoints")
  stopifnot(all(required %in% names(response_tab)))
  if (length(unique(response_tab$reporter)) > reporter_limit ||
      length(setdiff(unique(response_tab$perturbation), "Control")) > perturbation_limit ||
      any(abs(response_tab$pa_window_start_h - t_start_h) > 1e-8) ||
      any(abs(response_tab$pa_window_end_h - t_end_h) > 1e-8)) {
    stop("Cached response design/window does not match the requested run.", call. = FALSE)
  }
  response_tab$drug_factor <- factor(response_tab$drug_factor)
  response_tab$replicate <- factor(response_tab$replicate)
} else {
message("Loading PFIbox object: ", source_file)
load(source_file)
if (!exists("gfpData") || !is.list(gfpData) || is.null(gfpData$legend)) {
  stop("Expected a list object `gfpData` with a common `legend`.", call. = FALSE)
}

legend <- as.data.frame(gfpData$legend, stringsAsFactors = FALSE)
treatment_names <- setdiff(names(gfpData), "legend")
treatment_metadata <- do.call(rbind, lapply(treatment_names, parse_treatment))

perturbations <- sort(setdiff(treatment_metadata$perturbation, "Control"))
keep_perturbations <- perturbations[seq_len(min(perturbation_limit, length(perturbations)))]
keep_treatments <- treatment_metadata$treatment[
  treatment_metadata$perturbation %in% c("Control", keep_perturbations)
]

reporters <- sort(setdiff(unique(as.character(legend$Promoter)), c("Empty", "U66", "U139")))
keep_reporters <- reporters[seq_len(min(reporter_limit, length(reporters)))]

summarize_treatment_pa <- function(treatment_name) {
  tr <- gfpData[[treatment_name]]
  if (is.null(tr$promoter.activity)) {
    return(NULL)
  }
  pa <- as.matrix(tr$promoter.activity)
  meta <- treatment_metadata[treatment_metadata$treatment == treatment_name, , drop = FALSE]
  keep <- as.character(legend$Promoter) %in% keep_reporters
  pos <- which(keep)

  # The deposited promoter.activity matrix has fewer columns than raw.loess
  # because the published workflow drops an early prefix before derivative-based
  # activity analysis. We align its columns to the common five-minute acquisition
  # grid by accounting for that dropped prefix when raw.loess is available.
  n_raw <- if (!is.null(tr$raw.loess)) ncol(as.matrix(tr$raw.loess)) else ncol(pa)
  pa_offset <- max(0, n_raw - ncol(pa))
  time_pa_h <- (seq_len(ncol(pa)) + pa_offset - 1) * 5 / 60
  in_window <- time_pa_h >= t_start_h & time_pa_h <= t_end_h
  if (!any(in_window)) {
    stop(
      "No promoter.activity time points in requested window [",
      t_start_h, ", ", t_end_h, "] hours.",
      call. = FALSE
    )
  }

  data.frame(
    reporter = as.character(legend$Promoter[pos]),
    perturbation = meta$perturbation,
    treatment = treatment_name,
    drug = meta$drug,
    concentration_code = suppressWarnings(as.integer(meta$concentration_code)),
    is_control = meta$is_control,
    colony_position = pos,
    replicate = ave(
      seq_along(legend$Promoter[pos]),
      legend$Promoter[pos],
      FUN = seq_along
    ),
    y_pa = rowMeans(pa[pos, in_window, drop = FALSE], na.rm = TRUE),
    pa_window_start_h = min(time_pa_h[in_window]),
    pa_window_end_h = max(time_pa_h[in_window]),
    pa_n_timepoints = sum(in_window),
    stringsAsFactors = FALSE
  )
}

message(
  "Building PA response table for ", length(keep_reporters), " reporters and ",
  length(keep_perturbations), " non-control perturbations."
)
response_parts <- lapply(keep_treatments, summarize_treatment_pa)
response_parts <- Filter(Negate(is.null), response_parts)
response_tab <- do.call(rbind, response_parts)
response_tab <- response_tab[is.finite(response_tab$y_pa), , drop = FALSE]
response_tab$dose_code <- as.numeric(response_tab$concentration_code)
response_tab$dose_code[!is.finite(response_tab$dose_code)] <- 0
response_tab$dose_code_centered <- response_tab$dose_code - mean(response_tab$dose_code, na.rm = TRUE)
response_tab$drug_factor <- factor(response_tab$drug)
response_tab$drug_perturbation <- ifelse(response_tab$is_control, "Control", response_tab$drug)
response_tab$replicate <- factor(response_tab$replicate)
}
perturbation_column <- if (identical(perturbation_target, "drug")) "drug_perturbation" else "perturbation"

response_file <- file.path(out_dir, paste0("pfibox_pa_response_table_", run_tag, ".tsv"))
utils::write.table(response_tab, response_file, sep = "\t", quote = FALSE, row.names = FALSE)

expected_pairs <- length(unique(response_tab$reporter)) *
  length(setdiff(unique(response_tab[[perturbation_column]]), "Control"))
message("Subset contains about ", expected_pairs, " tested reporter-perturbation pairs.")
message("DStressR perturbation target: ", perturbation_column)

technical <- switch(
  design,
  none = character(),
  replicate = "replicate",
  dose = "dose_code_centered",
  dose_replicate = c("dose_code_centered", "replicate"),
  dose_replicate_drug = c("dose_code_centered", "replicate", "drug_factor")
)
message("DStressR technical design: ", if (length(technical)) paste(technical, collapse = " + ") else "none")
message("DStressR model path: ", model_path)
message("DStressR empirical-Bayes moderation: ", empirical_bayes)

assay <- prepare_assay(
  response_tab,
  reporter = "reporter",
  perturbation = perturbation_column,
  control = "Control",
  response = "y_pa",
  batch = "dose_code_centered",
  plate = "drug_factor",
  replicate = "replicate",
  numeric_covariates = "dose_code_centered"
)

fit <- fit_destress(
  assay,
  technical = technical,
  empirical_bayes = empirical_bayes,
  adjustment = "by_reporter",
  interaction = identical(model_path, "interaction"),
  verbose = TRUE
)
saveRDS(fit, file.path(out_dir, paste0("pfibox_pa_fit_", run_tag, ".rds")))
res <- results(fit, verbose = TRUE)
res$response <- "pa_mean_window"
res$t_start_h <- t_start_h
res$t_end_h <- t_end_h

result_file <- file.path(out_dir, paste0("pfibox_pa_destress_pair_results_", run_tag, ".tsv"))
utils::write.table(res, result_file, sep = "\t", quote = FALSE, row.names = FALSE)

diagnostics <- perturbation_diagnostics(
  res,
  mean_effect = "total_effect",
  variance_effect = "total_effect",
  reporter = "reporter",
  perturbation = "perturbation",
  min_reporters = 10
)
diagnostic_file <- file.path(out_dir, paste0("pfibox_pa_perturbation_diagnostics_", run_tag, ".tsv"))
utils::write.table(diagnostics, diagnostic_file, sep = "\t", quote = FALSE, row.names = FALSE)

metadata <- data.frame(
  run_tag = run_tag,
  reporter_limit = reporter_limit,
  perturbation_limit = perturbation_limit,
  n_reporters = length(unique(response_tab$reporter)),
  n_noncontrol_perturbations = length(setdiff(unique(response_tab[[perturbation_column]]), "Control")),
  n_response_rows = nrow(response_tab),
  n_result_pairs = nrow(res),
  requested_t_start_h = t_start_h,
  requested_t_end_h = t_end_h,
  realized_t_start_h = min(response_tab$pa_window_start_h, na.rm = TRUE),
  realized_t_end_h = max(response_tab$pa_window_end_h, na.rm = TRUE),
  pa_n_timepoints = unique(response_tab$pa_n_timepoints)[1],
  design = design,
  perturbation_target = perturbation_target,
  perturbation_column = perturbation_column,
  model_path = model_path,
  empirical_bayes = empirical_bayes,
  technical = if (length(technical)) paste(technical, collapse = " + ") else "none",
  response_source = if (nzchar(response_cache)) normalizePath(response_cache) else source_file,
  stringsAsFactors = FALSE
)
metadata_file <- file.path(out_dir, paste0("pfibox_pa_run_metadata_", run_tag, ".tsv"))
utils::write.table(metadata, metadata_file, sep = "\t", quote = FALSE, row.names = FALSE)

theme_paper <- function(base_size = 10) {
  ggplot2::theme_light(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(linewidth = 0.25, color = "#E5E7EB"),
      axis.title = ggplot2::element_text(color = "#111827"),
      axis.text = ggplot2::element_text(color = "#4B5563"),
      legend.position = "right",
      legend.justification = "center",
      legend.background = ggplot2::element_blank(),
      legend.key = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(size = base_size - 1),
      legend.text = ggplot2::element_text(size = base_size - 2),
      plot.title = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_blank()
    )
}

p_hist <- ggplot2::ggplot(res[is.finite(res$specific_pvalue), , drop = FALSE],
                          ggplot2::aes(specific_pvalue)) +
  ggplot2::geom_histogram(bins = 50, fill = "#4E79A7", color = "white", linewidth = 0.15) +
  ggplot2::labs(
    x = "Specific-effect p-value",
    y = "Number of reporter-perturbation pairs"
  ) +
  theme_paper()

p_volcano <- plot_volcano(
  res,
  effect = "specific_effect",
  pvalue = "specific_pvalue",
  padj = "specific_padj_by_reporter",
  y = "pvalue",
  reporter = "reporter",
  perturbation = "perturbation",
  fdr = 0.05,
  top_n = 14,
  title = NULL,
  subtitle = NULL,
  xlab = expression(hat(delta)[aj]^PA),
  ylab = expression(-log[10]("p-value")),
  call_labels = c("Not called", "FDR < 0.05"),
  point_size = 1.1,
  point_alpha = 0.72,
  hit_color = "#D95F02",
  background_color = "#CBD5E1"
) +
  ggplot2::theme(
    plot.title = ggplot2::element_blank(),
    plot.subtitle = ggplot2::element_blank(),
    legend.position = "right",
    legend.justification = "center",
    legend.background = ggplot2::element_blank(),
    legend.key = ggplot2::element_blank()
  ) +
  ggplot2::guides(
    color = ggplot2::guide_legend(
      title = NULL,
      override.aes = list(size = 2.2, alpha = 1)
    )
  )

p_rank <- plot_mean_variance_diagnostic(
  diagnostics = diagnostics,
  add_trend = TRUE,
  label_by = "variance",
  top_n = 8,
  title = NULL,
  subtitle = NULL,
  xlab = expression("Rank of " * abs(hat(Delta * bar(y))[j]^"tot, PA")),
  ylab = expression("Variance of " * hat(Delta * y)[aj]^"tot, PA")
) +
  ggplot2::theme(
    plot.title = ggplot2::element_blank(),
    plot.subtitle = ggplot2::element_blank(),
    legend.position = "right",
    legend.background = ggplot2::element_blank(),
    legend.key = ggplot2::element_blank()
  )

p_heatmap <- plot_response_heatmap(
  res,
  value = "total_effect",
  reporter = "reporter",
  perturbation = "perturbation",
  show_perturbation_ids = FALSE,
  top_n_perturbations = Inf,
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  show_perturbation_labels = TRUE,
  perturbation_label_angle = 45,
  title = NULL,
  subtitle = NULL,
  legend_title = expression(hat(Delta*y)[aj]^"tot, PA"),
  xlab = NULL,
  ylab = NULL
) +
  ggplot2::theme(
    plot.title = ggplot2::element_blank(),
    plot.subtitle = ggplot2::element_blank(),
    legend.position = "right",
    legend.background = ggplot2::element_blank(),
    legend.key = ggplot2::element_blank()
  ) +
  ggplot2::guides(
    fill = ggplot2::guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barheight = grid::unit(30, "mm"),
      barwidth = grid::unit(4, "mm")
    )
  )

summary_figure <- gridExtra::arrangeGrob(
  p_heatmap,
  gridExtra::arrangeGrob(p_rank, p_hist, p_volcano, ncol = 3),
  ncol = 1,
  heights = c(1.35, 1)
)

figure_base <- file.path(out_dir, paste0("pfibox_pa_destress_subset_summary_", run_tag))
ggplot2::ggsave(paste0(figure_base, ".png"), summary_figure, width = 13.2, height = 9.0, dpi = 280, bg = "white")
ggplot2::ggsave(paste0(figure_base, ".pdf"), summary_figure, width = 13.2, height = 9.0, bg = "white")

message("Wrote PA response table: ", response_file)
message("Wrote DStressR pair results: ", result_file)
message("Wrote perturbation diagnostics: ", diagnostic_file)
message("Wrote run metadata: ", metadata_file)
message("Wrote summary figure: ", paste0(figure_base, ".png"))
