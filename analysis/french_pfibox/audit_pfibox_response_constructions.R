#!/usr/bin/env Rscript

# Compare PFIbox response constructions in DStressR notation.
#
# The French et al. workflow naturally exposes three response objects:
# 1. endpoint fluorescence fold-change against the untreated/reference plates,
# 2. biomass-normalized fluorescence q_i(t),
# 3. promoter activity PA_i(t), approximately d q_i(t) / dt.
#
# This script builds compact colony-level summaries for a subset by default and
# fits response-only DStressR models to the biomass-normalized and promoter-
# activity summaries. It is an audit script, not a manuscript figure script.

source(file.path("analysis", "_helpers.R"))
load_destress_package()

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package `ggplot2` is required for PFIbox response-construction audit.", call. = FALSE)
}
if (!requireNamespace("gridExtra", quietly = TRUE)) {
  stop("Package `gridExtra` is required for PFIbox response-construction audit.", call. = FALSE)
}
if (!requireNamespace("ggrepel", quietly = TRUE)) {
  stop("Package `ggrepel` is required for PFIbox response-construction audit.", call. = FALSE)
}

ggplot2 <- asNamespace("ggplot2")
gridExtra <- asNamespace("gridExtra")

source_file <- analysis_path("data-raw", "french_pfibox", "source", "gfpData.Rdata")
out_dir <- analysis_output_dir("french_pfibox", "response_construction_audit")
if (!file.exists(source_file)) {
  stop("Missing PFIbox object: ", source_file, call. = FALSE)
}

reporter_limit <- suppressWarnings(as.integer(Sys.getenv("PFIBOX_REPORTER_LIMIT", "80")))
perturbation_limit <- suppressWarnings(as.integer(Sys.getenv("PFIBOX_PERTURBATION_LIMIT", "8")))
late_start_h <- suppressWarnings(as.numeric(Sys.getenv("PFIBOX_LATE_START_H", "2")))
run_tag <- paste0("reporters", reporter_limit, "_perturbations", perturbation_limit, "_late", gsub("[.]", "p", late_start_h))

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

load(source_file)
if (!exists("gfpData") || !is.list(gfpData) || is.null(gfpData$legend)) {
  stop("Expected a list object `gfpData` with a common `legend`.", call. = FALSE)
}

legend <- as.data.frame(gfpData$legend, stringsAsFactors = FALSE)
treatment_names <- setdiff(names(gfpData), "legend")
treatment_metadata <- do.call(rbind, lapply(treatment_names, parse_treatment))
treatment_metadata <- treatment_metadata[order(treatment_metadata$is_control, treatment_metadata$treatment), ]
perturbations <- sort(setdiff(treatment_metadata$perturbation, "Control"))
keep_perturbations <- perturbations[seq_len(min(perturbation_limit, length(perturbations)))]
keep_treatments <- treatment_metadata$treatment[
  treatment_metadata$perturbation %in% c("Control", keep_perturbations)
]

reporters <- sort(setdiff(unique(as.character(legend$Promoter)), c("Empty", "U66", "U139")))
keep_reporters <- reporters[seq_len(min(reporter_limit, length(reporters)))]

summarize_treatment <- function(treatment_name) {
  tr <- gfpData[[treatment_name]]
  if (is.null(tr$raw.loess) || is.null(tr$constit.normalized) || is.null(tr$promoter.activity)) {
    return(NULL)
  }
  meta <- treatment_metadata[treatment_metadata$treatment == treatment_name, , drop = FALSE]
  raw <- as.matrix(tr$raw.loess)
  q <- as.matrix(tr$constit.normalized)
  pa <- as.matrix(tr$promoter.activity)
  keep <- as.character(legend$Promoter) %in% keep_reporters
  pos <- which(keep)
  time_raw_h <- (seq_len(ncol(raw)) - 1) * 5 / 60
  time_q_h <- (seq_len(ncol(q)) - 1) * 5 / 60
  # promoter.activity has fewer columns because the original workflow excludes
  # early points before taking/using the derivative. We align its first column
  # to the end of the dropped prefix for diagnostics.
  pa_offset <- max(0, ncol(q) - ncol(pa))
  time_pa_h <- (seq_len(ncol(pa)) + pa_offset - 1) * 5 / 60
  late_q <- time_q_h >= late_start_h
  late_pa <- time_pa_h >= late_start_h

  data.frame(
    reporter = as.character(legend$Promoter[pos]),
    perturbation = meta$perturbation,
    treatment = treatment_name,
    drug = meta$drug,
    concentration_code = as.integer(meta$concentration_code),
    is_control = meta$is_control,
    colony_position = pos,
    raw_endpoint = rowMeans(raw[pos, seq(max(1, ncol(raw) - 11), ncol(raw)), drop = FALSE], na.rm = TRUE),
    q_auc_all = rowSums(q[pos, , drop = FALSE], na.rm = TRUE),
    q_auc_late = rowSums(q[pos, late_q, drop = FALSE], na.rm = TRUE),
    q_endpoint = rowMeans(q[pos, seq(max(1, ncol(q) - 11), ncol(q)), drop = FALSE], na.rm = TRUE),
    pa_auc_all = rowSums(pa[pos, , drop = FALSE], na.rm = TRUE),
    pa_auc_late = rowSums(pa[pos, late_pa, drop = FALSE], na.rm = TRUE),
    pa_mean_late = rowMeans(pa[pos, late_pa, drop = FALSE], na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

response_parts <- lapply(keep_treatments, summarize_treatment)
response_parts <- Filter(Negate(is.null), response_parts)
response_tab <- do.call(rbind, response_parts)
response_tab <- response_tab[is.finite(response_tab$q_auc_late), , drop = FALSE]

control_endpoint <- stats::aggregate(
  raw_endpoint ~ reporter,
  data = response_tab[response_tab$perturbation == "Control", , drop = FALSE],
  FUN = mean,
  na.rm = TRUE
)
names(control_endpoint)[2] <- "control_raw_endpoint"
response_tab <- merge(response_tab, control_endpoint, by = "reporter", all.x = TRUE, sort = FALSE)
response_tab$log2_raw_endpoint_fold_change <- log2(
  (response_tab$raw_endpoint + 1e-8) / (response_tab$control_raw_endpoint + 1e-8)
)
response_tab$log2_q_auc_late <- log2(response_tab$q_auc_late + 1e-8)
response_tab$log2_q_endpoint <- log2(response_tab$q_endpoint + 1e-8)

response_file <- file.path(out_dir, paste0("pfibox_response_construction_table_", run_tag, ".tsv"))
utils::write.table(response_tab, response_file, sep = "\t", quote = FALSE, row.names = FALSE)

fit_response <- function(value, label) {
  d <- response_tab[is.finite(response_tab[[value]]), , drop = FALSE]
  assay <- prepare_assay(
    d,
    reporter = "reporter",
    perturbation = "perturbation",
    control = "Control",
    response = value
  )
  fit <- fit_destress(
    assay,
    empirical_bayes = TRUE,
    adjustment = "by_reporter",
    interaction = FALSE
  )
  res <- results(fit)
  res$response_construction <- label
  res
}

results_list <- list(
  endpoint_fold_change = fit_response("log2_raw_endpoint_fold_change", "endpoint fluorescence fold-change"),
  q_auc_late = fit_response("log2_q_auc_late", "log2 late AUC of biomass-normalized fluorescence"),
  pa_mean_late = fit_response("pa_mean_late", "late mean promoter activity")
)
result_tab <- do.call(rbind, results_list)
result_file <- file.path(out_dir, paste0("pfibox_response_construction_results_", run_tag, ".tsv"))
utils::write.table(result_tab, result_file, sep = "\t", quote = FALSE, row.names = FALSE)

plot_data <- rbind(
  data.frame(
    response = "Endpoint fluorescence fold-change",
    value = response_tab$log2_raw_endpoint_fold_change,
    stringsAsFactors = FALSE
  ),
  data.frame(
    response = "Late AUC of q(t)",
    value = response_tab$log2_q_auc_late,
    stringsAsFactors = FALSE
  ),
  data.frame(
    response = "Late mean PA(t)",
    value = response_tab$pa_mean_late,
    stringsAsFactors = FALSE
  )
)
plot_data <- plot_data[is.finite(plot_data$value), , drop = FALSE]
plot_data$response <- factor(plot_data$response, levels = c(
  "Endpoint fluorescence fold-change",
  "Late AUC of q(t)",
  "Late mean PA(t)"
))

p_distributions <- ggplot2$ggplot(plot_data, ggplot2$aes(x = value)) +
  ggplot2$geom_histogram(bins = 60, fill = "#4E79A7", color = "white", linewidth = 0.15) +
  ggplot2$facet_wrap(ggplot2$vars(response), scales = "free", nrow = 1) +
  ggplot2$labs(x = NULL, y = "Number of colony observations") +
  ggplot2$theme_light(base_size = 10) +
  ggplot2$theme(
    panel.grid.minor = ggplot2$element_blank(),
    panel.grid.major = ggplot2$element_line(linewidth = 0.25, color = "#E5E7EB"),
    strip.background = ggplot2$element_rect(fill = "#111827", color = NA),
    strip.text = ggplot2$element_text(color = "white", face = "bold", size = 8.5),
    axis.title = ggplot2$element_text(color = "#111827", size = 9),
    axis.text = ggplot2$element_text(color = "#4B5563", size = 8)
  )

pvalues <- result_tab[
  is.finite(result_tab$specific_pvalue),
  c("response_construction", "specific_pvalue"),
  drop = FALSE
]
pvalues$response_construction <- factor(
  pvalues$response_construction,
  levels = c(
    "endpoint fluorescence fold-change",
    "log2 late AUC of biomass-normalized fluorescence",
    "late mean promoter activity"
  ),
  labels = c(
    "Endpoint fluorescence fold-change",
    "Late AUC of q(t)",
    "Late mean PA(t)"
  )
)
p_pvalues <- ggplot2$ggplot(pvalues, ggplot2$aes(x = specific_pvalue)) +
  ggplot2$geom_histogram(bins = 50, fill = "#2A9D8F", color = "white", linewidth = 0.15) +
  ggplot2$facet_wrap(ggplot2$vars(response_construction), nrow = 1) +
  ggplot2$labs(x = "Specific-effect p-value", y = "Number of reporter-perturbation pairs") +
  ggplot2$theme_light(base_size = 10) +
  ggplot2$theme(
    panel.grid.minor = ggplot2$element_blank(),
    panel.grid.major = ggplot2$element_line(linewidth = 0.25, color = "#E5E7EB"),
    strip.background = ggplot2$element_rect(fill = "#111827", color = NA),
    strip.text = ggplot2$element_text(color = "white", face = "bold", size = 8.5),
    axis.title = ggplot2$element_text(color = "#111827", size = 9),
    axis.text = ggplot2$element_text(color = "#4B5563", size = 8)
  )

combined <- gridExtra$arrangeGrob(
  p_distributions,
  p_pvalues,
  ncol = 1,
  heights = c(1, 1)
)

ggplot2$ggsave(
  file.path(out_dir, paste0("pfibox_response_construction_audit_", run_tag, ".png")),
  combined,
  width = 12.2,
  height = 6.6,
  dpi = 280,
  bg = "white"
)
ggplot2$ggsave(
  file.path(out_dir, paste0("pfibox_response_construction_audit_", run_tag, ".pdf")),
  combined,
  width = 12.2,
  height = 6.6,
  bg = "white"
)

message("Wrote response-construction table: ", response_file)
message("Wrote response-construction DStressR results: ", result_file)
