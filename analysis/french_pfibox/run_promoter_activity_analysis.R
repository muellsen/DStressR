#!/usr/bin/env Rscript
# Retained PFIbox workflow: source fluorescence -> unscaled PA integral -> fit -> tests.
# Run from the repository root. See README.md for source provenance and commands.
source(file.path("analysis", "_helpers.R"))
load_destress_package()

source_file <- analysis_path("data-raw", "french_pfibox", "source", "gfpData.Rdata")
out <- Sys.getenv("PFIBOX_OUTPUT_DIR", analysis_output_dir("french_pfibox", "temporal_unscaled_integral_analysis"))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
summary_type <- Sys.getenv("PFIBOX_RESPONSE_SUMMARY", "unscaled_integral")
if (summary_type != "unscaled_integral") stop("This workflow uses unscaled_integral only.")
window <- match.arg(Sys.getenv("PFIBOX_TIME_WINDOW", "all_available"), c("all_available", "requested"))
reporter_limit <- as.integer(Sys.getenv("PFIBOX_REPORTER_LIMIT", "2000"))
perturbation_limit <- as.integer(Sys.getenv("PFIBOX_PERTURBATION_LIMIT", "100"))
t_start <- as.numeric(Sys.getenv("PFIBOX_T_START_H", "1.6666666666666667"))
t_end <- as.numeric(Sys.getenv("PFIBOX_T_END_H", "18"))
stopifnot(is.finite(reporter_limit), reporter_limit > 0,
  is.finite(perturbation_limit), perturbation_limit > 0,
  is.finite(t_start), is.finite(t_end), t_start < t_end)
if (Sys.getenv("PFIBOX_DESIGN", "replicate") != "replicate" ||
    Sys.getenv("PFIBOX_PERTURBATION_TARGET", "treatment") != "treatment" ||
    Sys.getenv("PFIBOX_MODEL_PATH", "scalable") != "scalable") {
  stop("The retained model uses treatment levels, repeated colony position, and scalable reporter-wise fits.")
}
window_tag <- if (window == "all_available") "all_available" else
  paste0(gsub("[.]", "p", t_start), "-", gsub("[.]", "p", t_end), "h")
tag <- paste0("reporters", reporter_limit, "_perturbations", perturbation_limit,
  "_T", window_tag, "_treatment_replicate_scalable_temporal_unscaled_integral")

# 1. Read one source object and preserve its colony and treatment identifiers.
if (!file.exists(source_file)) stop("Missing source; see data-raw/french_pfibox/README.md.")
message("Loading deposited PFIbox source: ", source_file)
e <- new.env(parent = emptyenv())
load(source_file, envir = e)
legend <- as.data.frame(e$gfpData$legend, stringsAsFactors = FALSE)
write.table(cbind(colony_position = seq_len(nrow(legend)), legend),
  file.path(out, "pfibox_colony_layout.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
conditions <- setdiff(names(e$gfpData), "legend")
stopifnot(all(grepl("^.+-[0-9]+$", conditions)))
metadata <- data.frame(treatment = conditions,
  drug = sub("-[0-9]+$", "", conditions),
  concentration_code = as.integer(sub("^.*-", "", conditions)))
metadata$is_control <- metadata$concentration_code == 0L
metadata$perturbation <- ifelse(metadata$is_control, "Control", metadata$treatment)
perturbations <- head(sort(unique(metadata$perturbation[!metadata$is_control])), perturbation_limit)
metadata <- metadata[metadata$is_control | metadata$perturbation %in% perturbations, ]
reporters <- head(sort(setdiff(unique(as.character(legend$Promoter)),
  c("Empty", "U66", "U139"))), reporter_limit)
positions <- which(as.character(legend$Promoter) %in% reporters)
colony_index <- ave(seq_along(positions), legend$Promoter[positions], FUN = seq_along)

# 2. Reconstruct the audited unscaled PA and integrate it on the chosen grid.
# q has 217 samples. Differences 20:216 correspond to right endpoints at
# 100 minutes through 18 hours. Trimming precedes quadratic LOESS (span 0.5).
# The per-hour derivative divides by dt; trapezoidal integration multiplies
# by dt, so the two factors cancel in the final scalar response.
time_index <- seq_len(197L)
time_h <- (20:216) / 12
selected <- if (window == "all_available") rep(TRUE, 197L) else time_h >= t_start & time_h <= t_end
stopifnot(sum(selected) >= 2L)
S <- vapply(time_index, function(k) {
  y <- as.numeric(time_index == k)
  as.numeric(predict(stats::loess(y ~ time_index, span = .5, degree = 2)))
}, numeric(197L))
weights <- numeric(197L)
weights[selected] <- c(.5, rep(1, sum(selected) - 2L), .5)
integration_weights <- as.numeric(crossprod(S, weights))
parts <- lapply(seq_len(nrow(metadata)), function(k) {
  m <- metadata[k, ]
  message("Integrating colony responses: ", m$treatment)
  q <- as.matrix(e$gfpData[[m$treatment]]$constit.normalized)[positions, , drop = FALSE]
  stopifnot(ncol(q) == 217L)
  dq <- q[, 21:217, drop = FALSE] - q[, 20:216, drop = FALSE]
  value <- as.numeric(dq %*% integration_weights)
  value[rowSums(is.finite(dq)) != 197L] <- NA_real_
  data.frame(reporter = as.character(legend$Promoter[positions]),
    perturbation = m$perturbation, treatment = m$treatment, drug = m$drug,
    concentration_code = m$concentration_code, is_control = m$is_control,
    colony_position = positions, replicate = colony_index, y_pa = value,
    response_summary = "unscaled_integral", time_window = window,
    pa_time_alignment = "Audited difference indices 20:216; five-minute right-endpoint grid",
    pa_n_finite = rowSums(is.finite(dq[, selected, drop = FALSE])),
    pa_window_start_h = min(time_h[selected]), pa_window_end_h = max(time_h[selected]),
    pa_n_timepoints = sum(selected), stringsAsFactors = FALSE)
})
rm(e)
invisible(gc())
response <- do.call(rbind, parts)
response <- response[is.finite(response$y_pa), ]
# Keep the same metadata columns as the saved manuscript input. Dose codes are
# source identifiers, not assumed physical concentrations or model covariates.
response$dose_code <- as.numeric(response$concentration_code)
response$dose_code_centered <- response$dose_code - mean(response$dose_code)
response$drug_factor <- factor(response$drug)
response$drug_perturbation <- ifelse(response$is_control, "Control", response$drug)
response$replicate <- factor(response$replicate)
write.table(response, file.path(out, paste0("pfibox_pa_response_table_", tag, ".tsv")),
  sep = "\t", quote = FALSE, row.names = FALSE)
message("Retained ", nrow(response), " colony-condition observations.")
if (identical(Sys.getenv("PFIBOX_RESPONSE_ONLY"), "1")) quit(save = "no", status = 0)

# 3. Fit once, preserving the fit separately from inference. Each untreated
# source observation occurs once; repeated colony-position indices adjust
# reporter-wise fits and are not labelled independent biological replicates.
assay <- prepare_assay(response, reporter = "reporter", perturbation = "perturbation",
  control = "Control", response = "y_pa", batch = "dose_code_centered",
  plate = "drug_factor", replicate = "replicate", numeric_covariates = "dose_code_centered")
fit <- fit_destress(assay, technical = "replicate", empirical_bayes = TRUE,
  adjustment = "by_reporter", interaction = FALSE, verbose = TRUE)
saveRDS(fit, file.path(out, paste0("pfibox_pa_fit_", tag, ".rds")))

# 4. Extract both adjustment families. Figure 6 explicitly uses global BH
# via total_padj_global; the fit's reporting default does not select that plot.
res <- results(fit, verbose = TRUE)
res$response <- "pa_unscaled_integral_window"
res$t_start_h <- if (window == "all_available") NA_real_ else t_start
res$t_end_h <- if (window == "all_available") NA_real_ else t_end
res$time_window <- window
write.table(res, file.path(out, paste0("pfibox_pa_destress_pair_results_", tag, ".tsv")),
  sep = "\t", quote = FALSE, row.names = FALSE)
diagnostics <- perturbation_diagnostics(res, mean_effect = "total_effect",
  variance_effect = "total_effect", reporter = "reporter", perturbation = "perturbation", min_reporters = 10)
write.table(diagnostics, file.path(out, paste0("pfibox_pa_perturbation_diagnostics_", tag, ".tsv")),
  sep = "\t", quote = FALSE, row.names = FALSE)
run_metadata <- data.frame(run_tag = tag, response_summary = "unscaled_integral",
  time_window = window, pa_time_alignment = response$pa_time_alignment[1],
  reporter_limit = reporter_limit, perturbation_limit = perturbation_limit,
  n_reporters = length(unique(response$reporter)),
  n_noncontrol_perturbations = length(setdiff(unique(response$perturbation), "Control")),
  n_response_rows = nrow(response), n_result_pairs = nrow(res),
  requested_t_start_h = t_start, requested_t_end_h = t_end,
  realized_t_start_h = min(time_h[selected]), realized_t_end_h = max(time_h[selected]),
  pa_n_timepoints = sum(selected), design = "replicate", perturbation_target = "treatment",
  perturbation_column = "perturbation", model_path = "scalable", empirical_bayes = TRUE,
  technical = "replicate", response_source = source_file)
write.table(run_metadata, file.path(out, paste0("pfibox_pa_run_metadata_", tag, ".tsv")),
  sep = "\t", quote = FALSE, row.names = FALSE)
message("Saved response, fitted object, inference and diagnostics in ", out)
