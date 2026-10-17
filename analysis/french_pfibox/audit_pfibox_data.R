#!/usr/bin/env Rscript

# Audit the deposited French et al. PFIbox gfpData.Rdata object for DStressR
# suitability. This script records treatment structure, reporter replication,
# and matrix availability without running the full model.

source(file.path("analysis", "_helpers.R"))

source_file <- analysis_path("data-raw", "french_pfibox", "source", "gfpData.Rdata")
processed_dir <- analysis_path("data-raw", "french_pfibox", "processed")
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(source_file)) {
  stop(
    "Missing PFIbox object: ", source_file, "\n",
    "Download the full deposited object from Mendeley Data ",
    "doi:10.17632/fmv6vv4bsf.1.",
    call. = FALSE
  )
}

load(source_file)
if (!exists("gfpData") || !is.list(gfpData) || is.null(gfpData$legend)) {
  stop("Expected a list object `gfpData` with a `legend` element.", call. = FALSE)
}

legend <- as.data.frame(gfpData$legend, stringsAsFactors = FALSE)
treatments <- setdiff(names(gfpData), "legend")

parse_treatment <- function(x) {
  m <- regexec("^(.*)-([0-9]+)$", x)
  r <- regmatches(x, m)[[1]]
  if (length(r) == 3) {
    data.frame(
      treatment = x,
      drug = r[2],
      concentration_code = r[3],
      is_control = identical(r[3], "0"),
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      treatment = x,
      drug = x,
      concentration_code = NA_character_,
      is_control = FALSE,
      stringsAsFactors = FALSE
    )
  }
}

treatment_metadata <- do.call(rbind, lapply(treatments, parse_treatment))
treatment_metadata$analysis_perturbation <- ifelse(
  treatment_metadata$is_control,
  "Control",
  treatment_metadata$treatment
)

reporter_counts <- as.data.frame(table(legend$Promoter), stringsAsFactors = FALSE)
names(reporter_counts) <- c("reporter", "n_colony_positions")
reporter_counts <- reporter_counts[order(-reporter_counts$n_colony_positions, reporter_counts$reporter), ]

matrix_audit <- do.call(rbind, lapply(treatments, function(nm) {
  x <- gfpData[[nm]]
  data.frame(
    treatment = nm,
    has_raw = is.matrix(x$raw),
    raw_rows = if (is.matrix(x$raw)) nrow(x$raw) else NA_integer_,
    raw_cols = if (is.matrix(x$raw)) ncol(x$raw) else NA_integer_,
    has_raw_loess = is.matrix(x$raw.loess),
    has_constit_normalized = is.matrix(x$constit.normalized),
    has_promoter_activity = is.matrix(x$promoter.activity),
    promoter_activity_rows = if (is.matrix(x$promoter.activity)) nrow(x$promoter.activity) else NA_integer_,
    promoter_activity_cols = if (is.matrix(x$promoter.activity)) ncol(x$promoter.activity) else NA_integer_,
    raw_finite_fraction = if (is.matrix(x$raw)) mean(is.finite(x$raw)) else NA_real_,
    constit_normalized_finite_fraction =
      if (is.matrix(x$constit.normalized)) mean(is.finite(x$constit.normalized)) else NA_real_,
    promoter_activity_finite_fraction =
      if (is.matrix(x$promoter.activity)) mean(is.finite(x$promoter.activity)) else NA_real_,
    stringsAsFactors = FALSE
  )
}))

summary_tab <- data.frame(
  item = c(
    "treatments",
    "drugs",
    "control_treatments",
    "reporters_including_controls",
    "non_control_reporters",
    "positions",
    "raw_timepoints",
    "promoter_activity_timepoints"
  ),
  value = c(
    length(treatments),
    length(unique(treatment_metadata$drug)),
    sum(treatment_metadata$is_control),
    length(unique(legend$Promoter)),
    length(setdiff(unique(legend$Promoter), c("Empty", "U66", "U139"))),
    nrow(legend),
    unique(matrix_audit$raw_cols)[1],
    unique(matrix_audit$promoter_activity_cols)[1]
  )
)

utils::write.table(
  treatment_metadata,
  file.path(processed_dir, "pfibox_treatment_metadata.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
utils::write.table(
  reporter_counts,
  file.path(processed_dir, "pfibox_reporter_position_counts.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
utils::write.table(
  matrix_audit,
  file.path(processed_dir, "pfibox_object_audit.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
utils::write.table(
  summary_tab,
  file.path(processed_dir, "pfibox_audit_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

print(summary_tab)
message("Wrote PFIbox audit tables to: ", processed_dir)

