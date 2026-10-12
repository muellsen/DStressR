#!/usr/bin/env Rscript

# Download and parse the public PFIbox GitHub sample files. The sample contains
# one untreated time-course experiment, so it validates data extraction but does
# not by itself define a perturbation screen.

source(file.path("analysis", "_helpers.R"))

source_dir <- analysis_path("data-raw", "french_pfibox", "source", "github_sample")
processed_dir <- analysis_path("data-raw", "french_pfibox", "processed")
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

download_if_missing <- function(url, path) {
  if (!file.exists(path)) {
    utils::download.file(url, path, mode = "wb", quiet = FALSE)
  }
  invisible(path)
}

base_url <- "https://raw.githubusercontent.com/sfrench007/pfibox/master/Analysis%20scripts"
timecourse_file <- download_if_missing(
  paste0(base_url, "/timecourse_results.csv"),
  file.path(source_dir, "timecourse_results.csv")
)
legend_file <- download_if_missing(
  paste0(base_url, "/legend.txt"),
  file.path(source_dir, "legend.txt")
)

legend <- utils::read.delim(
  legend_file,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)
raw <- utils::read.csv(
  timecourse_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

extract_timepoint <- function(label) {
  as.integer(sub(".*:([0-9]{1,4})$", "\\1", label))
}

raw$timepoint <- extract_timepoint(raw$Label)
raw$position <- ave(raw$timepoint, raw$timepoint, FUN = seq_along)
raw$reporter <- legend$Gene[raw$position]
raw$plate6144 <- paste0(legend$Row6144[raw$position], legend$Column6144[raw$position])
raw$time_min <- (raw$timepoint - min(raw$timepoint, na.rm = TRUE)) * 5
raw$time_h <- raw$time_min / 60
raw$perturbation <- "Untreated"

keep <- is.finite(raw$timepoint) & !is.na(raw$reporter) & nzchar(raw$reporter)
raw <- raw[keep, c(
  "perturbation", "reporter", "plate6144", "position", "timepoint", "time_min",
  "time_h", "IntDen", "RawIntDen"
)]

compact_file <- file.path(processed_dir, "pfibox_github_sample_timecourse_long.tsv")
utils::write.table(
  raw,
  compact_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

summary_tab <- aggregate(
  RawIntDen ~ perturbation + reporter,
  data = raw[raw$time_min >= max(raw$time_min, na.rm = TRUE) - 60, ],
  FUN = function(x) c(mean = mean(x, na.rm = TRUE), sd = stats::sd(x, na.rm = TRUE), n = sum(is.finite(x)))
)
summary_out <- data.frame(
  perturbation = summary_tab$perturbation,
  reporter = summary_tab$reporter,
  endpoint_mean_raw_intensity = summary_tab$RawIntDen[, "mean"],
  endpoint_sd_raw_intensity = summary_tab$RawIntDen[, "sd"],
  n_colony_timepoints = summary_tab$RawIntDen[, "n"],
  stringsAsFactors = FALSE
)
summary_file <- file.path(processed_dir, "pfibox_github_sample_endpoint_summary.tsv")
utils::write.table(
  summary_out,
  summary_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Wrote sample long table: ", compact_file)
message("Wrote sample endpoint summary: ", summary_file)

