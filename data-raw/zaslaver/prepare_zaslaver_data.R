# Build the Zaslaver et al. promoter-activity data shipped with DStressR.
#
# Sources:
# - PLOS Computational Biology article: https://doi.org/10.1371/journal.pcbi.1000545
# - Dataset S1: promoter activity and OD by condition
# - Dataset S2: promoter activity summaries at two fixed growth rates
# - Dataset S3: promoter annotation classes
#
# The package data preserve the source reporter labels and condition names. The
# full S1 table is already processed by the original authors into promoter
# activity and OD matrices; it is not the raw plate-reader fluorescence signal.

stopifnot(requireNamespace("readxl", quietly = TRUE))

dir.create("data", showWarnings = FALSE)

source_dir <- file.path("data-raw", "zaslaver", "source")
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

download_caltech_supplement <- function(filename, destfile) {
  url <- paste0(
    "https://authors.library.caltech.edu/records/g76ja-xga32/files/",
    filename,
    "?download=1"
  )
  utils::download.file(url, destfile = destfile, mode = "wb", quiet = FALSE)
}

s1_xls <- file.path(source_dir, "journal.pcbi.1000545.s021.xls")
s2_xls <- file.path(source_dir, "journal.pcbi.1000545.s022.xls")
s3_xls <- file.path(source_dir, "journal.pcbi.1000545.s023.xls")

for (path in c(s1_xls, s2_xls, s3_xls)) {
  if (!file.exists(path)) {
    download_caltech_supplement(basename(path), path)
  }
}

clean_number <- function(x) {
  suppressWarnings(as.numeric(x))
}

make_condition <- function(x) {
  out <- tolower(x)
  out <- gsub("[^a-z0-9]+", "_", out)
  out <- gsub("^_|_$", "", out)
  out
}

parse_condition_sheet <- function(path, sheet) {
  raw <- suppressWarnings(
    as.data.frame(readxl::read_excel(path, sheet = sheet, col_names = FALSE),
      stringsAsFactors = FALSE
    )
  )
  first_col <- as.character(raw[[1]])
  pa_start <- which(first_col == "PA")
  od_start <- which(first_col == "OD")
  stopifnot(length(pa_start) == 1L, length(od_start) == 1L)

  pa_rows <- seq(pa_start + 2L, od_start - 1L)
  od_rows <- seq(od_start + 2L, nrow(raw))
  stopifnot(length(pa_rows) == length(od_rows))

  reporters <- as.character(raw[[1]][pa_rows])
  od_reporters <- as.character(raw[[1]][od_rows])
  stopifnot(identical(reporters, od_reporters))

  value_cols <- seq.int(2L, ncol(raw))
  pa <- as.matrix(raw[pa_rows, value_cols])
  od <- as.matrix(raw[od_rows, value_cols])
  storage.mode(pa) <- "numeric"
  storage.mode(od) <- "numeric"

  n_reporter <- length(reporters)
  n_time <- length(value_cols)
  idx <- expand.grid(
    reporter_index = seq_len(n_reporter),
    time_index = seq_len(n_time)
  )

  data.frame(
    condition = make_condition(sheet),
    condition_label = sheet,
    reporter_index = idx$reporter_index,
    reporter_id = sprintf("reporter_%04d", idx$reporter_index),
    reporter = reporters[idx$reporter_index],
    time_index = idx$time_index,
    time_min = (idx$time_index - 1L) * 16,
    promoter_activity = as.vector(pa),
    od = as.vector(od),
    stringsAsFactors = FALSE
  )
}

condition_sheets <- readxl::excel_sheets(s1_xls)
zaslaver_promoter_timecourse <- do.call(
  rbind,
  lapply(condition_sheets, function(sheet) parse_condition_sheet(s1_xls, sheet))
)

summary_raw <- suppressWarnings(
  as.data.frame(readxl::read_excel(s2_xls, sheet = "CV_table", col_names = FALSE),
    stringsAsFactors = FALSE
  )
)
summary_header <- as.character(unlist(summary_raw[1, ], use.names = FALSE))
summary_data <- summary_raw[-1, , drop = FALSE]

summary_specs <- list(
  list(growth_rate = 0.8, gene_col = 1L, mean_col = 2L, sd_col = 3L, cv_col = 4L, first_pa_col = 5L, last_pa_col = 10L),
  list(growth_rate = 0.25, gene_col = 11L, mean_col = 12L, sd_col = 13L, cv_col = 14L, first_pa_col = 15L, last_pa_col = 20L)
)

summary_long <- lapply(summary_specs, function(spec) {
  pa_cols <- seq.int(spec$first_pa_col, spec$last_pa_col)
  condition_labels <- sub("^PA ", "", summary_header[pa_cols])
  condition_labels <- sub(paste0(" at ", spec$growth_rate, "$"), "", condition_labels)
  condition_labels <- sub(paste0(" ", spec$growth_rate, "$"), "", condition_labels)
  condition_labels <- sub("phosphate limitation", "Phosphate limited", condition_labels)
  condition_labels <- sub("nitrogen limitation", "Nitrogen limited", condition_labels)
  condition_labels <- sub("^glucose$", "Glucose", condition_labels)
  condition_labels <- sub("^no glucose$", "no Glucose", condition_labels)
  condition_labels <- sub("^no AA$", "no AA", condition_labels)
  condition_labels <- sub("^ethanol$", "Ethanol", condition_labels)

  rows <- lapply(seq_along(pa_cols), function(k) {
    data.frame(
      reporter_index = seq_len(nrow(summary_data)),
      reporter_id = sprintf("reporter_%04d", seq_len(nrow(summary_data))),
      reporter = as.character(summary_data[[spec$gene_col]]),
      growth_rate = spec$growth_rate,
      condition = make_condition(condition_labels[k]),
      condition_label = condition_labels[k],
      promoter_activity = clean_number(summary_data[[pa_cols[k]]]),
      mean_promoter_activity = clean_number(summary_data[[spec$mean_col]]),
      sd_promoter_activity = clean_number(summary_data[[spec$sd_col]]),
      cv_promoter_activity = clean_number(summary_data[[spec$cv_col]]),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
})
zaslaver_promoter_activity <- do.call(rbind, summary_long)

annotation_raw <- suppressWarnings(
  as.data.frame(readxl::read_excel(s3_xls, sheet = "annotations", col_names = FALSE),
    stringsAsFactors = FALSE
  )
)
annotation_names <- as.character(unlist(annotation_raw[1, ], use.names = FALSE))
annotation_names <- make.names(annotation_names, unique = TRUE)
zaslaver_promoter_annotations <- annotation_raw[-1, , drop = FALSE]
names(zaslaver_promoter_annotations) <- annotation_names
names(zaslaver_promoter_annotations)[1] <- "reporter"
zaslaver_promoter_annotations$reporter_index <- seq_len(nrow(zaslaver_promoter_annotations))
zaslaver_promoter_annotations$reporter_id <- sprintf("reporter_%04d", zaslaver_promoter_annotations$reporter_index)
zaslaver_promoter_annotations <- zaslaver_promoter_annotations[
  c("reporter_index", "reporter_id", setdiff(names(zaslaver_promoter_annotations), c("reporter_index", "reporter_id")))
]
for (j in seq.int(2L, ncol(zaslaver_promoter_annotations))) {
  if (!names(zaslaver_promoter_annotations)[j] %in% c("reporter_id", "reporter")) {
    zaslaver_promoter_annotations[[j]] <- as.integer(clean_number(zaslaver_promoter_annotations[[j]]))
  }
}

unlink(file.path("data", "zaslaver_promoter_data.rda"))
save(zaslaver_promoter_activity,
  file = file.path("data", "zaslaver_promoter_activity.rda"),
  compress = "xz"
)
save(zaslaver_promoter_timecourse,
  file = file.path("data", "zaslaver_promoter_timecourse.rda"),
  compress = "xz"
)
save(zaslaver_promoter_annotations,
  file = file.path("data", "zaslaver_promoter_annotations.rda"),
  compress = "xz"
)
