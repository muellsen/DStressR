# Compact, attributed derived inputs for GitHub; the large source stays local.
source(file.path("analysis", "_helpers.R"))
input <- analysis_output_dir("french_pfibox", "temporal_unscaled_integral_analysis")
out <- analysis_path("data-raw", "french_pfibox", "derived")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
tag <- "reporters2000_perturbations100_Tall_available_treatment_replicate_scalable_temporal_unscaled_integral"
responses <- read.delim(file.path(input, paste0("pfibox_pa_response_table_", tag, ".tsv")))
layout_file <- file.path(input, "pfibox_colony_layout.tsv")
if (!file.exists(layout_file)) {
  e <- new.env(parent = emptyenv())
  message("Reading source colony layout for the derived-data export.")
  load(analysis_path("data-raw", "french_pfibox", "source", "gfpData.Rdata"), envir = e)
  layout <- as.data.frame(e$gfpData$legend)
  write.table(cbind(colony_position = seq_len(nrow(layout)), layout), layout_file,
    sep = "\t", quote = FALSE, row.names = FALSE)
  rm(e)
  invisible(gc())
}
layout <- read.delim(layout_file, check.names = FALSE)
stopifnot(nrow(responses) == 273792L, !anyDuplicated(layout$colony_position))
layout <- layout[match(responses$colony_position, layout$colony_position), ]
stopifnot(identical(as.character(layout$Promoter), as.character(responses$reporter)))
observation_id <- paste(responses$treatment, responses$colony_position, sep = ":")
stopifnot(!anyDuplicated(observation_id), all(is.finite(responses$y_pa)))
design <- data.frame(observation_id = observation_id, Promoter = responses$reporter,
  source_treatment = responses$treatment, Perturbation = responses$perturbation,
  antibiotic = responses$drug, source_dose_code = responses$concentration_code,
  is_reference = responses$is_control, colony_position = responses$colony_position,
  repeat_position = responses$replicate, stringsAsFactors = FALSE)
design <- cbind(design, layout[setdiff(names(layout), c("colony_position", "Promoter"))])
scalar <- data.frame(observation_id = observation_id, summary = "signed_pa_integral",
  window_start_h = responses$pa_window_start_h, window_end_h = responses$pa_window_end_h,
  n_timepoints = responses$pa_n_timepoints, value = responses$y_pa,
  units = "biomass_normalized_fluorescence_units", stringsAsFactors = FALSE)
write_gzip_tsv <- function(x, path) {
  con <- gzfile(path, "wt", compression = 9)
  on.exit(close(con))
  write.table(x, con, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}
paths <- file.path(out, c("experimental_design.tsv.gz", "scalar_responses.tsv.gz"))
write_gzip_tsv(design, paths[1])
write_gzip_tsv(scalar, paths[2])
manifest <- data.frame(file = basename(paths), rows = c(nrow(design), nrow(scalar)),
  bytes = file.info(paths)$size, md5 = unname(tools::md5sum(paths)),
  source_doi = "10.17632/fmv6vv4bsf.1", licence = "CC-BY-4.0",
  source_sha256 = "41b3f09bc7a8c3c435af9d03d0e33cede018ba7a5fd0ec900cf6f7d07079ed63")
write.table(manifest, file.path(out, "manifest.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
stopifnot(identical(read.delim(paths[1])$observation_id, read.delim(paths[2])$observation_id),
  isTRUE(all.equal(read.delim(paths[2])$value, responses$y_pa, tolerance = 1e-12)))
print(manifest[, 1:4])
