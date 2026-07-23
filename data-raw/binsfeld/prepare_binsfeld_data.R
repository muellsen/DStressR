# Build the Binsfeld et al. reporter-screen data shipped with DStressR.
#
# Sources:
# - PLOS Biology article: https://doi.org/10.1371/journal.pbio.3003260
# - S3 Data: reporter OD/LUX AUC table
# - S4 Data: author score/Z-score table
# - Zenodo code/data archive: https://doi.org/10.5281/zenodo.15600688
#
# The package ships the AUC-level reporter data rather than the larger plate
# reader time courses. S1/S2 time courses and the original analysis code remain
# public at the sources above and are sufficient to regenerate these AUCs.

stopifnot(requireNamespace("readxl", quietly = TRUE))

dir.create("data", showWarnings = FALSE)

download_plos_supplement <- function(suffix, destfile) {
  url <- paste0(
    "https://journals.plos.org/plosbiology/article/file?",
    "type=supplementary&id=10.1371/journal.pbio.3003260.",
    suffix
  )
  utils::download.file(url, destfile = destfile, mode = "wb", quiet = FALSE)
}

raw_dir <- file.path("data-raw", "binsfeld", "source")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

s3_xlsx <- file.path(raw_dir, "binsfeld_plos_s3_auc.xlsx")
s4_xlsx <- file.path(raw_dir, "binsfeld_plos_s4_scores.xlsx")
if (!file.exists(s3_xlsx)) {
  download_plos_supplement("s009", s3_xlsx)
}
if (!file.exists(s4_xlsx)) {
  download_plos_supplement("s010", s4_xlsx)
}

auc <- as.data.frame(readxl::read_excel(s3_xlsx), stringsAsFactors = FALSE)
names(auc) <- c(
  "strain", "promoter", "replicate", "well", "drug", "concentration_index",
  "concentration_ug_ml", "od_auc", "lux_auc", "od_auc_per_lux_auc", "removed"
)
auc$promoter[auc$promoter == "EVC (empty vector control)"] <- "EVC"
auc$compound <- ifelse(grepl("^Water_", auc$drug), "Water", auc$drug)
auc$replicate <- as.integer(auc$replicate)
auc$concentration_index <- as.integer(auc$concentration_index)
auc$od_auc <- as.numeric(auc$od_auc)
auc$lux_auc <- as.numeric(auc$lux_auc)
auc$od_auc_per_lux_auc <- as.numeric(auc$od_auc_per_lux_auc)
auc$concentration_ug_ml <- as.numeric(auc$concentration_ug_ml)
auc <- auc[, c(
  "strain", "promoter", "replicate", "well", "drug", "compound",
  "concentration_index", "concentration_ug_ml", "od_auc", "lux_auc",
  "od_auc_per_lux_auc", "removed"
)]

scores <- as.data.frame(readxl::read_excel(s4_xlsx), stringsAsFactors = FALSE)
names(scores) <- make.names(names(scores))
score_cols <- grep("_(1|2)$", names(scores), value = TRUE)
score_rows <- lapply(score_cols, function(col) {
  data.frame(
    well = scores$Well,
    drug = scores$Drug,
    compound = ifelse(grepl("^Water_", scores$Drug), "Water", scores$Drug),
    concentration_ug_ml = as.numeric(scores$Conc),
    strain = scores$Strain,
    statistic = scores$Variable,
    promoter = sub("_[12]$", "", col),
    replicate = as.integer(sub("^.*_", "", col)),
    value = suppressWarnings(as.numeric(scores[[col]])),
    stringsAsFactors = FALSE
  )
})
binsfeld_reporter_scores <- do.call(rbind, score_rows)

binsfeld_reporter_auc <- auc

save(
  binsfeld_reporter_auc,
  binsfeld_reporter_scores,
  file = file.path("data", "binsfeld_reporter_data.rda"),
  compress = "xz"
)
