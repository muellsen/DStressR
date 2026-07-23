#' Binsfeld et al. reporter screen AUC data
#'
#' A public *E. coli* reporter-screen data set from Binsfeld et al. (2025),
#' prepared as an AUC-level long table for DStressR examples and tests. The
#' rows are promoter/strain/replicate/well observations from the PLOS Biology
#' S3 Data supplement. `compound` collapses the water control wells
#' (`Water_1`, `Water_2`) to `Water`; the original label remains in `drug`.
#' `dose_level` is derived from `concentration_index` so that larger values
#' correspond to higher compound concentration.
#'
#' The source article is https://doi.org/10.1371/journal.pbio.3003260. The
#' associated Zenodo code/data archive is https://doi.org/10.5281/zenodo.15600688.
#'
#' @format A data frame with 24,576 rows and 12 columns:
#' \describe{
#'   \item{strain}{Reporter host strain.}
#'   \item{promoter}{Reporter promoter, including `EVC`.}
#'   \item{replicate}{Reporter replicate number.}
#'   \item{well}{384-well plate coordinate.}
#'   \item{drug}{Original compound/control label.}
#'   \item{compound}{DStressR compound label, with water controls collapsed.}
#'   \item{concentration_index}{Dose-series index from the source table.}
#'   \item{dose_level}{Dose-oriented serial dilution level; larger values
#'   correspond to higher concentration.}
#'   \item{concentration_ug_ml}{Compound concentration in micrograms per ml.}
#'   \item{od_auc}{Optical-density area under the curve.}
#'   \item{lux_auc}{Luminescence area under the curve.}
#'   \item{od_auc_per_lux_auc}{Source-table OD/LUX AUC ratio.}
#'   \item{removed}{Author quality-control flag.}
#' }
#' @source Binsfeld et al. (2025), PLOS Biology, S3 Data.
"binsfeld_reporter_auc"

#' Binsfeld et al. reporter scores and Z-scores
#'
#' Long-form version of the PLOS Biology S4 Data supplement from Binsfeld et al.
#' (2025). These values reproduce the authors' score/Z-score hit-calling
#' workflow and can be compared with DStressR model-based calls from
#' [binsfeld_reporter_auc].
#'
#' @format A data frame with one row per well, strain, statistic, promoter, and
#'   replicate. Columns are `well`, `drug`, `compound`, `concentration_ug_ml`,
#'   `strain`, `statistic`, `promoter`, `replicate`, and `value`.
#' @source Binsfeld et al. (2025), PLOS Biology, S4 Data.
"binsfeld_reporter_scores"
