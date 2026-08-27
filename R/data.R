#' Binsfeld et al. reporter screen AUC data
#'
#' A public *E. coli* reporter-screen data set from Binsfeld et al. (2025),
#' prepared as an AUC-level long table for DStressR examples and tests. The
#' rows are reporter/strain/replicate/well observations from the PLOS Biology
#' S3 Data supplement. `drug` keeps the original Binsfeld perturbation labels,
#' including `Water_1` and `Water_2`. `compound` is a derived grouping column
#' that collapses the water control wells to `Water` for summaries.
#' `dose_level` is derived from `concentration_index` so that larger values
#' correspond to higher compound concentration.
#'
#' The source article is https://doi.org/10.1371/journal.pbio.3003260. The
#' associated Zenodo code/data archive is https://doi.org/10.5281/zenodo.15600688.
#'
#' @format A data frame with 24,576 rows and 13 columns:
#' \describe{
#'   \item{strain}{Reporter host strain.}
#'   \item{reporter}{Reporter label, including `EVC`.}
#'   \item{replicate}{Reporter replicate number.}
#'   \item{well}{384-well plate coordinate.}
#'   \item{drug}{Original drug/control label from the source table.}
#'   \item{compound}{Derived drug/control grouping label, with water controls
#'   collapsed.}
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
#' @format A data frame with one row per well, strain, statistic, reporter, and
#'   replicate. Columns are `well`, `drug`, `compound`, `concentration_ug_ml`,
#'   `strain`, `statistic`, `reporter`, `replicate`, and `value`.
#' @source Binsfeld et al. (2025), PLOS Biology, S4 Data.
"binsfeld_reporter_scores"

#' Zaslaver et al. fixed-growth-rate promoter activity data
#'
#' Public *E. coli* promoter-activity data from Zaslaver et al. (2009),
#' prepared from Dataset S2 of the PLOS Computational Biology article. Each row
#' is one reporter/condition/growth-rate summary. The `reporter` column keeps
#' the original source label, while `reporter_id` combines the source row index
#' and label because several reporter labels occur more than once in the source
#' workbooks and two reporter labels vary across condition sheets.
#'
#' The source article is https://doi.org/10.1371/journal.pcbi.1000545. The
#' supplementary Excel files were obtained from the CaltechAUTHORS record
#' https://authors.library.caltech.edu/records/g76ja-xga32.
#'
#' @format A data frame with 23,040 rows and 10 columns:
#' \describe{
#'   \item{reporter_index}{Source row index of the reporter construct.}
#'   \item{reporter_id}{Unique reporter identifier derived from row index and
#'   source row index.}
#'   \item{reporter}{Original source reporter label.}
#'   \item{growth_rate}{Fixed growth rate at which promoter activity was
#'   summarized, in divisions per hour.}
#'   \item{condition}{Machine-readable condition label.}
#'   \item{condition_label}{Source condition label.}
#'   \item{promoter_activity}{Promoter activity for the reporter under the
#'   condition at the fixed growth rate.}
#'   \item{mean_promoter_activity}{Mean promoter activity across conditions at
#'   the fixed growth rate.}
#'   \item{sd_promoter_activity}{Standard deviation of promoter activity across
#'   conditions at the fixed growth rate.}
#'   \item{cv_promoter_activity}{Coefficient of variation of promoter activity
#'   across conditions at the fixed growth rate.}
#' }
#' @source Zaslaver et al. (2009), PLOS Computational Biology, Dataset S2.
"zaslaver_promoter_activity"

#' Zaslaver et al. promoter activity and OD time courses
#'
#' Long-form version of Dataset S1 from Zaslaver et al. (2009). The source file
#' contains, for each growth condition, a promoter-activity matrix followed by
#' the matching OD matrix. These values are author-processed promoter activity
#' and OD summaries, not the raw fluorescence plate-reader traces.
#'
#' @format A data frame with 622,080 rows and 9 columns:
#' \describe{
#'   \item{condition}{Machine-readable condition label.}
#'   \item{condition_label}{Source condition label.}
#'   \item{reporter_index}{Source row index of the reporter construct.}
#'   \item{reporter_id}{Unique reporter identifier derived from row index and
#'   source row index.}
#'   \item{reporter}{Original source reporter label.}
#'   \item{time_index}{Source measurement-column index.}
#'   \item{time_min}{Nominal time in minutes, assuming 16-minute sampling
#'   intervals as reported by Zaslaver et al.}
#'   \item{promoter_activity}{Author-processed promoter activity.}
#'   \item{od}{Optical-density value from the matching OD block.}
#' }
#' @source Zaslaver et al. (2009), PLOS Computational Biology, Dataset S1.
"zaslaver_promoter_timecourse"

#' Zaslaver et al. reporter annotation classes
#'
#' Reporter annotation class matrix prepared from Dataset S3 of Zaslaver et al.
#' (2009). Annotation columns are binary indicators from the source workbook.
#' As in [zaslaver_promoter_activity], `reporter_id` is the unique reporter
#' construct identifier and `reporter` preserves the source label.
#'
#' @format A data frame with 1,920 rows and 249 columns. The first columns are
#'   `reporter_index`, `reporter_id`, and `reporter`; remaining columns are
#'   binary annotation indicators.
#' @source Zaslaver et al. (2009), PLOS Computational Biology, Dataset S3.
"zaslaver_promoter_annotations"
