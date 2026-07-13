#' Add DGrowthR-derived growth metrics to a DStressR assay table
#'
#' DStressR is designed to work hand in hand with DGrowthR: DGrowthR models the
#' optical-density growth curves, and DStressR can use one chosen DGrowthR
#' growth metric as the growth covariate in [prepare_assay()].
#'
#' This function is deliberately opt-in. The current hit-determination scripts
#' use the original exported growth summaries unless the analyst explicitly
#' calls this helper and passes the resulting column to [prepare_assay()]. This
#' keeps the present analysis reproducible while making the DGrowthR bridge
#' available for future comparisons.
#'
#' This helper joins `object@growth_parameters` from a fitted DGrowthR object
#' onto a promoter/luminescence assay table. Run
#' `DGrowthR::estimate_growth_parameters()` first, optionally with
#' `od_auc_at_t` if you want time-specific columns such as `OD_16` or `AUC_16`.
#'
#' @param data A data frame with one row per DStressR assay observation.
#' @param object A DGrowthR object after `DGrowthR::estimate_growth_parameters()`.
#' @param by Column in `data` and `object@metadata` identifying the growth
#'   curve for each assay row. DGrowthR calls this `curve_id`.
#' @param model_covariate Metadata column used as DGrowthR's
#'   `estimate_growth_parameters(model_covariate = ...)`. Defaults to `by`,
#'   corresponding to one GP fit per growth curve.
#' @param growth_metric Column in `object@growth_parameters` to use as the
#'   growth measurement. Common choices are `max_growth`, `AUC`, `OD_16`, and
#'   `AUC_16`.
#' @param output Column name to create in `data`.
#' @param keep_dgrowthr_columns If `TRUE`, also keep the DGrowthR join key and
#'   the unrenamed metric column when possible.
#' @return `data` with an additional numeric growth column named by `output`.
#' @export
add_dgrowthr_growth <- function(data,
                                object,
                                by = "curve_id",
                                model_covariate = by,
                                growth_metric = "max_growth",
                                output = "growth",
                                keep_dgrowthr_columns = FALSE) {
  stopifnot(is.data.frame(data))
  if (!methods::is(object, "DGrowthR")) {
    stop("`object` must be a DGrowthR object.", call. = FALSE)
  }
  if (!by %in% names(data)) {
    stop("Column `", by, "` was not found in `data`.", call. = FALSE)
  }

  metadata <- methods::slot(object, "metadata")
  growth_parameters <- methods::slot(object, "growth_parameters")
  if (!is.data.frame(metadata) || nrow(metadata) == 0) {
    stop("DGrowthR object has an empty `metadata` slot.", call. = FALSE)
  }
  if (!is.data.frame(growth_parameters) || nrow(growth_parameters) == 0) {
    stop(
      "DGrowthR object has no growth parameters. Run ",
      "`DGrowthR::estimate_growth_parameters()` first.",
      call. = FALSE
    )
  }
  missing_metadata <- setdiff(c(by, model_covariate), names(metadata))
  if (length(missing_metadata) > 0) {
    stop(
      "DGrowthR metadata is missing required column(s): ",
      paste(missing_metadata, collapse = ", "),
      call. = FALSE
    )
  }
  missing_growth <- setdiff(c("gpfit_id", growth_metric), names(growth_parameters))
  if (length(missing_growth) > 0) {
    stop(
      "DGrowthR growth parameters are missing required column(s): ",
      paste(missing_growth, collapse = ", "),
      ". Available metrics include: ",
      paste(setdiff(names(growth_parameters), "gpfit_id"), collapse = ", "),
      call. = FALSE
    )
  }

  map <- unique(metadata[, c(by, model_covariate), drop = FALSE])
  names(map) <- c(".dstressr_curve_id", ".dstressr_gpfit_id")
  map$.dstressr_curve_id <- as.character(map$.dstressr_curve_id)
  map$.dstressr_gpfit_id <- as.character(map$.dstressr_gpfit_id)
  duplicated_curves <- duplicated(map$.dstressr_curve_id)
  if (any(duplicated_curves)) {
    stop("DGrowthR metadata maps some `", by, "` values to multiple GP fits.", call. = FALSE)
  }

  growth_map <- growth_parameters[, c("gpfit_id", growth_metric), drop = FALSE]
  names(growth_map) <- c(".dstressr_gpfit_id", ".dstressr_growth")
  growth_map$.dstressr_gpfit_id <- as.character(growth_map$.dstressr_gpfit_id)
  duplicated_gpfit <- duplicated(growth_map$.dstressr_gpfit_id)
  if (any(duplicated_gpfit)) {
    stop("DGrowthR growth parameters contain duplicated `gpfit_id` values.", call. = FALSE)
  }

  map <- merge(map, growth_map, by = ".dstressr_gpfit_id", all.x = TRUE, sort = FALSE)
  out <- data
  out$.dstressr_order <- seq_len(nrow(out))
  out$.dstressr_curve_id <- as.character(out[[by]])
  out <- merge(out, map, by = ".dstressr_curve_id", all.x = TRUE, sort = FALSE)
  out <- out[order(out$.dstressr_order), , drop = FALSE]
  if (any(!is.finite(out$.dstressr_growth))) {
    missing_n <- sum(!is.finite(out$.dstressr_growth))
    stop(
      "Could not assign finite DGrowthR growth values to ",
      missing_n,
      " row(s). Check `by`, `model_covariate`, and `growth_metric`.",
      call. = FALSE
    )
  }

  out[[output]] <- as.numeric(out$.dstressr_growth)
  if (!isTRUE(keep_dgrowthr_columns)) {
    out$.dstressr_order <- NULL
    out$.dstressr_curve_id <- NULL
    out$.dstressr_gpfit_id <- NULL
    out$.dstressr_growth <- NULL
  } else {
    out$.dstressr_order <- NULL
  }
  rownames(out) <- NULL
  out
}
