#' Prepare a chemical-genomics assay table
#'
#' Computes a growth-adjusted log2 promoter-activity response from luminescence
#' and growth measurements. By default, promoter-specific growth exponents are
#' estimated from control wells with available technical-factor adjustment and
#' shrunk toward a global control-well slope.
#' Set `growth_exponent = 1` to reproduce the current workflow's log2(LUX / OD)
#' response.
#'
#' @param data A data frame with one row per promoter-compound-replicate well.
#' @param promoter,compound Column names identifying promoter and compound.
#' @param control Label in `compound` for the negative control, usually DMSO.
#' @param lux,growth Column names for luminescence and growth summaries.
#' @param growth_exponent Fixed coefficient for growth normalization, a named
#'   vector keyed by promoter, or `"estimate"` to estimate promoter-specific
#'   exponents from controls.
#' @param control_values Values in `compound` used as controls for growth
#'   exponent estimation. Defaults to `control`.
#' @param response Optional existing response column. If supplied, `lux` and
#'   `growth` are not used to compute the response.
#' @param batch,plate,replicate Optional technical-factor column names.
#'   When `growth_exponent = "estimate"`, these columns are also used as
#'   covariates while estimating promoter-specific growth exponents unless
#'   `growth_covariates` is supplied.
#' @param growth_covariates Optional technical covariate column names used only
#'   while estimating promoter-specific growth exponents from control wells. If
#'   `NULL`, DStressR uses the supplied `batch`, `plate`, and `replicate`
#'   columns for backwards compatibility.
#' @param numeric_covariates Optional subset of technical covariate column names
#'   that should remain numeric in model matrices. Other optional covariates are
#'   converted to factors.
#' @param background_promoter Optional reporter promoter used as a background
#'   reference, e.g. an Empty Vector Control. When supplied, the default
#'   background method is `"huber"`. The background reporter is matched to other
#'   promoters by `background_by`, the response is calibrated, and the
#'   background reporter is excluded from model-based testing.
#' @param background_method One of `"none"`, `"subtract"`, `"lm"`, or
#'   `"huber"`. If omitted, DStressR uses `"none"` when no
#'   `background_promoter` is supplied and `"huber"` when one is supplied.
#'   `"subtract"` subtracts the matched background response. `"lm"` and
#'   `"huber"` replace each non-background promoter response by residuals from
#'   a promoter-wise calibration against the matched background.
#' @param background_by Columns used to match each observation to the
#'   background reporter. Defaults to `compound` plus the supplied technical
#'   columns.
#' @param pseudocount Added before log2 transformation.
#' @return A `destress_assay` data frame.
#' @export
prepare_assay <- function(data,
                          promoter = "promoter",
                          compound = "compound",
                          control = "DMSO",
                          lux = "lux",
                          growth = "growth",
                          growth_exponent = "estimate",
                          control_values = control,
                          response = NULL,
                          batch = NULL,
                          plate = NULL,
                          replicate = NULL,
                          growth_covariates = NULL,
                          numeric_covariates = NULL,
                          background_promoter = NULL,
                          background_method = c("none", "subtract", "lm", "huber"),
                          background_by = NULL,
                          pseudocount = 1e-8) {
  stopifnot(is.data.frame(data))
  required <- c(promoter, compound)
  if (!is.null(response)) {
    required <- c(required, response)
  } else {
    required <- c(required, lux, growth)
  }
  optional <- unique(c(batch, plate, replicate, background_by, growth_covariates))
  missing_cols <- setdiff(c(required, optional), names(data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  numeric_covariates <- numeric_covariates[!is.na(numeric_covariates) & nzchar(numeric_covariates)]
  missing_numeric_covariates <- setdiff(numeric_covariates, names(data))
  if (length(missing_numeric_covariates) > 0) {
    stop("Unknown numeric covariate columns: ",
         paste(missing_numeric_covariates, collapse = ", "), call. = FALSE)
  }

  out <- data
  out$.promoter <- factor(out[[promoter]])
  out$.compound <- factor(out[[compound]])
  if (!control %in% levels(out$.compound)) {
    stop("Control compound '", control, "' was not found in `", compound, "`.", call. = FALSE)
  }
  out$.compound <- stats::relevel(out$.compound, ref = control)

  if (!is.null(response)) {
    out$.response <- as.numeric(out[[response]])
  } else {
    lux_value <- as.numeric(out[[lux]])
    growth_value <- as.numeric(out[[growth]])
    if (any(lux_value + pseudocount <= 0, na.rm = TRUE)) {
      stop("Luminescence values must be positive after adding pseudocount.", call. = FALSE)
    }
    if (any(growth_value + pseudocount <= 0, na.rm = TRUE)) {
      stop("Growth values must be positive after adding pseudocount.", call. = FALSE)
    }
    if (is.character(growth_exponent) && identical(growth_exponent, "estimate")) {
      if (is.null(growth_covariates)) {
        growth_covariates <- c(batch, plate, replicate)
      }
      growth_fit <- estimate_growth_exponents(
        out,
        promoter = promoter,
        compound = compound,
        lux = lux,
        growth = growth,
        covariates = growth_covariates,
        numeric_covariates = numeric_covariates,
        controls = control_values,
        pseudocount = pseudocount
      )
      alpha <- growth_fit$alpha_shrunk[match(as.character(out[[promoter]]), growth_fit$promoter)]
    } else if (length(growth_exponent) == 1) {
      growth_fit <- NULL
      alpha <- rep(as.numeric(growth_exponent), nrow(out))
    } else {
      growth_fit <- NULL
      if (is.null(names(growth_exponent))) {
        stop("Vector `growth_exponent` must be named by promoter.", call. = FALSE)
      }
      alpha <- as.numeric(growth_exponent[as.character(out[[promoter]])])
    }
    if (any(!is.finite(alpha))) {
      stop("Could not assign a finite growth exponent to every row.", call. = FALSE)
    }
    out$.growth_exponent <- alpha
    out$.response <- log2(lux_value + pseudocount) -
      alpha * log2(growth_value + pseudocount)
  }

  numeric_optional <- intersect(optional, numeric_covariates)
  for (nm in numeric_optional) {
    out[[nm]] <- as.numeric(out[[nm]])
  }
  for (nm in setdiff(optional, numeric_covariates)) {
    out[[nm]] <- factor(out[[nm]])
  }

  if (missing(background_method)) {
    background_method <- if (is.null(background_promoter)) "none" else "huber"
  }
  background_method <- normalize_background_method(background_method)
  background_fit <- NULL
  if (!is.null(background_promoter) && !identical(background_method, "none")) {
    background_promoter <- as.character(background_promoter)
    if (length(background_promoter) != 1 || is.na(background_promoter) || !nzchar(background_promoter)) {
      stop("`background_promoter` must be one promoter label.", call. = FALSE)
    }
    if (!background_promoter %in% levels(out$.promoter)) {
      stop("Background promoter '", background_promoter, "' was not found in `", promoter, "`.", call. = FALSE)
    }
    if (is.null(background_by)) {
      background_by <- unique(c(compound, batch, plate, replicate))
      background_by <- background_by[!is.na(background_by) & nzchar(background_by)]
    }
    background <- calibrate_background_response(
      out,
      promoter = promoter,
      background_promoter = background_promoter,
      method = background_method,
      by = background_by
    )
    out <- background$data
    background_fit <- background$fit
  } else {
    background_promoter <- NULL
    background_method <- "none"
    background_by <- character()
  }

  attr(out, "destress") <- list(
    promoter = promoter,
    compound = compound,
    control = control,
    lux = lux,
    growth = growth,
    control_values = control_values,
    response = response,
    batch = batch,
    plate = plate,
    replicate = replicate,
    growth_covariates = growth_covariates,
    numeric_covariates = numeric_covariates,
    growth_exponent = growth_exponent,
    growth_exponent_fit = if (exists("growth_fit")) growth_fit else NULL,
    background_promoter = background_promoter,
    background_method = background_method,
    background_by = background_by,
    background_fit = background_fit
  )
  class(out) <- c("destress_assay", class(out))
  out
}

normalize_background_method <- function(method = c("none", "subtract", "lm", "huber")) {
  method <- match.arg(method)
  method
}

calibration_key <- function(data, by) {
  if (length(by) == 0) {
    return(rep("all", nrow(data)))
  }
  do.call(paste, c(lapply(data[by], as.character), sep = "\r"))
}

calibrate_background_response <- function(data,
                                          promoter,
                                          background_promoter,
                                          method,
                                          by) {
  bg <- data[as.character(data[[promoter]]) == background_promoter, , drop = FALSE]
  if (nrow(bg) == 0) {
    stop("No rows found for `background_promoter = \"", background_promoter, "\"`.", call. = FALSE)
  }

  bg_key <- calibration_key(bg, by)
  bg_mean <- stats::aggregate(
    bg$.response,
    by = as.list(bg[by]),
    FUN = function(x) mean(x, na.rm = TRUE)
  )
  names(bg_mean)[ncol(bg_mean)] <- ".background_response"
  data$.background_key <- calibration_key(data, by)
  bg_mean$.background_key <- calibration_key(bg_mean, by)
  lookup <- stats::setNames(bg_mean$.background_response, bg_mean$.background_key)
  data$.background_response <- as.numeric(lookup[data$.background_key])
  data$.response_uncalibrated <- data$.response
  is_background <- as.character(data[[promoter]]) == background_promoter

  fit_rows <- data.frame(
    promoter = character(),
    method = character(),
    n = integer(),
    intercept = numeric(),
    slope = numeric(),
    stringsAsFactors = FALSE
  )

  if (identical(method, "subtract")) {
    data$.response[!is_background] <- data$.response[!is_background] -
      data$.background_response[!is_background]
    fit_rows <- data.frame(
      promoter = setdiff(unique(as.character(data[[promoter]])), background_promoter),
      method = "subtract",
      n = NA_integer_,
      intercept = 0,
      slope = 1,
      stringsAsFactors = FALSE
    )
  } else {
    if (identical(method, "huber") && !requireNamespace("MASS", quietly = TRUE)) {
      stop("Package `MASS` is required for `background_method = \"huber\"`.", call. = FALSE)
    }
    for (p in setdiff(unique(as.character(data[[promoter]])), background_promoter)) {
      idx <- as.character(data[[promoter]]) == p
      complete <- idx & is.finite(data$.response) & is.finite(data$.background_response)
      if (sum(complete) < 3) {
        stop("Need at least three matched background observations for promoter '", p, "'.", call. = FALSE)
      }
      d <- data.frame(
        response = data$.response[complete],
        background = data$.background_response[complete]
      )
      fit <- if (identical(method, "huber")) {
        MASS::rlm(response ~ background, data = d)
      } else {
        stats::lm(response ~ background, data = d)
      }
      coef <- stats::coef(fit)
      pred_idx <- idx & is.finite(data$.background_response)
      data$.response[pred_idx] <- data$.response[pred_idx] -
        (unname(coef[[1]]) + unname(coef[[2]]) * data$.background_response[pred_idx])
      data$.response[idx & !is.finite(data$.background_response)] <- NA_real_
      fit_rows <- rbind(
        fit_rows,
        data.frame(
          promoter = p,
          method = method,
          n = sum(complete),
          intercept = unname(coef[[1]]),
          slope = unname(coef[[2]]),
          stringsAsFactors = FALSE
        )
      )
    }
  }

  data$.background_key <- NULL
  list(data = data, fit = fit_rows[order(fit_rows$promoter), , drop = FALSE])
}

#' Read exported Campylobacter expression data
#'
#' @param expression_file Path to `expression_values.tsv.gz`.
#' @param libmap_file Path to `LibMap.txt`.
#' @return A data frame with `srn_code` and `ProductName` joined in.
#' @export
read_campylobacter_expression <- function(expression_file, libmap_file) {
  expr <- utils::read.delim(expression_file, check.names = FALSE)
  libmap <- utils::read.delim(libmap_file, check.names = FALSE)
  libmap$lib_plate <- paste0("lp", libmap[["Library plate"]])
  libmap$srn_code <- paste(libmap$lib_plate, libmap[["Well"]], sep = "_")
  if (!"srn_code" %in% names(expr)) {
    if (!all(c("libplate", "well") %in% names(expr))) {
      stop("Expression file needs either `srn_code` or `libplate` + `well`.", call. = FALSE)
    }
    expr$srn_code <- paste(expr$libplate, expr$well, sep = "_")
  }
  merge(expr, libmap[, c("srn_code", "ProductName", "Catalog Number")],
        by = "srn_code", all.x = TRUE)
}
