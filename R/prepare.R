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
#'   covariates while estimating promoter-specific growth exponents.
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
                          pseudocount = 1e-8) {
  stopifnot(is.data.frame(data))
  required <- c(promoter, compound)
  if (!is.null(response)) {
    required <- c(required, response)
  } else {
    required <- c(required, lux, growth)
  }
  optional <- c(batch, plate, replicate)
  missing_cols <- setdiff(c(required, optional), names(data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
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
      growth_covariates <- c(batch, plate, replicate)
      growth_fit <- estimate_growth_exponents(
        out,
        promoter = promoter,
        compound = compound,
        lux = lux,
        growth = growth,
        covariates = growth_covariates,
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

  for (nm in optional) {
    out[[nm]] <- factor(out[[nm]])
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
    growth_exponent = growth_exponent,
    growth_exponent_fit = if (exists("growth_fit")) growth_fit else NULL
  )
  class(out) <- c("destress_assay", class(out))
  out
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
