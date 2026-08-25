#' Estimate reporter-specific growth normalization exponents
#'
#' Estimates how luminescence scales with growth in negative-control wells,
#' optionally adjusting for technical factors such as batch, plate, or
#' replicate:
#'
#' \deqn{\log_2(LUX_i) = a_g + \alpha_g \log_2(growth_i) + X_i\theta + e_i}
#'
#' Raw reporter-specific slopes are then shrunk toward a global control-well
#' slope using an empirical-Bayes normal prior. The shrunken
#' \eqn{\alpha_g} values can be used in [prepare_assay()] to compute:
#'
#' \deqn{y_i = \log_2(LUX_i) - \alpha_g \log_2(growth_i)}
#'
#' @param data A data frame with one row per well.
#' @param reporter,perturbation,lux,growth Column names.
#' @param covariates Optional technical covariate column names to include as
#'   additive adjustment terms when estimating growth slopes. Only covariates
#'   with more than one observed value in the relevant control subset are used.
#' @param numeric_covariates Optional subset of `covariates` that should enter
#'   the growth-slope model as numeric variables. Other covariates are converted
#'   to factors.
#' @param controls Control values in `perturbation`, usually DMSO wells.
#' @param pseudocount Added before log2 transformation.
#' @param min_control_n Minimum control wells needed for a reporter-specific
#'   raw slope. Reporters with fewer controls use the global slope.
#' @param shrink If `TRUE`, shrink reporter-specific slopes toward the global
#'   control slope.
#' @param alpha_bounds Optional numeric length-2 bounds for the final exponent.
#'   Use `NULL` for no clipping.
#' @return A data frame with raw reporter intercepts and raw and shrunken growth
#'   exponents per reporter.
#' @export
estimate_growth_exponents <- function(data,
                                      reporter = "reporter",
                                      perturbation = "perturbation",
                                      lux = "lux",
                                      growth = "growth",
                                      covariates = NULL,
                                      numeric_covariates = NULL,
                                      controls = "DMSO",
                                      pseudocount = 1e-8,
                                      min_control_n = 8,
                                      shrink = TRUE,
                                      alpha_bounds = c(-2, 3)) {
  stopifnot(is.data.frame(data))
  required <- c(reporter, perturbation, lux, growth, covariates)
  missing_cols <- setdiff(required, names(data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  controls_df <- data[data[[perturbation]] %in% controls, , drop = FALSE]
  if (nrow(controls_df) < 3) {
    stop("Need at least three control rows to estimate growth exponent.", call. = FALSE)
  }

  controls_df$.log_lux <- log2(as.numeric(controls_df[[lux]]) + pseudocount)
  controls_df$.log_growth <- log2(as.numeric(controls_df[[growth]]) + pseudocount)
  controls_df <- controls_df[is.finite(controls_df$.log_lux) & is.finite(controls_df$.log_growth), ]
  if (nrow(controls_df) < 3) {
    stop("Control rows have too few finite log-transformed values.", call. = FALSE)
  }

  fit_growth_slope <- function(d, adjustment_cols) {
    usable_covariates <- adjustment_cols[vapply(adjustment_cols, function(nm) {
      length(unique(stats::na.omit(d[[nm]]))) > 1
    }, logical(1))]
    rhs <- c(".log_growth", usable_covariates)
    form <- stats::as.formula(paste(".log_lux ~", paste(rhs, collapse = " + ")))
    fit <- stats::lm(form, data = d)
    coefs <- summary(fit)$coefficients
    estimate <- NA_real_
    se <- NA_real_
    intercept <- NA_real_
    intercept_se <- NA_real_
    if (".log_growth" %in% rownames(coefs)) {
      estimate <- unname(coefs[".log_growth", "Estimate"])
      se <- unname(coefs[".log_growth", "Std. Error"])
    }
    if ("(Intercept)" %in% rownames(coefs)) {
      intercept <- unname(coefs["(Intercept)", "Estimate"])
      intercept_se <- unname(coefs["(Intercept)", "Std. Error"])
    }
    list(
      estimate = estimate,
      se = se,
      intercept = intercept,
      intercept_se = intercept_se,
      df = stats::df.residual(fit),
      covariates = paste(usable_covariates, collapse = ";")
    )
  }

  covariates <- unique(stats::na.omit(covariates))
  numeric_covariates <- unique(stats::na.omit(numeric_covariates))
  numeric_covariates <- intersect(covariates, numeric_covariates)
  controls_df[[reporter]] <- factor(controls_df[[reporter]])
  for (nm in covariates) {
    if (nm %in% numeric_covariates) {
      controls_df[[nm]] <- as.numeric(controls_df[[nm]])
    } else {
      controls_df[[nm]] <- factor(controls_df[[nm]])
    }
  }
  global_covariates <- unique(c(reporter, covariates))
  global <- fit_growth_slope(controls_df, global_covariates)
  if (!is.finite(global$estimate) || !is.finite(global$se) || global$se <= 0) {
    global <- fit_growth_slope(controls_df, character())
  }
  if (!is.finite(global$estimate)) {
    stop("Could not estimate a finite global growth exponent.", call. = FALSE)
  }
  if (!is.finite(global$se) || global$se <= 0) {
    global$se <- 1
  }
  global_alpha <- global$estimate
  global_se <- global$se

  by_reporter <- split(controls_df, controls_df[[reporter]])
  estimates <- lapply(by_reporter, function(d) {
    n <- nrow(d)
    out <- data.frame(
      reporter = as.character(d[[reporter]][1]),
      control_n = n,
      log_growth_sd = stats::sd(d$.log_growth),
      a_raw = NA_real_,
      a_raw_se = NA_real_,
      a_raw_df = NA_real_,
      alpha_raw = NA_real_,
      alpha_raw_se = NA_real_,
      alpha_raw_df = NA_real_,
      alpha_covariates = "",
      stringsAsFactors = FALSE
    )
    if (n >= min_control_n && is.finite(out$log_growth_sd) && out$log_growth_sd > 0) {
      slope <- fit_growth_slope(d, covariates)
      if (!is.finite(slope$estimate) || !is.finite(slope$se) || slope$se <= 0 ||
          !is.finite(slope$df) || slope$df <= 0) {
        slope <- fit_growth_slope(d, character())
      }
      out$a_raw <- slope$intercept
      out$a_raw_se <- slope$intercept_se
      out$a_raw_df <- slope$df
      out$alpha_raw <- slope$estimate
      out$alpha_raw_se <- slope$se
      out$alpha_raw_df <- slope$df
      out$alpha_covariates <- slope$covariates
    }
    out
  })
  estimates <- do.call(rbind, estimates)

  usable <- is.finite(estimates$alpha_raw) & is.finite(estimates$alpha_raw_se) &
    estimates$alpha_raw_se > 0
  if (sum(usable) >= 3) {
    prior_var <- max(stats::var(estimates$alpha_raw[usable], na.rm = TRUE) -
      stats::median(estimates$alpha_raw_se[usable]^2, na.rm = TRUE), 1e-4)
  } else {
    prior_var <- max(global_se^2, 1e-4)
  }

  estimates$alpha_global <- global_alpha
  estimates$alpha_global_se <- global_se
  estimates$alpha_global_covariates <- global$covariates
  estimates$alpha_prior_var <- prior_var
  estimates$alpha_prior_sd <- sqrt(prior_var)
  estimates$alpha_shrunk <- estimates$alpha_raw
  estimates$alpha_shrunk_se <- estimates$alpha_raw_se

  if (isTRUE(shrink)) {
    raw_var <- estimates$alpha_raw_se^2
    estimates$alpha_shrunk[usable] <- (
      estimates$alpha_raw[usable] / raw_var[usable] + global_alpha / prior_var
    ) / (1 / raw_var[usable] + 1 / prior_var)
    estimates$alpha_shrunk_se[usable] <- sqrt(1 / (1 / raw_var[usable] + 1 / prior_var))
  }

  estimates$alpha_shrunk[!usable] <- global_alpha
  estimates$alpha_shrunk_se[!usable] <- global_se
  if (!is.null(alpha_bounds)) {
    estimates$alpha_shrunk <- pmin(pmax(estimates$alpha_shrunk, alpha_bounds[1]), alpha_bounds[2])
  }

  estimates$alpha_fixed_one <- 1
  estimates$alpha_diff_from_one <- estimates$alpha_shrunk - 1
  rownames(estimates) <- NULL
  estimates[order(estimates$reporter), ]
}
