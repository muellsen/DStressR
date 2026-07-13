technical_formula <- function(technical) {
  if (length(technical) == 0) {
    "1"
  } else {
    paste(c("1", technical), collapse = " + ")
  }
}

make_formulas <- function(technical) {
  tech <- technical_formula(technical)
  list(
    total = stats::as.formula(paste(".response ~ .promoter + .compound +", tech)),
    full = stats::as.formula(paste(".response ~ .promoter * .compound +", tech))
  )
}

safe_vcov <- function(fit) {
  vc <- tryCatch(stats::vcov(fit), error = function(e) NULL)
  if (is.null(vc)) {
    matrix(NA_real_, nrow = length(stats::coef(fit)), ncol = length(stats::coef(fit)))
  } else {
    vc[is.na(vc)] <- 0
    vc
  }
}

contrast_estimate <- function(fit, newdata_a, newdata_b) {
  terms_obj <- stats::delete.response(stats::terms(fit))
  x_a <- stats::model.matrix(terms_obj, newdata_a, contrasts.arg = fit$contrasts)
  x_b <- stats::model.matrix(terms_obj, newdata_b, contrasts.arg = fit$contrasts)
  contrast <- drop(x_a - x_b)
  beta <- stats::coef(fit)
  beta[is.na(beta)] <- 0
  estimate <- sum(contrast * beta)
  vc <- safe_vcov(fit)
  se <- sqrt(drop(t(contrast) %*% vc %*% contrast))
  df <- stats::df.residual(fit)
  statistic <- estimate / se
  pvalue <- 2 * stats::pt(abs(statistic), df = df, lower.tail = FALSE)
  c(estimate = estimate, std_error = se, statistic = statistic, pvalue = pvalue)
}

representative_rows <- function(assay, technical) {
  reps <- assay[!duplicated(assay$.promoter), c(".promoter", ".compound", technical), drop = FALSE]
  template <- assay[1, c(".promoter", ".compound", technical), drop = FALSE]
  template <- template[rep(1, length(levels(assay$.promoter))), , drop = FALSE]
  template$.promoter <- factor(levels(assay$.promoter), levels = levels(assay$.promoter))
  for (col in technical) {
    level <- names(sort(table(assay[[col]]), decreasing = TRUE))[1]
    template[[col]] <- factor(level, levels = levels(assay[[col]]))
  }
  template
}

eb_shrink <- function(se, residual_sd, df) {
  var_raw <- se^2
  prior_var <- stats::median(var_raw[is.finite(var_raw) & var_raw > 0], na.rm = TRUE)
  if (!is.finite(prior_var)) {
    return(se)
  }
  weight <- df / (df + 4)
  sqrt(weight * var_raw + (1 - weight) * prior_var)
}

#' Fit the DStressR model
#'
#' Fits two linear models. The additive model estimates DMSO-relative total
#' compound effects after promoter and technical adjustment. The interaction
#' model estimates promoter-specific deviations from the compound-wide effect.
#'
#' @param assay A `destress_assay` produced by [prepare_assay()].
#' @param technical Character vector of batch, plate, replicate, or other
#'   technical-factor columns to include.
#' @param empirical_bayes If `TRUE`, lightly shrinks standard errors toward a
#'   common prior variance.
#' @return A `destress_fit` object.
#' @export
fit_destress <- function(assay, technical = NULL, empirical_bayes = TRUE) {
  if (!inherits(assay, "destress_assay")) {
    stop("`assay` must be produced by prepare_assay().", call. = FALSE)
  }
  technical <- technical[!is.na(technical) & nzchar(technical)]
  missing_technical <- setdiff(technical, names(assay))
  if (length(missing_technical) > 0) {
    stop("Unknown technical columns: ", paste(missing_technical, collapse = ", "), call. = FALSE)
  }
  formulas <- make_formulas(technical)
  total_fit <- stats::lm(formulas$total, data = assay, na.action = stats::na.exclude)
  full_fit <- stats::lm(formulas$full, data = assay, na.action = stats::na.exclude)

  structure(
    list(
      total_fit = total_fit,
      full_fit = full_fit,
      assay_info = attr(assay, "destress"),
      levels = list(
        promoter = levels(assay$.promoter),
        compound = levels(assay$.compound)
      ),
      technical = technical,
      empirical_bayes = empirical_bayes
    ),
    class = "destress_fit"
  )
}

#' Extract model results
#'
#' @param fit A `destress_fit` object.
#' @param compounds Optional compound subset.
#' @param promoters Optional promoter subset.
#' @return A data frame with total and promoter-specific effects.
#' @export
results <- function(fit, compounds = NULL, promoters = NULL) {
  if (!inherits(fit, "destress_fit")) {
    stop("`fit` must be a destress_fit.", call. = FALSE)
  }
  control <- fit$assay_info$control
  all_promoters <- fit$levels$promoter
  all_compounds <- setdiff(fit$levels$compound, control)
  promoters <- if (is.null(promoters)) all_promoters else intersect(promoters, all_promoters)
  compounds <- if (is.null(compounds)) all_compounds else intersect(compounds, all_compounds)

  grid <- expand.grid(
    promoter = promoters,
    compound = compounds,
    stringsAsFactors = FALSE
  )
  if (nrow(grid) == 0) {
    return(data.frame())
  }

  base <- data.frame(
    .promoter = factor(grid$promoter, levels = fit$levels$promoter),
    .compound = factor(control, levels = fit$levels$compound)
  )
  comp <- base
  comp$.compound <- factor(grid$compound, levels = fit$levels$compound)

  for (col in fit$technical) {
    level <- names(sort(table(stats::model.frame(fit$full_fit)[[col]]), decreasing = TRUE))[1]
    base[[col]] <- factor(level, levels = levels(stats::model.frame(fit$full_fit)[[col]]))
    comp[[col]] <- base[[col]]
  }

  total <- t(vapply(seq_len(nrow(grid)), function(i) {
    contrast_estimate(fit$total_fit, comp[i, , drop = FALSE], base[i, , drop = FALSE])
  }, numeric(4)))

  full_total <- t(vapply(seq_len(nrow(grid)), function(i) {
    contrast_estimate(fit$full_fit, comp[i, , drop = FALSE], base[i, , drop = FALSE])
  }, numeric(4)))

  global <- t(vapply(compounds, function(cmp) {
    rows <- comp[grid$compound == cmp, , drop = FALSE]
    refs <- base[grid$compound == cmp, , drop = FALSE]
    vals <- t(vapply(seq_len(nrow(rows)), function(i) {
      contrast_estimate(fit$full_fit, rows[i, , drop = FALSE], refs[i, , drop = FALSE])
    }, numeric(4)))
    colMeans(vals, na.rm = TRUE)
  }, numeric(4)))
  global_df <- data.frame(compound = compounds, global_effect = global[, "estimate"])

  out <- data.frame(
    promoter = grid$promoter,
    compound = grid$compound,
    total_effect = full_total[, "estimate"],
    total_se = full_total[, "std_error"],
    total_statistic = full_total[, "statistic"],
    total_pvalue = full_total[, "pvalue"],
    additive_total_effect = total[, "estimate"],
    additive_total_se = total[, "std_error"],
    stringsAsFactors = FALSE
  )
  out <- merge(out, global_df, by = "compound", sort = FALSE)
  out$specific_effect <- out$total_effect - out$global_effect
  out$specific_se <- out$total_se
  if (isTRUE(fit$empirical_bayes)) {
    out$specific_se <- eb_shrink(out$specific_se, stats::sigma(fit$full_fit), stats::df.residual(fit$full_fit))
  }
  out$specific_statistic <- out$specific_effect / out$specific_se
  out$specific_pvalue <- 2 * stats::pt(abs(out$specific_statistic),
                                       df = stats::df.residual(fit$full_fit),
                                       lower.tail = FALSE)
  out$total_padj <- stats::p.adjust(out$total_pvalue, method = "BH")
  out$specific_padj <- stats::p.adjust(out$specific_pvalue, method = "BH")
  out[order(out$promoter, out$compound), ]
}

#' Adjust p-values within promoter
#'
#' @param table Result table from [results()].
#' @param pvalue P-value column.
#' @param output Name of adjusted p-value column.
#' @param method Passed to [stats::p.adjust()].
#' @export
adjust_pvalues <- function(table, pvalue = "specific_pvalue", output = "specific_padj_by_promoter",
                           method = "BH") {
  split_idx <- split(seq_len(nrow(table)), table$promoter)
  table[[output]] <- NA_real_
  for (idx in split_idx) {
    table[[output]][idx] <- stats::p.adjust(table[[pvalue]][idx], method = method)
  }
  table
}

#' Call differential stress hits
#'
#' @param table Result table from [results()].
#' @param fdr FDR threshold.
#' @param lfc Minimum absolute effect size.
#' @param effect Effect column, usually `specific_effect` or `total_effect`.
#' @param padj Adjusted p-value column.
#' @export
call_hits <- function(table, fdr = 0.05, lfc = 0, effect = "specific_effect",
                      padj = "specific_padj") {
  hit <- rep("Not DE", nrow(table))
  sig <- is.finite(table[[padj]]) & table[[padj]] < fdr & abs(table[[effect]]) >= lfc
  hit[sig & table[[effect]] > 0] <- "Upregulated"
  hit[sig & table[[effect]] < 0] <- "Downregulated"
  table$hit <- hit
  table
}

#' Summarize model dimensions
#'
#' @param fit A `destress_fit`.
#' @export
model_matrix_report <- function(fit) {
  data.frame(
    model = c("additive", "interaction"),
    n_observations = c(stats::nobs(fit$total_fit), stats::nobs(fit$full_fit)),
    n_coefficients = c(length(stats::coef(fit$total_fit)), length(stats::coef(fit$full_fit))),
    residual_df = c(stats::df.residual(fit$total_fit), stats::df.residual(fit$full_fit)),
    sigma = c(stats::sigma(fit$total_fit), stats::sigma(fit$full_fit))
  )
}
