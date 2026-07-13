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

#' List available DStressR presets
#'
#' @return A character vector of preset names accepted by [fit_destress()].
#' @export
destress_presets <- function() {
  c("model", "median_polish_legacy", "empty_vector_control")
}

normalize_destress_preset <- function(preset) {
  if (is.null(preset)) {
    return(NULL)
  }
  if (length(preset) != 1 || is.na(preset) || !nzchar(preset)) {
    stop("`preset` must be one preset name.", call. = FALSE)
  }
  preset <- gsub("-", "_", tolower(preset), fixed = TRUE)
  aliases <- c(
    destress = "model",
    model_based = "model",
    median_polish = "median_polish_legacy",
    medianpolish = "median_polish_legacy",
    legacy = "median_polish_legacy",
    empty_vector = "empty_vector_control",
    evc = "empty_vector_control"
  )
  if (preset %in% names(aliases)) {
    preset <- aliases[[preset]]
  }
  choices <- destress_presets()
  if (!preset %in% choices) {
    stop(
      "Unknown preset `", preset, "`. Available presets are: ",
      paste(choices, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  preset
}

normalize_stage_choice <- function(value, choices, aliases, name) {
  if (length(value) != 1 || is.na(value) || !nzchar(value)) {
    stop("`", name, "` must be one choice.", call. = FALSE)
  }
  value <- gsub("-", "_", tolower(value), fixed = TRUE)
  if (value %in% names(aliases)) {
    value <- aliases[[value]]
  }
  if (!value %in% choices) {
    stop(
      "Unknown `", name, "` choice `", value, "`. Available choices are: ",
      paste(choices, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  value
}

destress_preset_stages <- function(preset, empirical_bayes = TRUE) {
  preset <- normalize_destress_preset(preset)
  switch(
    preset,
    model = list(
      normalization = "linear_model",
      testing = if (isTRUE(empirical_bayes)) "moderated_t" else "student_t",
      aggregation = "none",
      adjustment = "global"
    ),
    median_polish_legacy = list(
      normalization = "median_polish",
      testing = "gaussian_z",
      aggregation = "max_p",
      adjustment = "by_promoter"
    ),
    empty_vector_control = list(
      normalization = "empty_vector",
      testing = "gaussian_z",
      aggregation = "max_p",
      adjustment = "by_promoter"
    )
  )
}

resolve_destress_stages <- function(preset,
                                    normalization,
                                    testing,
                                    aggregation,
                                    adjustment,
                                    empirical_bayes) {
  stages <- if (is.null(preset)) {
    list(
      normalization = "linear_model",
      testing = if (isTRUE(empirical_bayes)) "moderated_t" else "student_t",
      aggregation = "none",
      adjustment = "global"
    )
  } else {
    destress_preset_stages(preset, empirical_bayes = empirical_bayes)
  }

  if (!is.null(normalization)) {
    stages$normalization <- normalize_stage_choice(
      normalization,
      choices = c("linear_model", "median_polish", "empty_vector"),
      aliases = c(model = "linear_model", lm = "linear_model", evc = "empty_vector"),
      name = "normalization"
    )
  }
  if (!is.null(testing)) {
    stages$testing <- normalize_stage_choice(
      testing,
      choices = c("student_t", "moderated_t", "gaussian_z"),
      aliases = c(t = "student_t", model_t = "student_t", eb = "moderated_t", empirical_bayes = "moderated_t", z = "gaussian_z"),
      name = "testing"
    )
  }
  if (!is.null(aggregation)) {
    stages$aggregation <- normalize_stage_choice(
      aggregation,
      choices = c("none", "max_p"),
      aliases = c(max = "max_p", conservative = "max_p"),
      name = "aggregation"
    )
  }
  if (!is.null(adjustment)) {
    stages$adjustment <- normalize_stage_choice(
      adjustment,
      choices = c("global", "by_promoter", "none"),
      aliases = c(promoter = "by_promoter", within_promoter = "by_promoter"),
      name = "adjustment"
    )
  }

  if (stages$normalization == "linear_model") {
    if (!stages$testing %in% c("student_t", "moderated_t")) {
      stop("`normalization = \"linear_model\"` currently supports `testing = \"student_t\"` or `\"moderated_t\"`.", call. = FALSE)
    }
    if (stages$aggregation != "none") {
      stop("`normalization = \"linear_model\"` currently supports `aggregation = \"none\"`.", call. = FALSE)
    }
  } else {
    expected <- if (stages$normalization == "median_polish") "median-polish" else "empty-vector"
    if (stages$testing != "gaussian_z" || stages$aggregation != "max_p" || stages$adjustment != "by_promoter") {
      stop(
        "The ", expected, " compatibility path currently requires ",
        "`testing = \"gaussian_z\"`, `aggregation = \"max_p\"`, and ",
        "`adjustment = \"by_promoter\"`.",
        call. = FALSE
      )
    }
  }
  stages
}

#' Fit DStressR with staged statistical options
#'
#' `fit_destress()` is the main DStressR entry point. By default it fits the
#' model-based analysis, but it can also run named compatibility presets for the
#' legacy median-polish and Empty Vector Control analyses.
#'
#' The staged options make the major statistical choices explicit:
#' normalization, test statistic and p-value calculation, replicate aggregation,
#' and p-value adjustment. Only implemented combinations are accepted. For the
#' model-based path, growth-response normalization is performed upstream by
#' [prepare_assay()], where `growth_exponent` can be fixed, estimated, or
#' supplied as promoter-specific values.
#'
#' @param assay A `destress_assay` produced by [prepare_assay()] or a raw assay
#'   data frame for `normalization = "linear_model"`, or a long expression
#'   table for the compatibility presets.
#' @param technical Character vector of batch, plate, replicate, or other
#'   technical-factor columns to include.
#' @param empirical_bayes If `TRUE`, lightly shrinks standard errors toward a
#'   common prior variance. This maps to `testing = "moderated_t"` for the
#'   model path; `FALSE` maps to `testing = "student_t"`.
#' @param normalization One of `"linear_model"`, `"median_polish"`, or
#'   `"empty_vector"`. `"model"` and `"evc"` are accepted aliases.
#' @param testing One of `"student_t"`, `"moderated_t"`, or `"gaussian_z"`.
#' @param aggregation One of `"none"` or `"max_p"`.
#' @param adjustment One of `"global"`, `"by_promoter"`, or `"none"`.
#' @param preset Optional named preset: `"model"`, `"median_polish_legacy"`, or
#'   `"empty_vector_control"`. Common aliases such as `"median_polish"` and
#'   `"evc"` are accepted.
#' @param ... For `normalization = "linear_model"` with a raw data frame,
#'   arguments passed to [prepare_assay()], including `growth_exponent`. For
#'   compatibility presets, arguments passed to the selected engine.
#' @return A fitted DStressR object. The model path returns a `destress_fit`;
#'   compatibility presets return their corresponding legacy result objects.
#' @export
fit_destress <- function(assay,
                         technical = NULL,
                         empirical_bayes = TRUE,
                         normalization = NULL,
                         testing = NULL,
                         aggregation = NULL,
                         adjustment = NULL,
                         preset = NULL,
                         ...) {
  preset <- normalize_destress_preset(preset)
  stages <- resolve_destress_stages(
    preset = preset,
    normalization = normalization,
    testing = testing,
    aggregation = aggregation,
    adjustment = adjustment,
    empirical_bayes = empirical_bayes
  )

  if (stages$normalization == "median_polish") {
    fit <- fit_median_polish(assay, ...)
    attr(fit, "destress_preset") <- if (is.null(preset)) "median_polish_legacy" else preset
    attr(fit, "destress_stages") <- stages
    return(fit)
  }
  if (stages$normalization == "empty_vector") {
    fit <- fit_empty_vector_control(assay, ...)
    attr(fit, "destress_preset") <- if (is.null(preset)) "empty_vector_control" else preset
    attr(fit, "destress_stages") <- stages
    return(fit)
  }

  empirical_bayes <- identical(stages$testing, "moderated_t")
  if (!inherits(assay, "destress_assay")) {
    if (!is.data.frame(assay)) {
      stop("`assay` must be a data frame or be produced by prepare_assay().", call. = FALSE)
    }
    assay <- prepare_assay(assay, ...)
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
      empirical_bayes = empirical_bayes,
      stages = stages,
      preset = if (is.null(preset)) "model" else preset,
      adjustment = stages$adjustment
    ),
    class = "destress_fit"
  )
}

adjust_destress_pvalues <- function(pvalue, groups = NULL, adjustment = "global") {
  out <- rep(NA_real_, length(pvalue))
  finite <- is.finite(pvalue)
  if (adjustment == "none") {
    out[finite] <- pvalue[finite]
  } else if (adjustment == "global") {
    out[finite] <- stats::p.adjust(pvalue[finite], method = "BH")
  } else if (adjustment == "by_promoter") {
    split_idx <- split(seq_along(pvalue), groups)
    for (idx in split_idx) {
      finite_idx <- idx[is.finite(pvalue[idx])]
      out[finite_idx] <- stats::p.adjust(pvalue[finite_idx], method = "BH")
    }
  }
  out
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
  adjustment <- if (is.null(fit$adjustment)) "global" else fit$adjustment
  out$total_padj <- adjust_destress_pvalues(out$total_pvalue, out$promoter, adjustment)
  out$specific_padj <- adjust_destress_pvalues(out$specific_pvalue, out$promoter, adjustment)
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
