technical_formula <- function(technical) {
  if (length(technical) == 0) {
    "1"
  } else {
    paste(c("1", technical), collapse = " + ")
  }
}

destress_progress <- function(verbose, ...) {
  if (isTRUE(verbose)) {
    message(...)
  }
}

make_formulas <- function(technical) {
  tech <- technical_formula(technical)
  list(
    total = stats::as.formula(paste(".response ~ .reporter + .perturbation +", tech)),
    full = stats::as.formula(paste(".response ~ .reporter * .perturbation +", tech)),
    technical = stats::as.formula(paste(".response ~", tech))
  )
}

safe_vcov <- function(fit) {
  if (inherits(fit, "destress_sparse_lm")) {
    return(fit$vcov)
  }
  vc <- tryCatch(stats::vcov(fit), error = function(e) NULL)
  if (is.null(vc)) {
    matrix(NA_real_, nrow = length(stats::coef(fit)), ncol = length(stats::coef(fit)))
  } else {
    vc[is.na(vc)] <- 0
    vc
  }
}

fit_coef <- function(fit) {
  if (inherits(fit, "destress_sparse_lm")) {
    return(fit$coefficients)
  }
  stats::coef(fit)
}

fit_df_residual <- function(fit) {
  if (inherits(fit, "destress_sparse_lm")) {
    return(fit$df.residual)
  }
  stats::df.residual(fit)
}

fit_sigma <- function(fit) {
  if (inherits(fit, "destress_sparse_lm")) {
    return(fit$sigma)
  }
  stats::sigma(fit)
}

fit_nobs <- function(fit) {
  if (inherits(fit, "destress_sparse_lm")) {
    return(fit$nobs)
  }
  stats::nobs(fit)
}

fit_terms <- function(fit) {
  if (inherits(fit, "destress_sparse_lm")) {
    return(fit$terms)
  }
  stats::terms(fit)
}

fit_model_frame <- function(fit) {
  if (inherits(fit, "destress_sparse_lm")) {
    return(fit$model)
  }
  stats::model.frame(fit)
}

fit_contrasts <- function(fit) {
  if (inherits(fit, "destress_sparse_lm")) {
    return(fit$contrasts)
  }
  fit$contrasts
}

contrast_estimate <- function(fit, newdata_a, newdata_b) {
  terms_obj <- stats::delete.response(fit_terms(fit))
  x_a <- stats::model.matrix(terms_obj, newdata_a, contrasts.arg = fit_contrasts(fit))
  x_b <- stats::model.matrix(terms_obj, newdata_b, contrasts.arg = fit_contrasts(fit))
  contrast <- drop(x_a - x_b)
  beta <- fit_coef(fit)
  keep <- intersect(names(beta), names(contrast))
  contrast <- contrast[keep]
  beta <- beta[keep]
  beta[is.na(beta)] <- 0
  estimate <- sum(contrast * beta)
  vc <- safe_vcov(fit)
  vc <- vc[keep, keep, drop = FALSE]
  se <- sqrt(drop(t(contrast) %*% vc %*% contrast))
  df <- fit_df_residual(fit)
  statistic <- estimate / se
  pvalue <- 2 * stats::pt(abs(statistic), df = df, lower.tail = FALSE)
  c(estimate = estimate, std_error = se, statistic = statistic, pvalue = pvalue)
}

contrast_estimates <- function(fit, newdata_a, newdata_b, center_by = NULL) {
  terms_obj <- stats::delete.response(fit_terms(fit))
  x_a <- stats::model.matrix(terms_obj, newdata_a, contrasts.arg = fit_contrasts(fit))
  x_b <- stats::model.matrix(terms_obj, newdata_b, contrasts.arg = fit_contrasts(fit))
  contrast <- x_a - x_b
  if (!is.null(center_by)) {
    stopifnot(length(center_by) == nrow(contrast))
    for (idx in split(seq_len(nrow(contrast)), center_by)) {
      contrast[idx, ] <- sweep(contrast[idx, , drop = FALSE], 2,
                              colMeans(contrast[idx, , drop = FALSE]), "-")
    }
  }
  beta <- fit_coef(fit)
  keep <- intersect(names(beta), colnames(contrast))
  contrast <- contrast[, keep, drop = FALSE]
  beta <- beta[keep]
  beta[is.na(beta)] <- 0
  estimate <- as.numeric(contrast %*% beta)
  vc <- safe_vcov(fit)
  vc <- vc[keep, keep, drop = FALSE]
  se <- sqrt(rowSums((contrast %*% vc) * contrast))
  df <- fit_df_residual(fit)
  statistic <- estimate / se
  pvalue <- 2 * stats::pt(abs(statistic), df = df, lower.tail = FALSE)
  cbind(estimate = estimate, std_error = se, statistic = statistic, pvalue = pvalue)
}

representative_rows <- function(assay, technical) {
  reps <- assay[!duplicated(assay$.reporter), c(".reporter", ".perturbation", technical), drop = FALSE]
  template <- assay[1, c(".reporter", ".perturbation", technical), drop = FALSE]
  template <- template[rep(1, length(levels(assay$.reporter))), , drop = FALSE]
  template$.reporter <- factor(levels(assay$.reporter), levels = levels(assay$.reporter))
  for (col in technical) {
    level <- names(sort(table(assay[[col]]), decreasing = TRUE))[1]
    template[[col]] <- factor(level, levels = levels(assay[[col]]))
  }
  template
}

trigamma_inverse <- function(y) {
  if (!is.finite(y) || y <= 0) {
    return(Inf)
  }

  lo <- .Machine$double.eps
  hi <- 1
  while (psigamma(hi, deriv = 1) > y && hi < 1e12) {
    hi <- hi * 2
  }
  if (hi >= 1e12) {
    return(Inf)
  }

  for (i in seq_len(80)) {
    mid <- (lo + hi) / 2
    if (psigamma(mid, deriv = 1) > y) {
      lo <- mid
    } else {
      hi <- mid
    }
  }
  hi
}

estimate_eb_prior_variance <- function(var_raw, df, max_prior_df = 1e6) {
  if (length(df) == 1L) df <- rep(df, length(var_raw))
  if (length(df) != length(var_raw)) {
    stop("`df` must be scalar or have the same length as the variances.", call. = FALSE)
  }
  ok <- is.finite(var_raw) & var_raw > 0 & is.finite(df) & df > 0
  if (sum(ok) < 2) {
    return(list(var = NA_real_, df = 0))
  }

  df <- df[ok]
  # Remove the known log-chi-square mean separately for each contrast.
  log_var <- log(var_raw[ok]) - digamma(df / 2) + log(df / 2)
  log_var_var <- stats::var(log_var)
  if (!is.finite(log_var_var)) {
    return(list(var = stats::median(var_raw[ok], na.rm = TRUE), df = 0))
  }

  sampling_var <- mean(psigamma(df / 2, deriv = 1))
  excess_var <- log_var_var - sampling_var
  prior_df <- if (is.finite(excess_var) && excess_var > 0) {
    2 * trigamma_inverse(excess_var)
  } else {
    max_prior_df
  }
  prior_df <- min(prior_df, max_prior_df)

  if (!is.finite(prior_df) || prior_df <= 0) {
    return(list(var = stats::median(var_raw[ok], na.rm = TRUE), df = 0))
  }

  log_prior_var <- mean(log_var) +
    digamma(prior_df / 2) - log(prior_df / 2)
  prior_var <- exp(log_prior_var)
  if (!is.finite(prior_var) || prior_var <= 0) {
    prior_var <- stats::median(var_raw[ok], na.rm = TRUE)
  }
  list(var = prior_var, df = prior_df)
}

eb_moderate_se <- function(se, df) {
  var_raw <- se^2
  df <- as.numeric(df)
  prior <- estimate_eb_prior_variance(var_raw, df)
  if (!is.finite(prior$var) || !is.finite(prior$df) || prior$df <= 0) {
    return(list(se = se, df = df, prior_var = prior$var, prior_df = prior$df))
  }

  moderated_var <- (df * var_raw + prior$df * prior$var) / (df + prior$df)
  moderated_se <- sqrt(moderated_var)
  moderated_se[!is.finite(se)] <- se[!is.finite(se)]
  list(
    se = moderated_se,
    df = df + prior$df,
    prior_var = prior$var,
    prior_df = prior$df
  )
}

centered_contrast_df <- function(variance, df) {
  n <- length(variance)
  if (length(df) == 1L) df <- rep(df, n)
  if (length(df) != n) {
    stop("`df` must be scalar or have the same length as the variances.", call. = FALSE)
  }
  if (n <= 1L || any(!is.finite(variance) | variance < 0) ||
      any(is.na(df) | df <= 0) || max(variance) <= 0) {
    return(rep(NA_real_, n))
  }
  # Rescaling avoids overflow; it cancels from the moment-matched df.
  v <- variance / max(variance)
  own_weight <- ((n - 1) / n)^2
  other_weight <- 1 / n^2
  combined <- own_weight * v + other_weight * (sum(v) - v)
  component <- v^2 / df
  denominator <- own_weight^2 * component +
    other_weight^2 * (sum(component) - component)
  out <- combined^2 / denominator
  out[combined <= 0] <- NA_real_
  out
}

wald_t_test <- function(estimate, se, df) {
  statistic <- estimate / se
  zero_se <- is.finite(se) & se == 0
  statistic[zero_se & abs(estimate) < sqrt(.Machine$double.eps)] <- 0
  statistic[zero_se & abs(estimate) >= sqrt(.Machine$double.eps)] <-
    sign(estimate[zero_se & abs(estimate) >= sqrt(.Machine$double.eps)]) * Inf
  pvalue <- 2 * stats::pt(abs(statistic), df = df, lower.tail = FALSE)
  list(statistic = statistic, pvalue = pvalue)
}

validate_background_rank <- function(background_rank) {
  if (is.null(background_rank)) {
    return(0L)
  }
  if (length(background_rank) != 1 || is.na(background_rank)) {
    stop("`background_rank` must be one non-negative integer.", call. = FALSE)
  }
  if (!is.numeric(background_rank) || background_rank < 0 || background_rank != floor(background_rank)) {
    stop("`background_rank` must be one non-negative integer.", call. = FALSE)
  }
  as.integer(background_rank)
}

validate_background_rank_request <- function(background_rank) {
  if (is.character(background_rank) && length(background_rank) == 1 && identical(tolower(background_rank), "auto")) {
    return(list(auto = TRUE, rank = NA_integer_))
  }
  list(auto = FALSE, rank = validate_background_rank(background_rank))
}

validate_background_rank_threshold <- function(threshold) {
  threshold <- as.numeric(threshold)
  if (length(threshold) != 1 || !is.finite(threshold) || threshold <= 0 || threshold >= 1) {
    stop("`background_rank_threshold` must be one finite number between 0 and 1.", call. = FALSE)
  }
  threshold
}

validate_background_impute <- function(impute) {
  if (is.null(impute)) {
    return("column_mean")
  }
  match.arg(impute, c("column_mean", "global_mean", "zero"))
}

impute_effect_matrix <- function(mat, method = "column_mean") {
  method <- validate_background_impute(method)
  observed <- is.finite(mat)
  out <- mat
  if (method == "zero") {
    out[!observed] <- 0
    return(out)
  }
  finite <- mat[observed]
  global_mean <- if (length(finite) > 0) mean(finite) else 0
  if (!is.finite(global_mean)) {
    global_mean <- 0
  }
  if (method == "global_mean") {
    out[!observed] <- global_mean
    return(out)
  }
  for (j in seq_len(ncol(out))) {
    miss_j <- !observed[, j]
    if (!any(miss_j)) {
      next
    }
    mean_j <- mean(out[observed[, j], j], na.rm = TRUE)
    if (!is.finite(mean_j)) {
      mean_j <- global_mean
    }
    out[miss_j, j] <- mean_j
  }
  out
}

select_background_rank <- function(diagnostics) {
  if (is.null(diagnostics) || nrow(diagnostics) == 0) {
    return(0L)
  }
  reference <- if ("null_threshold" %in% names(diagnostics)) {
    diagnostics$null_threshold
  } else if ("null_q99" %in% names(diagnostics)) {
    diagnostics$null_q99
  } else if ("null_q95" %in% names(diagnostics)) {
    diagnostics$null_q95
  } else {
    rep(Inf, nrow(diagnostics))
  }
  pass <- is.finite(diagnostics$observed) & is.finite(reference) & diagnostics$observed > reference
  if (!any(pass)) {
    return(0L)
  }
  fail <- which(!pass)
  if (length(fail) == 0) {
    return(max(diagnostics$component))
  }
  first_fail <- fail[1]
  if (first_fail == 1) {
    0L
  } else {
    max(diagnostics$component[seq_len(first_fail - 1L)])
  }
}

low_rank_background_effect <- function(table, effect, rank, impute = "column_mean") {
  rank <- validate_background_rank(rank)
  impute <- validate_background_impute(impute)
  if (rank == 0 || nrow(table) == 0) {
    return(rep(0, nrow(table)))
  }

  required <- c("reporter", "perturbation", effect)
  missing <- setdiff(required, names(table))
  if (length(missing) > 0) {
    stop("Cannot estimate low-rank background; missing columns: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }

  reporters <- unique(as.character(table$reporter))
  perturbations <- unique(as.character(table$perturbation))
  mat <- matrix(
    NA_real_,
    nrow = length(reporters),
    ncol = length(perturbations),
    dimnames = list(reporters, perturbations)
  )
  idx <- cbind(
    match(as.character(table$reporter), reporters),
    match(as.character(table$perturbation), perturbations)
  )
  mat[idx] <- as.numeric(table[[effect]])

  observed <- is.finite(mat)
  decomp <- impute_effect_matrix(mat, method = impute)

  rank <- min(rank, nrow(decomp), ncol(decomp))
  if (rank == 0 || all(abs(decomp[observed]) < sqrt(.Machine$double.eps))) {
    return(rep(0, nrow(table)))
  }

  sv <- svd(decomp, nu = rank, nv = rank)
  keep <- seq_len(rank)
  low_rank <- sv$u[, keep, drop = FALSE] %*%
    (diag(sv$d[keep], nrow = rank, ncol = rank) %*% t(sv$v[, keep, drop = FALSE]))

  low_rank[!observed] <- NA_real_
  as.numeric(low_rank[idx])
}

effect_matrix_from_table <- function(table, effect, reporter = "reporter", perturbation = "perturbation") {
  required <- c(reporter, perturbation, effect)
  missing <- setdiff(required, names(table))
  if (length(missing) > 0) {
    stop("Cannot build effect matrix; missing columns: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }

  reporters <- unique(as.character(table[[reporter]]))
  perturbations <- unique(as.character(table[[perturbation]]))
  mat <- matrix(
    NA_real_,
    nrow = length(reporters),
    ncol = length(perturbations),
    dimnames = list(reporters, perturbations)
  )
  idx <- cbind(
    match(as.character(table[[reporter]]), reporters),
    match(as.character(table[[perturbation]]), perturbations)
  )
  mat[idx] <- as.numeric(table[[effect]])
  mat
}

#' Diagnose low-rank background structure
#'
#' Computes singular values of a reporter-by-perturbation effect matrix and compares
#' them with a permutation null. The default permutation shuffles reporter
#' labels within each perturbation, preserving the perturbation-wise marginal
#' distribution while breaking shared reporter-loading structure.
#'
#' @param table Data frame with reporter, perturbation, and effect columns.
#' @param effect Numeric effect column to decompose, usually `total_effect` or
#'   `background_adjusted_effect`.
#' @param reporter,perturbation Column names identifying reporters and perturbations.
#' @param rank_max Maximum component index to report.
#' @param permutations Number of null permutations. Use `0` to skip the null.
#' @param threshold Permutation reference quantile used for automatic rank
#'   selection and reported as `null_threshold`.
#' @param impute Method used to fill missing matrix entries before
#'   singular-value decomposition.
#' @param seed Optional random seed for reproducible permutations.
#' @return A data frame with observed singular values, variance fractions, and
#'   optional permutation summaries.
#' @export
background_rank_diagnostics <- function(table,
                                        effect = "total_effect",
                                        reporter = "reporter",
                                        perturbation = "perturbation",
                                        rank_max = 10,
                                        permutations = 100,
                                        threshold = 0.99,
                                        impute = "column_mean",
                                        seed = NULL) {
  rank_max <- validate_background_rank(rank_max)
  permutations <- validate_background_rank(permutations)
  threshold <- validate_background_rank_threshold(threshold)
  impute <- validate_background_impute(impute)
  mat <- effect_matrix_from_table(
    table,
    effect = effect,
    reporter = reporter,
    perturbation = perturbation
  )
  observed <- is.finite(mat)
  decomp <- impute_effect_matrix(mat, method = impute)

  rank_max <- min(rank_max, nrow(decomp), ncol(decomp))
  if (rank_max == 0) {
    return(data.frame())
  }

  singular_values <- svd(decomp, nu = 0, nv = 0)$d
  total_ss <- sum(singular_values^2)
  component <- seq_len(rank_max)
  observed_sv <- singular_values[component]
  prop_var <- if (total_ss > 0) observed_sv^2 / total_ss else rep(NA_real_, rank_max)

  null_median <- null_q95 <- null_q99 <- null_threshold <- rep(NA_real_, rank_max)
  if (permutations > 0) {
    if (!is.null(seed)) {
      old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
      } else {
        NULL
      }
      on.exit({
        if (is.null(old_seed)) {
          rm(".Random.seed", envir = .GlobalEnv)
        } else {
          assign(".Random.seed", old_seed, envir = .GlobalEnv)
        }
      }, add = TRUE)
      set.seed(seed)
    }

    null_sv <- matrix(NA_real_, nrow = permutations, ncol = rank_max)
    for (b in seq_len(permutations)) {
      permuted <- mat
      for (j in seq_len(ncol(permuted))) {
        obs_j <- observed[, j]
        if (sum(obs_j) > 1) {
          permuted[obs_j, j] <- sample(permuted[obs_j, j])
        }
      }
      permuted <- impute_effect_matrix(permuted, method = impute)
      sv_b <- svd(permuted, nu = 0, nv = 0)$d
      null_sv[b, ] <- sv_b[component]
    }
    null_median <- apply(null_sv, 2, stats::median, na.rm = TRUE)
    null_q95 <- apply(null_sv, 2, stats::quantile, probs = 0.95, na.rm = TRUE)
    null_q99 <- apply(null_sv, 2, stats::quantile, probs = 0.99, na.rm = TRUE)
    null_threshold <- apply(null_sv, 2, stats::quantile, probs = threshold, na.rm = TRUE)
  }

  data.frame(
    component = component,
    observed = observed_sv,
    prop_variance = prop_var,
    cumulative_prop_variance = cumsum(prop_var),
    null_median = null_median,
    null_q95 = null_q95,
    null_q99 = null_q99,
    null_threshold = null_threshold,
    threshold = threshold,
    impute = impute,
    n_reporters = nrow(decomp),
    n_perturbations = ncol(decomp),
    permutations = permutations,
    stringsAsFactors = FALSE
  )
}

modal_factor_level <- function(x) {
  names(sort(table(x), decreasing = TRUE))[1]
}

technical_adjusted_response <- function(fit, assay, technical) {
  if (length(technical) == 0) {
    return(assay$.response)
  }

  reference <- assay
  for (col in technical) {
    reference[[col]] <- factor(modal_factor_level(assay[[col]]), levels = levels(assay[[col]]))
  }
  observed_fit <- suppressWarnings(stats::predict(fit, newdata = assay))
  reference_fit <- suppressWarnings(stats::predict(fit, newdata = reference))
  assay$.response - (observed_fit - reference_fit)
}

independent_columns <- function(x) {
  qr_x <- qr(x)
  x[, sort(qr_x$pivot[seq_len(qr_x$rank)]), drop = FALSE]
}

technical_reference_row <- function(d, technical, n) {
  out <- data.frame(row.names = seq_len(n))
  for (col in technical) {
    x <- d[[col]]
    if (is.factor(x)) {
      level <- modal_factor_level(x)
      out[[col]] <- factor(level, levels = levels(x), ordered = is.ordered(x))
    } else if (is.numeric(x) || is.integer(x)) {
      out[[col]] <- stats::median(as.numeric(x), na.rm = TRUE)
    } else {
      level <- names(sort(table(x), decreasing = TRUE))[1]
      out[[col]] <- rep(level, n)
    }
  }
  out
}

active_reporter_technical <- function(d, technical) {
  keep <- vapply(technical, function(col) {
    x <- d[[col]]
    if (is.factor(x)) {
      nlevels(droplevels(x)) > 1
    } else if (is.numeric(x) || is.integer(x)) {
      vals <- unique(as.numeric(x[is.finite(x)]))
      length(vals) > 1
    } else {
      length(unique(x[!is.na(x)])) > 1
    }
  }, logical(1))
  technical[keep]
}

fit_reporter_effects_lm <- function(d, reporter, perturbations, control, technical) {
  technical <- active_reporter_technical(d, technical)
  d <- droplevels(d)
  form <- stats::as.formula(paste(".response ~ .perturbation +", technical_formula(technical)))
  reporter_fit <- stats::lm(form, data = d, na.action = stats::na.exclude)
  perturbation_levels <- levels(d$.perturbation)

  base <- data.frame(.perturbation = factor(control, levels = perturbation_levels))
  base <- base[rep(1, length(perturbations)), , drop = FALSE]
  comp <- base
  comp$.perturbation <- factor(perturbations, levels = perturbation_levels)
  tech_ref <- technical_reference_row(d, technical, length(perturbations))
  base <- cbind(base, tech_ref)
  comp <- cbind(comp, tech_ref)

  total <- contrast_estimates(reporter_fit, comp, base)
  data.frame(
    reporter = reporter,
    perturbation = perturbations,
    total_effect = total[, "estimate"],
    total_se = total[, "std_error"],
    total_statistic = total[, "statistic"],
    total_pvalue = total[, "pvalue"],
    residual_df = fit_df_residual(reporter_fit),
    sigma = fit_sigma(reporter_fit),
    n_observations = fit_nobs(reporter_fit),
    n_coefficients = length(stats::coef(reporter_fit)),
    stringsAsFactors = FALSE
  )
}

fit_reporter_effects <- function(assay, technical, control, verbose = FALSE, progress_every = 100L) {
  technical_terms <- technical_formula(technical)
  technical_design_formula <- stats::as.formula(paste("~", technical_terms))
  by_reporter <- split(assay, assay$.reporter)

  reporter_names <- names(by_reporter)
  rows <- lapply(seq_along(reporter_names), function(k) {
    reporter <- reporter_names[k]
    if (isTRUE(verbose) &&
        (k == 1L || k == length(reporter_names) || k %% progress_every == 0L)) {
      destress_progress(
        TRUE,
        "DStressR: reporter-specific fits ",
        k,
        "/",
        length(reporter_names),
        "."
      )
    }
    d <- by_reporter[[reporter]]
    d <- d[is.finite(d$.response), , drop = FALSE]
    if (nrow(d) == 0) {
      return(NULL)
    }
    active_technical <- active_reporter_technical(d, technical)
    d <- droplevels(d)

    perturbation <- as.character(d$.perturbation)
    perturbations <- sort(setdiff(unique(perturbation), control))
    if (length(perturbations) == 0 || !any(perturbation == control)) {
      return(NULL)
    }

    closed_form <- tryCatch({
      y <- d$.response
      technical_terms <- technical_formula(active_technical)
      technical_design_formula <- stats::as.formula(paste("~", technical_terms))
      z <- stats::model.matrix(technical_design_formula, data = d)
      z <- independent_columns(z)
      ztz_inv <- solve(crossprod(z))
      mz_y <- as.numeric(stats::lm.fit(z, y)$residuals)
      y_mz_y <- sum(mz_y^2)

      counts_all <- table(factor(perturbation, levels = c(control, perturbations)))
      n_perturbation <- as.numeric(counts_all[perturbations])
      names(n_perturbation) <- perturbations

      z_sum_all <- rowsum(z, factor(perturbation, levels = c(control, perturbations)), reorder = FALSE)
      b <- z_sum_all[perturbations, , drop = FALSE]
      u_all <- rowsum(mz_y, factor(perturbation, levels = c(control, perturbations)), reorder = FALSE)
      u <- as.numeric(u_all[perturbations, , drop = TRUE])
      names(u) <- perturbations

      n_inv <- 1 / n_perturbation
      s <- solve(solve(ztz_inv) - crossprod(b, b * n_inv))
      w <- n_inv * u
      beta <- w + n_inv * as.numeric(b %*% (s %*% crossprod(b, w)))
      a_inv_diag <- n_inv + n_inv^2 * rowSums((b %*% s) * b)

      sse <- y_mz_y - sum(beta * u)
      if (is.finite(sse) && sse < 0 && abs(sse) < sqrt(.Machine$double.eps)) {
        sse <- 0
      }
      df <- length(y) - ncol(z) - length(perturbations)
      sigma2 <- sse / df
      se <- sqrt(sigma2 * a_inv_diag)
      test <- wald_t_test(beta, se, df)

      data.frame(
        reporter = reporter,
        perturbation = perturbations,
        total_effect = beta,
        total_se = se,
        total_statistic = test$statistic,
        total_pvalue = test$pvalue,
        residual_df = df,
        sigma = sqrt(sigma2),
        n_observations = length(y),
        n_coefficients = ncol(z) + length(perturbations),
        stringsAsFactors = FALSE
      )
    }, error = function(e) NULL)

    if (!is.null(closed_form) &&
        nrow(closed_form) > 0 &&
        all(is.finite(closed_form$total_effect)) &&
        all(is.finite(closed_form$total_se))) {
      return(closed_form)
    }

    if (isTRUE(verbose)) {
      destress_progress(
        TRUE,
        "DStressR: reporter-specific QR fallback for reporter ",
        reporter,
        "."
      )
    }
    fit_reporter_effects_lm(d, reporter, perturbations, control, active_technical)
  })

  out <- do.call(rbind, rows)
  if (is.null(out)) {
    out <- data.frame()
  }
  out
}

reporter_effect_results <- function(fit, perturbations = NULL, reporters = NULL) {
  out <- fit$reporter_effects
  if (!is.null(perturbations)) {
    out <- out[out$perturbation %in% perturbations, , drop = FALSE]
  }
  if (nrow(out) == 0) {
    return(data.frame())
  }

  valid_inference <- is.finite(out$total_effect) &
    is.finite(out$total_se) &
    is.finite(out$residual_df) &
    out$residual_df > 0
  out <- out[valid_inference, , drop = FALSE]
  if (nrow(out) == 0) {
    return(data.frame())
  }

  out$additive_total_effect <- NA_real_
  out$additive_total_se <- NA_real_
  out$total_var_raw <- out$total_se^2
  df <- stats::median(out$residual_df[out$residual_df > 0], na.rm = TRUE)
  total_df <- out$residual_df
  if (isTRUE(fit$empirical_bayes)) {
    moderated <- eb_moderate_se(out$total_se, out$residual_df)
    out$total_se <- moderated$se
    total_test <- wald_t_test(out$total_effect, out$total_se, moderated$df)
    out$total_statistic <- total_test$statistic
    out$total_pvalue <- total_test$pvalue
    total_df <- moderated$df
  }

  out$total_var <- out$total_se^2
  out$total_df_for_test <- total_df
  evc <- fit$empty_vector_reporter
  if (!is.null(evc) && nzchar(evc)) {
    evc_rows <- out[out$reporter == evc, c(
      "perturbation", "total_effect", "total_var", "total_var_raw", "residual_df"
    ), drop = FALSE]
    if (nrow(evc_rows) == 0) {
      stop("No fitted effects found for `empty_vector_reporter = \"", evc, "\"`.", call. = FALSE)
    }
    names(evc_rows) <- c("perturbation", "empty_vector_effect", "empty_vector_var",
                         "empty_vector_var_raw", "empty_vector_df")
    out <- merge(out, evc_rows, by = "perturbation", all.x = TRUE, sort = FALSE)
    out$background_adjusted_effect <- out$total_effect - out$empty_vector_effect
    out$background_adjusted_var <- out$total_var + out$empty_vector_var
    out$background_adjusted_var_raw <- out$total_var_raw + out$empty_vector_var_raw
    out$background_adjusted_df_raw <- vapply(seq_len(nrow(out)), function(i) {
      centered_contrast_df(
        c(out$total_var_raw[i], out$empty_vector_var_raw[i]),
        c(out$residual_df[i], out$empty_vector_df[i])
      )[1]
    }, numeric(1))
    out <- out[out$reporter != evc, , drop = FALSE]
  } else {
    out$empty_vector_effect <- NA_real_
    out$empty_vector_var <- NA_real_
    out$background_adjusted_effect <- out$total_effect
    out$background_adjusted_var <- out$total_var
    out$background_adjusted_var_raw <- out$total_var_raw
    out$background_adjusted_df_raw <- out$residual_df
  }

  if (!is.null(reporters)) {
    out <- out[out$reporter %in% reporters, , drop = FALSE]
  }
  if (nrow(out) == 0) {
    return(data.frame())
  }

  global <- stats::aggregate(background_adjusted_effect ~ perturbation, out, mean, na.rm = TRUE)
  names(global)[2] <- "global_effect"
  out <- merge(out, global, by = "perturbation", all.x = TRUE, sort = FALSE)
  out$centering_effect <- out$global_effect
  out$low_rank_effect <- low_rank_background_effect(
    out,
    effect = "background_adjusted_effect",
    rank = fit$background_rank,
    impute = fit$background_impute
  )
  out$rank_adjusted_total_effect <- out$background_adjusted_effect - out$low_rank_effect

  variance_perturbations <- unique(out$perturbation)
  total_variance_split <- split(out$total_var, out$perturbation)
  raw_total_variance_split <- split(out$total_var_raw, out$perturbation)
  evc_variance_split <- split(out$empty_vector_var, out$perturbation)
  background_variance_split <- split(out$background_adjusted_var, out$perturbation)
  variance_summary <- data.frame(
    perturbation = variance_perturbations,
    sum_total_var = vapply(total_variance_split[variance_perturbations], sum, numeric(1), na.rm = TRUE),
    sum_total_var_raw = vapply(raw_total_variance_split[variance_perturbations], sum, numeric(1), na.rm = TRUE),
    sum_background_adjusted_var = vapply(background_variance_split[variance_perturbations], sum, numeric(1), na.rm = TRUE),
    empty_vector_var_for_perturbation = vapply(evc_variance_split[variance_perturbations], function(x) {
      vals <- unique(x[is.finite(x)])
      if (length(vals) == 0) {
        0
      } else {
        vals[1]
      }
    }, numeric(1)),
    n_reporters_for_perturbation = vapply(total_variance_split[variance_perturbations], function(x) sum(is.finite(x)), numeric(1)),
    stringsAsFactors = FALSE
  )
  out <- merge(out, variance_summary, by = "perturbation", all.x = TRUE, sort = FALSE)
  m <- out$n_reporters_for_perturbation
  row_df <- out$residual_df
  out$global_se <- sqrt(out$sum_total_var / m^2 + out$empty_vector_var_for_perturbation)
  global_test <- wald_t_test(out$global_effect, out$global_se, df)
  out$global_statistic <- global_test$statistic
  out$global_pvalue <- global_test$pvalue
  out$rank_adjusted_total_se <- sqrt(out$background_adjusted_var)
  rank_adjusted_total_df <- if (isTRUE(fit$empirical_bayes)) out$total_df_for_test else row_df
  if (!is.null(evc) && nzchar(evc)) {
    # Construct the EVC-adjusted contrast uncertainty before moderating it once.
    out$rank_adjusted_total_se <- sqrt(out$background_adjusted_var_raw)
    rank_adjusted_total_df <- out$background_adjusted_df_raw
    if (isTRUE(fit$empirical_bayes)) {
      adjusted <- eb_moderate_se(out$rank_adjusted_total_se, rank_adjusted_total_df)
      out$rank_adjusted_total_se <- adjusted$se
      rank_adjusted_total_df <- adjusted$df
    }
  }
  out$rank_adjusted_total_df_for_test <- rank_adjusted_total_df
  rank_adjusted_total_test <- wald_t_test(
    out$rank_adjusted_total_effect,
    out$rank_adjusted_total_se,
    rank_adjusted_total_df
  )
  out$rank_adjusted_total_statistic <- rank_adjusted_total_test$statistic
  out$rank_adjusted_total_pvalue <- rank_adjusted_total_test$pvalue
  rank_adjusted_global <- stats::aggregate(rank_adjusted_total_effect ~ perturbation, out, mean, na.rm = TRUE)
  names(rank_adjusted_global)[2] <- "rank_adjusted_global_effect"
  out <- merge(out, rank_adjusted_global, by = "perturbation", all.x = TRUE, sort = FALSE)
  out$specific_effect <- out$rank_adjusted_total_effect - out$rank_adjusted_global_effect
  out$specific_var <- ((m - 1) / m)^2 * out$total_var_raw +
    (out$sum_total_var_raw - out$total_var_raw) / m^2
  out$specific_se <- sqrt(out$specific_var)
  out$specific_se[!is.finite(out$specific_se) | m <= 1] <- NA_real_
  specific_df <- rep(NA_real_, nrow(out))
  for (idx in split(seq_len(nrow(out)), out$perturbation)) {
    specific_df[idx] <- centered_contrast_df(out$total_var_raw[idx], row_df[idx])
  }
  out$specific_df_raw <- specific_df
  if (isTRUE(fit$empirical_bayes)) {
    moderated <- eb_moderate_se(out$specific_se, specific_df)
    out$specific_se <- moderated$se
    specific_df <- moderated$df
  }
  out$specific_df_for_test <- specific_df
  specific_test <- wald_t_test(out$specific_effect, out$specific_se, specific_df)
  out$specific_statistic <- specific_test$statistic
  out$specific_pvalue <- specific_test$pvalue

  out$total_padj_global <- adjust_destress_pvalues(out$total_pvalue, out$reporter, "global")
  out$total_padj_by_reporter <- adjust_destress_pvalues(out$total_pvalue, out$reporter, "by_reporter")
  out$rank_adjusted_total_padj_global <- adjust_destress_pvalues(out$rank_adjusted_total_pvalue, out$reporter, "global")
  out$rank_adjusted_total_padj_by_reporter <- adjust_destress_pvalues(out$rank_adjusted_total_pvalue, out$reporter, "by_reporter")
  out$specific_padj_global <- adjust_destress_pvalues(out$specific_pvalue, out$reporter, "global")
  out$specific_padj_by_reporter <- adjust_destress_pvalues(out$specific_pvalue, out$reporter, "by_reporter")
  adjustment <- if (is.null(fit$adjustment)) "global" else fit$adjustment
  out$total_padj <- if (adjustment == "by_reporter") out$total_padj_by_reporter else if (adjustment == "none") out$total_pvalue else out$total_padj_global
  out$specific_padj <- if (adjustment == "by_reporter") out$specific_padj_by_reporter else if (adjustment == "none") out$specific_pvalue else out$specific_padj_global

  out <- out[, c(
    "reporter", "perturbation",
    "total_effect", "total_se", "total_statistic", "total_pvalue",
    "additive_total_effect", "additive_total_se",
    "empty_vector_effect", "background_adjusted_effect",
    "background_adjusted_var_raw", "background_adjusted_df_raw",
    "global_effect", "global_se", "global_statistic", "global_pvalue",
    "low_rank_effect",
    "rank_adjusted_total_effect", "rank_adjusted_total_se",
    "rank_adjusted_total_statistic", "rank_adjusted_total_pvalue",
    "rank_adjusted_total_df_for_test",
    "rank_adjusted_global_effect",
    "specific_effect", "specific_se", "specific_statistic", "specific_pvalue",
    "total_var_raw", "specific_var", "residual_df",
    "specific_df_raw", "specific_df_for_test",
    "total_padj_global", "total_padj_by_reporter",
    "rank_adjusted_total_padj_global", "rank_adjusted_total_padj_by_reporter",
    "specific_padj_global",
    "specific_padj_by_reporter", "total_padj", "specific_padj"
  ), drop = FALSE]
  out[order(out$reporter, out$perturbation), ]
}

reporter_lm_results <- function(fit, perturbations = NULL, reporters = NULL) {
  control <- fit$assay_info$control
  all_reporters <- fit$levels$reporter
  all_perturbations <- setdiff(fit$levels$perturbation, control)
  reporters <- if (is.null(reporters)) all_reporters else intersect(reporters, all_reporters)
  perturbations <- if (is.null(perturbations)) all_perturbations else intersect(perturbations, all_perturbations)

  if (length(reporters) == 0 || length(perturbations) == 0) {
    return(data.frame())
  }

  rows <- lapply(reporters, function(reporter) {
    reporter_fit <- fit$reporter_fits[[reporter]]
    reporter_data <- fit_model_frame(reporter_fit)
    perturbation_levels <- levels(reporter_data$.perturbation)
    reporter_perturbations <- intersect(perturbations, setdiff(perturbation_levels, control))
    if (length(reporter_perturbations) == 0) {
      return(NULL)
    }

    base <- data.frame(.perturbation = factor(control, levels = perturbation_levels))
    base <- base[rep(1, length(reporter_perturbations)), , drop = FALSE]
    comp <- base
    comp$.perturbation <- factor(reporter_perturbations, levels = perturbation_levels)
    for (col in fit$technical) {
      level <- modal_factor_level(reporter_data[[col]])
      base[[col]] <- factor(level, levels = levels(reporter_data[[col]]))
      comp[[col]] <- base[[col]]
    }

    total <- contrast_estimates(reporter_fit, comp, base)

    data.frame(
      reporter = reporter,
      perturbation = reporter_perturbations,
      total_effect = total[, "estimate"],
      total_se = total[, "std_error"],
      total_statistic = total[, "statistic"],
      total_pvalue = total[, "pvalue"],
      additive_total_effect = NA_real_,
      additive_total_se = NA_real_,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  if (is.null(out) || nrow(out) == 0) {
    return(data.frame())
  }

  global <- stats::aggregate(total_effect ~ perturbation, out, mean, na.rm = TRUE)
  names(global)[2] <- "global_effect"
  out <- merge(out, global, by = "perturbation", all.x = TRUE, sort = FALSE)
  out$low_rank_effect <- low_rank_background_effect(
    out,
    effect = "total_effect",
    rank = fit$background_rank,
    impute = fit$background_impute
  )
  out$rank_adjusted_total_effect <- out$total_effect - out$low_rank_effect

  df <- min(vapply(fit$reporter_fits[reporters], fit_df_residual, numeric(1)), na.rm = TRUE)
  if (isTRUE(fit$empirical_bayes)) {
    moderated <- eb_moderate_se(out$total_se, df)
    out$total_se <- moderated$se
    total_test <- wald_t_test(out$total_effect, out$total_se, moderated$df)
    out$total_statistic <- total_test$statistic
    out$total_pvalue <- total_test$pvalue
  }

  out$total_var <- out$total_se^2
  out$rank_adjusted_total_se <- out$total_se
  rank_adjusted_total_test <- wald_t_test(out$rank_adjusted_total_effect, out$rank_adjusted_total_se, df)
  out$rank_adjusted_total_statistic <- rank_adjusted_total_test$statistic
  out$rank_adjusted_total_pvalue <- rank_adjusted_total_test$pvalue
  variance_summary <- stats::aggregate(
    total_var ~ perturbation,
    out,
    function(x) c(sum = sum(x, na.rm = TRUE), m = sum(is.finite(x)))
  )
  variance_summary <- do.call(data.frame, variance_summary)
  names(variance_summary) <- c("perturbation", "sum_total_var", "n_reporters_for_perturbation")
  out <- merge(out, variance_summary, by = "perturbation", all.x = TRUE, sort = FALSE)

  m <- out$n_reporters_for_perturbation
  out$global_se <- sqrt(out$sum_total_var) / m
  global_test <- wald_t_test(out$global_effect, out$global_se, df)
  out$global_statistic <- global_test$statistic
  out$global_pvalue <- global_test$pvalue
  rank_adjusted_global <- stats::aggregate(rank_adjusted_total_effect ~ perturbation, out, mean, na.rm = TRUE)
  names(rank_adjusted_global)[2] <- "rank_adjusted_global_effect"
  out <- merge(out, rank_adjusted_global, by = "perturbation", all.x = TRUE, sort = FALSE)
  out$specific_effect <- out$rank_adjusted_total_effect - out$rank_adjusted_global_effect
  out$specific_var <- ((m - 1) / m)^2 * out$total_var +
    (out$sum_total_var - out$total_var) / m^2
  out$specific_se <- sqrt(out$specific_var)
  out$specific_se[!is.finite(out$specific_se) | m <= 1] <- NA_real_

  specific_df <- df
  if (isTRUE(fit$empirical_bayes)) {
    moderated <- eb_moderate_se(out$specific_se, df)
    out$specific_se <- moderated$se
    specific_df <- moderated$df
  }
  specific_test <- wald_t_test(out$specific_effect, out$specific_se, specific_df)
  out$specific_statistic <- specific_test$statistic
  out$specific_pvalue <- specific_test$pvalue

  out$total_padj_global <- adjust_destress_pvalues(out$total_pvalue, out$reporter, "global")
  out$total_padj_by_reporter <- adjust_destress_pvalues(out$total_pvalue, out$reporter, "by_reporter")
  out$rank_adjusted_total_padj_global <- adjust_destress_pvalues(out$rank_adjusted_total_pvalue, out$reporter, "global")
  out$rank_adjusted_total_padj_by_reporter <- adjust_destress_pvalues(out$rank_adjusted_total_pvalue, out$reporter, "by_reporter")
  out$specific_padj_global <- adjust_destress_pvalues(out$specific_pvalue, out$reporter, "global")
  out$specific_padj_by_reporter <- adjust_destress_pvalues(out$specific_pvalue, out$reporter, "by_reporter")
  adjustment <- if (is.null(fit$adjustment)) "global" else fit$adjustment
  out$total_padj <- if (adjustment == "by_reporter") out$total_padj_by_reporter else if (adjustment == "none") out$total_pvalue else out$total_padj_global
  out$specific_padj <- if (adjustment == "by_reporter") out$specific_padj_by_reporter else if (adjustment == "none") out$specific_pvalue else out$specific_padj_global

  out <- out[, c(
    "reporter", "perturbation",
    "total_effect", "total_se", "total_statistic", "total_pvalue",
    "additive_total_effect", "additive_total_se",
    "global_effect", "global_se", "global_statistic", "global_pvalue",
    "low_rank_effect",
    "rank_adjusted_total_effect", "rank_adjusted_total_se",
    "rank_adjusted_total_statistic", "rank_adjusted_total_pvalue",
    "rank_adjusted_global_effect",
    "specific_effect", "specific_se", "specific_statistic", "specific_pvalue",
    "total_padj_global", "total_padj_by_reporter",
    "rank_adjusted_total_padj_global", "rank_adjusted_total_padj_by_reporter",
    "specific_padj_global",
    "specific_padj_by_reporter", "total_padj", "specific_padj"
  ), drop = FALSE]
  out[order(out$reporter, out$perturbation), ]
}

observed_mean_results <- function(fit, perturbations = NULL, reporters = NULL) {
  assay <- fit$assay_data
  control <- fit$assay_info$control
  all_reporters <- fit$levels$reporter
  reporters <- if (is.null(reporters)) all_reporters else intersect(reporters, all_reporters)

  d <- assay[assay$.reporter %in% reporters, , drop = FALSE]
  if (!is.null(perturbations)) {
    d <- d[d$.perturbation %in% c(control, perturbations), , drop = FALSE]
  }
  if (nrow(d) == 0) {
    return(data.frame())
  }

  d$.adjusted_response <- technical_adjusted_response(fit$total_fit, d, fit$technical)
  cell_mean <- stats::aggregate(
    .adjusted_response ~ .reporter + .perturbation,
    d,
    function(x) c(mean = mean(x, na.rm = TRUE), n = sum(is.finite(x)))
  )
  cell_mean <- do.call(data.frame, cell_mean)
  names(cell_mean) <- c("reporter", "perturbation", "mean_response", "n")
  cell_mean$reporter <- as.character(cell_mean$reporter)
  cell_mean$perturbation <- as.character(cell_mean$perturbation)
  cell_mean$n <- as.numeric(cell_mean$n)

  cell_key <- paste(cell_mean$reporter, cell_mean$perturbation, sep = "\r")
  mean_by_cell <- stats::setNames(cell_mean$mean_response, cell_key)
  d_key <- paste(as.character(d$.reporter), as.character(d$.perturbation), sep = "\r")
  within_residual <- d$.adjusted_response - mean_by_cell[d_key]
  residual_df <- sum(is.finite(within_residual)) - nrow(cell_mean)
  if (!is.finite(residual_df) || residual_df <= 0) {
    residual_df <- fit_df_residual(fit$total_fit)
  }
  residual_sigma <- fit_sigma(fit$total_fit)

  control_mean <- cell_mean[cell_mean$perturbation == control, c("reporter", "mean_response", "n"), drop = FALSE]
  names(control_mean) <- c("reporter", "control_mean_response", "control_n")

  out <- merge(
    cell_mean[cell_mean$perturbation != control, , drop = FALSE],
    control_mean,
    by = "reporter",
    all.x = TRUE,
    sort = FALSE
  )
  out <- out[is.finite(out$mean_response) & is.finite(out$control_mean_response), , drop = FALSE]
  if (!is.null(perturbations)) {
    out <- out[out$perturbation %in% perturbations, , drop = FALSE]
  }
  if (nrow(out) == 0) {
    return(data.frame())
  }

  out$total_effect <- out$mean_response - out$control_mean_response
  out$total_var <- residual_sigma^2 * (1 / out$n + 1 / out$control_n)
  global <- stats::aggregate(total_effect ~ perturbation, out, mean, na.rm = TRUE)
  names(global)[2] <- "global_effect"
  out <- merge(out, global, by = "perturbation", all.x = TRUE, sort = FALSE)
  out$low_rank_effect <- low_rank_background_effect(
    out,
    effect = "total_effect",
    rank = fit$background_rank,
    impute = fit$background_impute
  )
  out$rank_adjusted_total_effect <- out$total_effect - out$low_rank_effect

  variance_summary <- stats::aggregate(
    total_var ~ perturbation,
    out,
    function(x) c(sum = sum(x, na.rm = TRUE), m = sum(is.finite(x)))
  )
  variance_summary <- do.call(data.frame, variance_summary)
  names(variance_summary) <- c("perturbation", "sum_total_var", "n_reporters_for_perturbation")
  out <- merge(out, variance_summary, by = "perturbation", all.x = TRUE, sort = FALSE)

  out$total_se <- sqrt(out$total_var)
  out$rank_adjusted_total_se <- out$total_se
  rank_adjusted_total_test <- wald_t_test(out$rank_adjusted_total_effect, out$rank_adjusted_total_se, residual_df)
  out$rank_adjusted_total_statistic <- rank_adjusted_total_test$statistic
  out$rank_adjusted_total_pvalue <- rank_adjusted_total_test$pvalue
  m <- out$n_reporters_for_perturbation
  rank_adjusted_global <- stats::aggregate(rank_adjusted_total_effect ~ perturbation, out, mean, na.rm = TRUE)
  names(rank_adjusted_global)[2] <- "rank_adjusted_global_effect"
  out <- merge(out, rank_adjusted_global, by = "perturbation", all.x = TRUE, sort = FALSE)
  out$specific_effect <- out$rank_adjusted_total_effect - out$rank_adjusted_global_effect
  out$specific_var <- ((m - 1) / m)^2 * out$total_var +
    (out$sum_total_var - out$total_var) / m^2
  out$specific_se <- sqrt(out$specific_var)
  out$specific_se[!is.finite(out$specific_se) | m <= 1] <- NA_real_
  specific_df <- residual_df
  if (isTRUE(fit$empirical_bayes)) {
    moderated <- eb_moderate_se(out$specific_se, residual_df)
    out$specific_se <- moderated$se
    specific_df <- moderated$df
  }
  out$specific_statistic <- out$specific_effect / out$specific_se
  out$specific_pvalue <- 2 * stats::pt(
    abs(out$specific_statistic),
    df = specific_df,
    lower.tail = FALSE
  )
  out$total_statistic <- out$total_effect / out$total_se
  out$total_pvalue <- 2 * stats::pt(
    abs(out$total_statistic),
    df = residual_df,
    lower.tail = FALSE
  )
  out$additive_total_effect <- NA_real_
  out$additive_total_se <- NA_real_

  out$total_padj_global <- adjust_destress_pvalues(out$total_pvalue, out$reporter, "global")
  out$total_padj_by_reporter <- adjust_destress_pvalues(out$total_pvalue, out$reporter, "by_reporter")
  out$rank_adjusted_total_padj_global <- adjust_destress_pvalues(out$rank_adjusted_total_pvalue, out$reporter, "global")
  out$rank_adjusted_total_padj_by_reporter <- adjust_destress_pvalues(out$rank_adjusted_total_pvalue, out$reporter, "by_reporter")
  out$specific_padj_global <- adjust_destress_pvalues(out$specific_pvalue, out$reporter, "global")
  out$specific_padj_by_reporter <- adjust_destress_pvalues(out$specific_pvalue, out$reporter, "by_reporter")
  adjustment <- if (is.null(fit$adjustment)) "global" else fit$adjustment
  out$total_padj <- if (adjustment == "by_reporter") out$total_padj_by_reporter else if (adjustment == "none") out$total_pvalue else out$total_padj_global
  out$specific_padj <- if (adjustment == "by_reporter") out$specific_padj_by_reporter else if (adjustment == "none") out$specific_pvalue else out$specific_padj_global

  out <- out[, c(
    "reporter", "perturbation",
    "total_effect", "total_se", "total_statistic", "total_pvalue",
    "additive_total_effect", "additive_total_se",
    "global_effect", "low_rank_effect",
    "rank_adjusted_total_effect", "rank_adjusted_total_se",
    "rank_adjusted_total_statistic", "rank_adjusted_total_pvalue",
    "rank_adjusted_global_effect",
    "specific_effect", "specific_se", "specific_statistic",
    "specific_pvalue",
    "total_padj_global", "total_padj_by_reporter",
    "rank_adjusted_total_padj_global", "rank_adjusted_total_padj_by_reporter",
    "specific_padj_global",
    "specific_padj_by_reporter", "total_padj", "specific_padj"
  ), drop = FALSE]
  out[order(out$reporter, out$perturbation), ]
}

#' List available DStressR presets
#'
#' @return A character vector of preset names accepted by [fit_destress()].
#' @export
destress_presets <- function() {
  c("model", "median_polish", "empty_vector_control")
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
    median_polish_legacy = "median_polish",
    medianpolish = "median_polish",
    legacy = "median_polish",
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
    median_polish = list(
      normalization = "median_polish",
      testing = "gaussian_z",
      aggregation = "max_p",
      adjustment = "by_reporter"
    ),
    empty_vector_control = list(
      normalization = "empty_vector",
      testing = "gaussian_z",
      aggregation = "max_p",
      adjustment = "by_reporter"
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
      choices = c("global", "by_reporter", "none"),
      aliases = c(reporter = "by_reporter", within_reporter = "by_reporter"),
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
    if (stages$testing != "gaussian_z" || stages$aggregation != "max_p" || stages$adjustment != "by_reporter") {
      stop(
        "The ", expected, " compatibility path currently requires ",
        "`testing = \"gaussian_z\"`, `aggregation = \"max_p\"`, and ",
        "`adjustment = \"by_reporter\"`.",
        call. = FALSE
      )
    }
  }
  stages
}

#' Fit DStressR with staged statistical options
#'
#' `fit_destress()` is the main DStressR entry point. By default it fits the
#' model-based analysis, but it can also run named presets for the
#' screen-specific median-polish and empty-vector-control analyses.
#'
#' The staged options make the major statistical choices explicit:
#' normalization, test statistic and p-value calculation, replicate aggregation,
#' and p-value adjustment. Only implemented combinations are accepted. For the
#' model-based path, growth-response normalization is performed upstream by
#' [prepare_assay()], where `growth_exponent` can be fixed, estimated, or
#' supplied as reporter-specific values.
#'
#' @param assay A `destress_assay` produced by [prepare_assay()] or a raw assay
#'   data frame for `normalization = "linear_model"`, or a long expression
#'   table for the compatibility presets.
#' @param technical Character vector of batch, plate, replicate, or other
#'   technical-factor columns to include.
#' @param empirical_bayes If `TRUE`, moderates standard errors toward an
#'   empirical prior variance with prior degrees of freedom estimated from the
#'   observed variance distribution, using each contrast's degrees of freedom.
#'   Reporter-wise specific effects combine raw total-effect variances with
#'   Welch-Satterthwaite degrees of freedom and are moderated once. This maps to `testing = "moderated_t"` for
#'   the model path; `FALSE` maps to `testing = "student_t"`.
#' @param empty_vector_reporter Optional reporter/control strain used as an
#'   empty-vector reporter in the model-based path. When supplied, its
#'   reference-relative perturbation effect is subtracted from every reporter's
#'   reference-relative perturbation effect before reporter-library centering.
#'   Its raw contrast variance is added to the reporter's raw contrast variance,
#'   with Welch-Satterthwaite degrees of freedom; the combined contrast variance
#'   is moderated once when requested. This effect-level adjustment is an
#'   alternative to response-level background calibration in [prepare_assay()].
#'   The two routes cannot be combined. Supported for `interaction = FALSE`.
#' @param background_rank Non-negative integer. The default `0` removes no
#'   latent background. Values `1` or `2` additionally subtract a low-rank
#'   background term from the reference-relative total-effect matrix before
#'   testing rank-adjusted total and reporter-specific residual effects. Use
#'   `"auto"` to choose the rank by permutation parallel analysis.
#' @param background_rank_max Maximum rank considered when
#'   `background_rank = "auto"`.
#' @param background_rank_permutations Number of permutations used when
#'   `background_rank = "auto"`.
#' @param background_rank_threshold Permutation reference quantile used when
#'   `background_rank = "auto"`.
#' @param background_rank_seed Optional random seed used when
#'   `background_rank = "auto"`.
#' @param background_impute Method used to fill missing entries of the effect
#'   matrix before singular-value decomposition. The default `"column_mean"`
#'   replaces missing entries by the observed mean of the corresponding
#'   perturbation column. `"global_mean"` and `"zero"` are also available.
#' @param normalization One of `"linear_model"`, `"median_polish"`, or
#'   `"empty_vector"`. `"model"` and `"evc"` are accepted aliases.
#' @param testing One of `"student_t"`, `"moderated_t"`, or `"gaussian_z"`.
#' @param aggregation One of `"none"` or `"max_p"`.
#' @param adjustment One of `"global"`, `"by_reporter"`, or `"none"`.
#' @param interaction If `FALSE`, fit one Gaussian linear
#'   model per reporter with the control perturbation as reference and the supplied
#'   technical covariates as design terms. The latter is the scalable path for
#'   reporter-specific perturbation effects. If `TRUE`, fit the full
#'   reporter-by-perturbation interaction model.
#' @param preset Optional named preset: `"model"`, `"median_polish"`, or
#'   `"empty_vector_control"`. The model-based workflow remains the default
#'   when omitted. `"median_polish_legacy"` is accepted as an alias for
#'   `"median_polish"`, and `"evc"` as an alias for `"empty_vector_control"`.
#' @param verbose If `TRUE`, print stage-level progress messages during model
#'   fitting. This is useful for large screens and full interaction models.
#' @param ... For `normalization = "linear_model"` with a raw data frame,
#'   arguments passed to [prepare_assay()], including `growth_exponent`. For
#'   compatibility presets, arguments passed to the selected engine.
#' @return A fitted DStressR object. The model path returns a `destress_fit`;
#'   screen-specific presets return their corresponding workflow result objects.
#' @export
fit_destress <- function(assay,
                         technical = NULL,
                         empirical_bayes = TRUE,
                         empty_vector_reporter = NULL,
                         background_rank = 0,
                         background_rank_max = 5,
                         background_rank_permutations = 100,
                         background_rank_threshold = 0.99,
                         background_rank_seed = NULL,
                         background_impute = c("column_mean", "global_mean", "zero"),
                         normalization = NULL,
                         testing = NULL,
                         aggregation = NULL,
                         adjustment = NULL,
                         interaction = FALSE,
                         preset = NULL,
                         verbose = FALSE,
                         ...) {
  destress_progress(verbose, "DStressR: resolving staged model options.")
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
    destress_progress(verbose, "DStressR: running median-polish preset.")
    fit <- fit_median_polish(assay, ...)
    attr(fit, "destress_preset") <- if (is.null(preset)) "median_polish" else preset
    attr(fit, "destress_stages") <- stages
    return(fit)
  }
  if (stages$normalization == "empty_vector") {
    destress_progress(verbose, "DStressR: running empty-vector-control preset.")
    fit <- fit_empty_vector_control(assay, ...)
    attr(fit, "destress_preset") <- if (is.null(preset)) "empty_vector_control" else preset
    attr(fit, "destress_stages") <- stages
    return(fit)
  }

  empirical_bayes <- identical(stages$testing, "moderated_t")
  background_rank_request <- validate_background_rank_request(background_rank)
  background_rank_max <- validate_background_rank(background_rank_max)
  background_rank_permutations <- validate_background_rank(background_rank_permutations)
  background_rank_threshold <- validate_background_rank_threshold(background_rank_threshold)
  background_impute <- validate_background_impute(background_impute)
  if (!inherits(assay, "destress_assay")) {
    if (!is.data.frame(assay)) {
      stop("`assay` must be a data frame or be produced by prepare_assay().", call. = FALSE)
    }
    destress_progress(verbose, "DStressR: preparing assay object.")
    assay <- prepare_assay(assay, ...)
  }
  background_reporter <- attr(assay, "destress")$background_reporter
  if (!is.null(empty_vector_reporter)) {
    empty_vector_reporter <- as.character(empty_vector_reporter)
    if (length(empty_vector_reporter) != 1 || is.na(empty_vector_reporter) || !nzchar(empty_vector_reporter)) {
      stop("`empty_vector_reporter` must be one reporter label.", call. = FALSE)
    }
    if (!is.null(background_reporter) && nzchar(background_reporter)) {
      stop(
        "Conflicting EVC settings: response-level background calibration in ",
        "prepare_assay() cannot be combined with `empty_vector_reporter` in ",
        "fit_destress(). Choose response-level calibration/subtraction or ",
        "fitted-effect subtraction, not both.",
        call. = FALSE
      )
    }
    if (isTRUE(interaction)) {
      stop("`empty_vector_reporter` requires `interaction = FALSE`.", call. = FALSE)
    }
    if (!empty_vector_reporter %in% levels(assay$.reporter)) {
      stop("Empty-vector reporter '", empty_vector_reporter, "' was not found in the assay.", call. = FALSE)
    }
  }
  if (!is.null(background_reporter) && nzchar(background_reporter)) {
    assay_for_model <- assay[as.character(assay$.reporter) != background_reporter, , drop = FALSE]
    assay_for_model$.reporter <- droplevels(assay_for_model$.reporter)
    assay_for_model$.perturbation <- droplevels(assay_for_model$.perturbation)
  } else {
    assay_for_model <- assay
  }
  technical <- technical[!is.na(technical) & nzchar(technical)]
  missing_technical <- setdiff(technical, names(assay))
  if (length(missing_technical) > 0) {
    stop("Unknown technical columns: ", paste(missing_technical, collapse = ", "), call. = FALSE)
  }
  interaction <- isTRUE(interaction)
  destress_progress(
    verbose,
    "DStressR: fitting ",
    if (interaction) "full reporter-by-perturbation interaction model" else "scalable reporter-specific model",
    " with ",
    nlevels(assay_for_model$.reporter),
    " reporters, ",
    max(0, nlevels(assay_for_model$.perturbation) - 1),
    " non-reference perturbations, and ",
    nrow(assay_for_model),
    " observations."
  )
  if (isTRUE(background_rank_request$auto) && interaction) {
    stop("`background_rank = \"auto\"` is currently supported for the scalable reporter-specific model only.", call. = FALSE)
  }
  formulas <- make_formulas(technical)
  total_fit <- if (interaction) {
    destress_progress(verbose, "DStressR: fitting additive linear model.")
    stats::lm(formulas$total, data = assay_for_model, na.action = stats::na.exclude)
  } else {
    NULL
  }
  full_fit <- if (interaction) {
    destress_progress(verbose, "DStressR: fitting full interaction linear model.")
    stats::lm(formulas$full, data = assay_for_model, na.action = stats::na.exclude)
  } else {
    NULL
  }
  assay_data <- if (interaction) {
    NULL
  } else {
    assay_for_model
  }
  reporter_formula <- stats::as.formula(paste(".response ~ .perturbation +", technical_formula(technical)))
  reporter_fits <- if (interaction) {
    NULL
  } else {
    NULL
  }
  reporter_effects <- if (interaction) {
    NULL
  } else {
    destress_progress(verbose, "DStressR: estimating reporter-specific perturbation effects.")
    fit_reporter_effects(
      assay_for_model,
      technical,
      attr(assay, "destress")$control,
      verbose = verbose
    )
  }
  background_rank_diagnostics_table <- NULL
  selected_background_rank <- background_rank_request$rank
  if (isTRUE(background_rank_request$auto)) {
    destress_progress(verbose, "DStressR: selecting low-rank background by permutation diagnostics.")
    background_rank_diagnostics_table <- background_rank_diagnostics(
      reporter_effects,
      effect = "total_effect",
      reporter = "reporter",
      perturbation = "perturbation",
      rank_max = background_rank_max,
      permutations = background_rank_permutations,
      threshold = background_rank_threshold,
      impute = background_impute,
      seed = background_rank_seed
    )
    selected_background_rank <- select_background_rank(background_rank_diagnostics_table)
    destress_progress(verbose, "DStressR: selected low-rank background rank ", selected_background_rank, ".")
  }

  destress_progress(verbose, "DStressR: assembling fitted object.")
  structure(
    list(
      total_fit = total_fit,
      full_fit = full_fit,
      interaction = interaction,
      assay_data = assay_data,
      reporter_fits = reporter_fits,
      reporter_effects = reporter_effects,
      growth_exponents = attr(assay, "destress")$growth_exponent_fit,
      assay_info = attr(assay, "destress"),
      levels = list(
        reporter = levels(assay_for_model$.reporter),
        perturbation = levels(assay_for_model$.perturbation)
      ),
      technical = technical,
      empirical_bayes = empirical_bayes,
      empty_vector_reporter = empty_vector_reporter,
      background_reporter = background_reporter,
      background_rank = selected_background_rank,
      background_rank_requested = background_rank,
      background_rank_diagnostics = background_rank_diagnostics_table,
      background_rank_threshold = background_rank_threshold,
      background_rank_seed = background_rank_seed,
      background_impute = background_impute,
      stages = stages,
      preset = if (is.null(preset)) "model" else preset,
      adjustment = stages$adjustment
    ),
    class = "destress_fit"
  )
}

#' Extract estimated model parameters
#'
#' @param fit A `destress_fit` or `destress_summary_fit` object.
#' @return A named list of estimated parameter tables available for the fitted
#'   model. The scalable model path includes reporter-specific growth
#'   normalization estimates and reporter-perturbation effect estimates.
#' @export
model_parameters <- function(fit) {
  if (inherits(fit, "destress_summary_fit")) {
    out <- list(
      background = data.frame(
        background_rank = validate_background_rank(fit$background_rank),
        background_rank_requested = if (is.null(fit$background_rank_requested)) fit$background_rank else as.character(fit$background_rank_requested),
        background_rank_threshold = if (is.null(fit$background_rank_threshold)) NA_real_ else fit$background_rank_threshold,
        background_impute = if (is.null(fit$background_impute)) "column_mean" else fit$background_impute
      ),
      summary_level = data.frame(
        reference_perturbation = fit$reference_perturbation,
        inference = fit$inference,
        global_method = fit$global_method,
        adjustment = fit$adjustment
      ),
      effects = fit$effects
    )
    if (!is.null(fit$background_rank_diagnostics)) {
      out$background_rank_diagnostics <- fit$background_rank_diagnostics
    }
    return(out)
  }
  if (!inherits(fit, "destress_fit")) {
    stop("`fit` must be a destress_fit or destress_summary_fit.", call. = FALSE)
  }

  out <- list(
    background = data.frame(
      background_rank = validate_background_rank(fit$background_rank),
      background_rank_requested = if (is.null(fit$background_rank_requested)) fit$background_rank else as.character(fit$background_rank_requested),
      background_rank_threshold = if (is.null(fit$background_rank_threshold)) NA_real_ else fit$background_rank_threshold,
      background_impute = if (is.null(fit$background_impute)) "column_mean" else fit$background_impute,
      background_reporter = if (is.null(fit$background_reporter)) NA_character_ else fit$background_reporter,
      background_method = if (is.null(fit$assay_info$background_method)) "none" else fit$assay_info$background_method
    ),
    growth_exponents = fit$growth_exponents,
    background_calibration = fit$assay_info$background_fit,
    reporter_effects = fit$reporter_effects
  )
  if (!is.null(fit$background_rank_diagnostics)) {
    out$background_rank_diagnostics <- fit$background_rank_diagnostics
  }

  if (isTRUE(fit$interaction)) {
    coef_table <- function(model) {
      coefs <- summary(model)$coefficients
      data.frame(
        term = rownames(coefs),
        estimate = coefs[, "Estimate"],
        std_error = coefs[, "Std. Error"],
        statistic = coefs[, "t value"],
        pvalue = coefs[, "Pr(>|t|)"],
        row.names = NULL,
        check.names = FALSE
      )
    }
    out$additive_coefficients <- coef_table(fit$total_fit)
    out$interaction_coefficients <- coef_table(fit$full_fit)
  }

  out
}

adjust_destress_pvalues <- function(pvalue, groups = NULL, adjustment = "global") {
  out <- rep(NA_real_, length(pvalue))
  finite <- is.finite(pvalue)
  if (adjustment == "none") {
    out[finite] <- pvalue[finite]
  } else if (adjustment == "global") {
    out[finite] <- stats::p.adjust(pvalue[finite], method = "BH")
  } else if (adjustment == "by_reporter") {
    split_idx <- split(seq_along(pvalue), groups)
    for (idx in split_idx) {
      finite_idx <- idx[is.finite(pvalue[idx])]
      out[finite_idx] <- stats::p.adjust(pvalue[finite_idx], method = "BH")
    }
  }
  out
}

summary_level_effects <- function(d,
                                  reference_perturbation,
                                  reporter,
                                  perturbation,
                                  response,
                                  response_variance,
                                  global_method,
                                  background_rank,
                                  background_impute,
                                  adjustment) {
  d <- data.frame(
    reporter = as.character(d[[reporter]]),
    perturbation = as.character(d[[perturbation]]),
    response = as.numeric(d[[response]]),
    response_variance = as.numeric(d[[response_variance]]),
    stringsAsFactors = FALSE
  )
  d$response_variance[!is.finite(d$response_variance) | d$response_variance < 0] <- NA_real_

  ref <- d[d$perturbation == reference_perturbation, c(
    "reporter", "response", "response_variance"
  )]
  names(ref) <- c("reporter", "reference_response", "reference_response_variance")
  if (anyDuplicated(ref$reporter)) {
    stop("Each reporter must have at most one reference response.", call. = FALSE)
  }
  out <- merge(d, ref, by = "reporter", all.x = TRUE, sort = FALSE)
  out <- out[out$perturbation != reference_perturbation, , drop = FALSE]
  out$total_effect <- out$response - out$reference_response
  out$total_variance <- out$response_variance + out$reference_response_variance
  out$total_se <- sqrt(out$total_variance)
  out$total_statistic <- out$total_effect / out$total_se
  out$total_pvalue <- 2 * stats::pnorm(-abs(out$total_statistic))

  compute_center <- function(effect, variance) {
    finite <- is.finite(effect) & is.finite(variance) & variance > 0
    if (!any(finite)) {
      return(list(center = NA_real_, center_variance = NA_real_, weights = rep(NA_real_, length(effect))))
    }
    weights <- rep(NA_real_, length(effect))
    if (global_method == "weighted_mean") {
      w <- 1 / variance[finite]
      w <- w / sum(w)
      weights[finite] <- w
      center <- sum(w * effect[finite])
      center_variance <- sum((w^2) * variance[finite])
    } else if (global_method == "mean") {
      w <- rep(1 / sum(finite), sum(finite))
      weights[finite] <- w
      center <- mean(effect[finite])
      center_variance <- sum((w^2) * variance[finite])
    } else {
      stop("Unknown global method.", call. = FALSE)
    }
    list(center = center, center_variance = center_variance, weights = weights)
  }

  split_idx <- split(seq_len(nrow(out)), out$perturbation)
  out$global_effect <- NA_real_
  out$global_variance <- NA_real_
  out$global_weight <- NA_real_
  out$specific_effect <- NA_real_
  out$specific_variance <- NA_real_
  for (idx in split_idx) {
    center <- compute_center(out$total_effect[idx], out$total_variance[idx])
    out$global_effect[idx] <- center$center
    out$global_variance[idx] <- center$center_variance
    out$global_weight[idx] <- center$weights
    out$specific_effect[idx] <- out$total_effect[idx] - center$center
    out$specific_variance[idx] <-
      out$total_variance[idx] - 2 * center$weights * out$total_variance[idx] +
      center$center_variance
  }
  out$specific_se <- sqrt(out$specific_variance)
  out$specific_statistic <- out$specific_effect / out$specific_se
  out$specific_pvalue <- 2 * stats::pnorm(-abs(out$specific_statistic))

  out$low_rank_effect <- low_rank_background_effect(
    out,
    effect = "total_effect",
    rank = background_rank,
    impute = background_impute
  )
  out$rank_adjusted_total_effect <- out$total_effect - out$low_rank_effect
  out$rank_adjusted_total_variance <- out$total_variance
  out$rank_adjusted_total_se <- out$total_se
  out$rank_adjusted_total_statistic <-
    out$rank_adjusted_total_effect / out$rank_adjusted_total_se
  out$rank_adjusted_total_pvalue <-
    2 * stats::pnorm(-abs(out$rank_adjusted_total_statistic))

  out$rank_adjusted_global_effect <- NA_real_
  out$rank_adjusted_global_variance <- NA_real_
  out$rank_adjusted_global_weight <- NA_real_
  out$rank_adjusted_specific_effect <- NA_real_
  out$rank_adjusted_specific_variance <- NA_real_
  for (idx in split_idx) {
    center <- compute_center(out$rank_adjusted_total_effect[idx], out$total_variance[idx])
    out$rank_adjusted_global_effect[idx] <- center$center
    out$rank_adjusted_global_variance[idx] <- center$center_variance
    out$rank_adjusted_global_weight[idx] <- center$weights
    out$rank_adjusted_specific_effect[idx] <-
      out$rank_adjusted_total_effect[idx] - center$center
    out$rank_adjusted_specific_variance[idx] <-
      out$total_variance[idx] - 2 * center$weights * out$total_variance[idx] +
      center$center_variance
  }
  out$rank_adjusted_specific_se <- sqrt(out$rank_adjusted_specific_variance)
  out$rank_adjusted_specific_statistic <-
    out$rank_adjusted_specific_effect / out$rank_adjusted_specific_se
  out$rank_adjusted_specific_pvalue <-
    2 * stats::pnorm(-abs(out$rank_adjusted_specific_statistic))

  out$total_padj_global <- adjust_destress_pvalues(out$total_pvalue, out$reporter, "global")
  out$total_padj_by_reporter <- adjust_destress_pvalues(out$total_pvalue, out$reporter, "by_reporter")
  out$specific_padj_global <- adjust_destress_pvalues(out$specific_pvalue, out$reporter, "global")
  out$specific_padj_by_reporter <- adjust_destress_pvalues(out$specific_pvalue, out$reporter, "by_reporter")
  out$rank_adjusted_total_padj_global <-
    adjust_destress_pvalues(out$rank_adjusted_total_pvalue, out$reporter, "global")
  out$rank_adjusted_total_padj_by_reporter <-
    adjust_destress_pvalues(out$rank_adjusted_total_pvalue, out$reporter, "by_reporter")
  out$rank_adjusted_specific_padj_global <-
    adjust_destress_pvalues(out$rank_adjusted_specific_pvalue, out$reporter, "global")
  out$rank_adjusted_specific_padj_by_reporter <-
    adjust_destress_pvalues(out$rank_adjusted_specific_pvalue, out$reporter, "by_reporter")

  out$total_padj <- if (adjustment == "by_reporter") {
    out$total_padj_by_reporter
  } else if (adjustment == "none") {
    out$total_pvalue
  } else {
    out$total_padj_global
  }
  out$specific_padj <- if (adjustment == "by_reporter") {
    out$specific_padj_by_reporter
  } else if (adjustment == "none") {
    out$specific_pvalue
  } else {
    out$specific_padj_global
  }
  out$rank_adjusted_total_padj <- if (adjustment == "by_reporter") {
    out$rank_adjusted_total_padj_by_reporter
  } else if (adjustment == "none") {
    out$rank_adjusted_total_pvalue
  } else {
    out$rank_adjusted_total_padj_global
  }
  out$rank_adjusted_specific_padj <- if (adjustment == "by_reporter") {
    out$rank_adjusted_specific_padj_by_reporter
  } else if (adjustment == "none") {
    out$rank_adjusted_specific_pvalue
  } else {
    out$rank_adjusted_specific_padj_global
  }
  out[order(out$reporter, out$perturbation), , drop = FALSE]
}

#' Fit DStressR from supplied response estimates and variances
#'
#' `fit_destress_summary()` analyzes data where each reporter-perturbation
#' combination is already summarized by a response estimate and an estimated
#' sampling variance. This is useful for public datasets that release processed
#' promoter activities with uncertainty estimates, but not replicate-level wells.
#' The resulting p-values are summary-level Gaussian approximations rather than
#' replicate-level moderated DStressR tests.
#'
#' @param data A data frame with one row per reporter and perturbation.
#' @param reporter,perturbation Columns identifying reporters and perturbations.
#' @param response Numeric response estimate column.
#' @param response_variance Numeric response-variance column. This is the
#'   primary uncertainty input.
#' @param response_se Optional response-standard-error column. Used only when
#'   `response_variance` is `NULL`, in which case it is squared internally.
#' @param reference_perturbation Reference perturbation/condition.
#' @param global_method How to estimate perturbation-level global effects.
#'   `"weighted_mean"` uses inverse-variance weights; `"mean"` uses equal
#'   weights.
#' @param background_rank Non-negative integer or `"auto"` for permutation-based
#'   rank selection.
#' @inheritParams fit_destress
#' @return A `destress_summary_fit` object. Use [results()] to extract the
#'   effect table.
#' @export
fit_destress_summary <- function(data,
                                 reporter = "reporter",
                                 perturbation = "perturbation",
                                 response = "response",
                                 response_variance = "response_variance",
                                 response_se = NULL,
                                 reference_perturbation,
                                 global_method = c("weighted_mean", "mean"),
                                 background_rank = 0,
                                 background_rank_max = 5,
                                 background_rank_permutations = 100,
                                 background_rank_threshold = 0.99,
                                 background_rank_seed = NULL,
                                 background_impute = c("column_mean", "global_mean", "zero"),
                                 adjustment = c("global", "by_reporter", "none")) {
  stopifnot(is.data.frame(data))
  global_method <- match.arg(global_method)
  adjustment <- match.arg(adjustment)
  background_rank_request <- validate_background_rank_request(background_rank)
  background_rank_max <- validate_background_rank(background_rank_max)
  background_rank_permutations <- validate_background_rank(background_rank_permutations)
  background_rank_threshold <- validate_background_rank_threshold(background_rank_threshold)
  background_impute <- validate_background_impute(background_impute)
  if (missing(reference_perturbation) || length(reference_perturbation) != 1 ||
      is.na(reference_perturbation) || !nzchar(as.character(reference_perturbation))) {
    stop("`reference_perturbation` must be one non-empty label.", call. = FALSE)
  }
  required <- c(reporter, perturbation, response)
  if (is.null(response_variance)) {
    if (is.null(response_se)) {
      stop("Supply `response_variance`, or supply `response_se` and set `response_variance = NULL`.", call. = FALSE)
    }
    required <- c(required, response_se)
  } else {
    required <- c(required, response_variance)
  }
  missing_cols <- setdiff(required, names(data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  d <- data
  if (is.null(response_variance)) {
    d$.destress_response_variance <- as.numeric(d[[response_se]])^2
    response_variance <- ".destress_response_variance"
  }

  initial_rank <- if (isTRUE(background_rank_request$auto)) 0L else background_rank_request$rank
  effects0 <- summary_level_effects(
    d,
    reference_perturbation = as.character(reference_perturbation),
    reporter = reporter,
    perturbation = perturbation,
    response = response,
    response_variance = response_variance,
    global_method = global_method,
    background_rank = initial_rank,
    background_impute = background_impute,
    adjustment = adjustment
  )
  background_rank_diagnostics_table <- NULL
  selected_background_rank <- initial_rank
  if (isTRUE(background_rank_request$auto)) {
    background_rank_diagnostics_table <- background_rank_diagnostics(
      effects0,
      effect = "total_effect",
      reporter = "reporter",
      perturbation = "perturbation",
      rank_max = background_rank_max,
      permutations = background_rank_permutations,
      threshold = background_rank_threshold,
      impute = background_impute,
      seed = background_rank_seed
    )
    selected_background_rank <- select_background_rank(background_rank_diagnostics_table)
  }
  effects <- if (selected_background_rank == initial_rank) {
    effects0
  } else {
    summary_level_effects(
      d,
      reference_perturbation = as.character(reference_perturbation),
      reporter = reporter,
      perturbation = perturbation,
      response = response,
      response_variance = response_variance,
      global_method = global_method,
      background_rank = selected_background_rank,
      background_impute = background_impute,
      adjustment = adjustment
    )
  }

  structure(
    list(
      data = data,
      effects = effects,
      reference_perturbation = as.character(reference_perturbation),
      columns = list(
        reporter = reporter,
        perturbation = perturbation,
        response = response,
        response_variance = response_variance,
        response_se = response_se
      ),
      inference = "summary_level_gaussian",
      global_method = global_method,
      background_rank = selected_background_rank,
      background_rank_requested = background_rank,
      background_rank_diagnostics = background_rank_diagnostics_table,
      background_rank_threshold = background_rank_threshold,
      background_rank_seed = background_rank_seed,
      background_impute = background_impute,
      adjustment = adjustment
    ),
    class = "destress_summary_fit"
  )
}

#' Extract model results
#'
#' @param fit A `destress_fit` or `destress_summary_fit` object.
#' @param perturbations Optional perturbation subset.
#' @param reporters Optional reporter subset.
#' @param verbose If `TRUE`, print stage-level progress messages while
#'   extracting model results.
#' @return A data frame with total and reporter-specific effects. Reporter-wise
#'   results also include `total_var_raw` and `specific_var` (before moderation),
#'   `residual_df`, `specific_df_raw` (Welch-Satterthwaite effective degrees of
#'   freedom), and `specific_df_for_test` (including prior degrees of freedom
#'   when moderation is enabled). Full interaction-model results include
#'   `total_var_raw`, `residual_df`, and `total_df_for_test`. Their
#'   `specific_var` uses the centered contrast and full coefficient covariance;
#'   `specific_df_raw` and `specific_df_for_test` record its ordinary and
#'   moderated degrees of freedom. Both targets honor `empirical_bayes`.
#' @export
results <- function(fit, perturbations = NULL, reporters = NULL, verbose = FALSE) {
  if (inherits(fit, "destress_summary_fit")) {
    destress_progress(verbose, "DStressR: extracting summary-level results.")
    out <- fit$effects
    if (!is.null(perturbations)) {
      out <- out[out$perturbation %in% perturbations, , drop = FALSE]
    }
    if (!is.null(reporters)) {
      out <- out[out$reporter %in% reporters, , drop = FALSE]
    }
    return(out[order(out$reporter, out$perturbation), , drop = FALSE])
  }
  if (!inherits(fit, "destress_fit")) {
    stop("`fit` must be a destress_fit or destress_summary_fit.", call. = FALSE)
  }
  if (!isTRUE(fit$interaction)) {
    destress_progress(verbose, "DStressR: extracting scalable reporter-specific results.")
    return(reporter_effect_results(fit, perturbations = perturbations, reporters = reporters))
  }
  if (!is.null(fit$empty_vector_reporter)) {
    stop("Model-based empty-vector adjustment is currently implemented for the scalable reporter-specific path.", call. = FALSE)
  }
  control <- fit$assay_info$control
  all_reporters <- fit$levels$reporter
  all_perturbations <- setdiff(fit$levels$perturbation, control)
  reporters <- if (is.null(reporters)) all_reporters else intersect(reporters, all_reporters)
  perturbations <- if (is.null(perturbations)) all_perturbations else intersect(perturbations, all_perturbations)

  grid <- expand.grid(
    reporter = reporters,
    perturbation = perturbations,
    stringsAsFactors = FALSE
  )
  if (nrow(grid) == 0) {
    return(data.frame())
  }
  destress_progress(
    verbose,
    "DStressR: extracting interaction-model contrasts for ",
    length(reporters),
    " reporters and ",
    length(perturbations),
    " perturbations."
  )

  base <- data.frame(
    .reporter = factor(grid$reporter, levels = fit$levels$reporter),
    .perturbation = factor(control, levels = fit$levels$perturbation)
  )
  comp <- base
  comp$.perturbation <- factor(grid$perturbation, levels = fit$levels$perturbation)

  full_model_frame <- fit_model_frame(fit$full_fit)
  for (col in fit$technical) {
    model_col <- full_model_frame[[col]]
    if (is.factor(model_col)) {
      level <- names(sort(table(model_col), decreasing = TRUE))[1]
      base[[col]] <- factor(level, levels = levels(model_col), ordered = is.ordered(model_col))
    } else if (is.numeric(model_col) || is.integer(model_col)) {
      base[[col]] <- stats::median(as.numeric(model_col), na.rm = TRUE)
    } else {
      level <- names(sort(table(model_col), decreasing = TRUE))[1]
      base[[col]] <- level
    }
    comp[[col]] <- base[[col]]
  }

  destress_progress(verbose, "DStressR: extracting additive total-effect contrasts.")
  total <- contrast_estimates(fit$total_fit, comp, base)

  destress_progress(verbose, "DStressR: extracting full interaction contrasts.")
  full_total <- contrast_estimates(fit$full_fit, comp, base)
  full_specific <- contrast_estimates(fit$full_fit, comp, base,
                                      center_by = grid$perturbation)

  destress_progress(verbose, "DStressR: estimating perturbation-level global effects.")
  global <- t(vapply(perturbations, function(cmp) {
    rows <- comp[grid$perturbation == cmp, , drop = FALSE]
    refs <- base[grid$perturbation == cmp, , drop = FALSE]
    vals <- contrast_estimates(fit$full_fit, rows, refs)
    colMeans(vals, na.rm = TRUE)
  }, numeric(4)))
  global_df <- data.frame(perturbation = perturbations, global_effect = global[, "estimate"])

  out <- data.frame(
    reporter = grid$reporter,
    perturbation = grid$perturbation,
    total_effect = full_total[, "estimate"],
    total_se = full_total[, "std_error"],
    total_statistic = full_total[, "statistic"],
    total_pvalue = full_total[, "pvalue"],
    specific_var = full_specific[, "std_error"]^2,
    additive_total_effect = total[, "estimate"],
    additive_total_se = total[, "std_error"],
    stringsAsFactors = FALSE
  )
  out$residual_df <- fit_df_residual(fit$full_fit)
  out$total_var_raw <- out$total_se^2
  total_df <- out$residual_df
  if (isTRUE(fit$empirical_bayes)) {
    moderated_total <- eb_moderate_se(out$total_se, total_df)
    out$total_se <- moderated_total$se
    total_df <- moderated_total$df
    total_test <- wald_t_test(out$total_effect, out$total_se, total_df)
    out$total_statistic <- total_test$statistic
    out$total_pvalue <- total_test$pvalue
  }
  out$total_df_for_test <- total_df
  out <- merge(out, global_df, by = "perturbation", sort = FALSE)
  destress_progress(verbose, "DStressR: applying low-rank background and p-value adjustment.")
  out$low_rank_effect <- low_rank_background_effect(
    out,
    effect = "total_effect",
    rank = fit$background_rank,
    impute = fit$background_impute
  )
  out$rank_adjusted_total_effect <- out$total_effect - out$low_rank_effect
  out$rank_adjusted_total_se <- out$total_se
  rank_adjusted_total_test <- wald_t_test(out$rank_adjusted_total_effect, out$rank_adjusted_total_se, out$total_df_for_test)
  out$rank_adjusted_total_statistic <- rank_adjusted_total_test$statistic
  out$rank_adjusted_total_pvalue <- rank_adjusted_total_test$pvalue
  rank_adjusted_global <- stats::aggregate(rank_adjusted_total_effect ~ perturbation, out, mean, na.rm = TRUE)
  names(rank_adjusted_global)[2] <- "rank_adjusted_global_effect"
  out <- merge(out, rank_adjusted_global, by = "perturbation", all.x = TRUE, sort = FALSE)
  out$specific_effect <- out$rank_adjusted_total_effect - out$rank_adjusted_global_effect
  out$specific_se <- sqrt(out$specific_var)
  out$specific_se[out$specific_var <= 0] <- NA_real_
  specific_df <- fit_df_residual(fit$full_fit)
  out$specific_df_raw <- rep(specific_df, nrow(out))
  if (isTRUE(fit$empirical_bayes)) {
    moderated <- eb_moderate_se(out$specific_se, specific_df)
    out$specific_se <- moderated$se
    specific_df <- moderated$df
  }
  out$specific_df_for_test <- rep_len(specific_df, nrow(out))
  out$specific_statistic <- out$specific_effect / out$specific_se
  out$specific_pvalue <- 2 * stats::pt(abs(out$specific_statistic),
                                       df = specific_df,
                                       lower.tail = FALSE)
  out$total_padj_global <- adjust_destress_pvalues(out$total_pvalue, out$reporter, "global")
  out$total_padj_by_reporter <- adjust_destress_pvalues(out$total_pvalue, out$reporter, "by_reporter")
  out$rank_adjusted_total_padj_global <- adjust_destress_pvalues(out$rank_adjusted_total_pvalue, out$reporter, "global")
  out$rank_adjusted_total_padj_by_reporter <- adjust_destress_pvalues(out$rank_adjusted_total_pvalue, out$reporter, "by_reporter")
  out$specific_padj_global <- adjust_destress_pvalues(out$specific_pvalue, out$reporter, "global")
  out$specific_padj_by_reporter <- adjust_destress_pvalues(out$specific_pvalue, out$reporter, "by_reporter")
  adjustment <- if (is.null(fit$adjustment)) "global" else fit$adjustment
  out$total_padj <- if (adjustment == "by_reporter") out$total_padj_by_reporter else if (adjustment == "none") out$total_pvalue else out$total_padj_global
  out$specific_padj <- if (adjustment == "by_reporter") out$specific_padj_by_reporter else if (adjustment == "none") out$specific_pvalue else out$specific_padj_global
  out <- out[, c(
    "reporter", "perturbation",
    "total_effect", "total_se", "total_statistic", "total_pvalue",
    "additive_total_effect", "additive_total_se",
    "total_var_raw", "residual_df", "total_df_for_test",
    "specific_var", "specific_df_raw", "specific_df_for_test",
    "global_effect", "low_rank_effect",
    "rank_adjusted_total_effect", "rank_adjusted_total_se",
    "rank_adjusted_total_statistic", "rank_adjusted_total_pvalue",
    "rank_adjusted_global_effect",
    "specific_effect", "specific_se", "specific_statistic", "specific_pvalue",
    "total_padj_global", "total_padj_by_reporter",
    "rank_adjusted_total_padj_global", "rank_adjusted_total_padj_by_reporter",
    "specific_padj_global",
    "specific_padj_by_reporter", "total_padj", "specific_padj"
  ), drop = FALSE]
  destress_progress(verbose, "DStressR: finished extracting model results.")
  out[order(out$reporter, out$perturbation), ]
}

#' Adjust p-values within reporter
#'
#' @param table Result table from [results()].
#' @param pvalue P-value column.
#' @param output Name of adjusted p-value column.
#' @param method Passed to [stats::p.adjust()].
#' @export
adjust_pvalues <- function(table, pvalue = "total_pvalue",
                           output = sub("_pvalue$", "_padj_by_reporter", pvalue),
                           method = "BH") {
  split_idx <- split(seq_len(nrow(table)), table$reporter)
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
#' @param effect Effect column. Defaults to `total_effect`; use
#'   `specific_effect` for deviations from the reporter-panel average.
#' @param padj Adjusted p-value column.
#' @export
call_hits <- function(table, fdr = 0.05, lfc = 0, effect = "total_effect",
                      padj = sub("_effect$", "_padj", effect)) {
  missing_cols <- setdiff(c(effect, padj), names(table))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
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
  if (!isTRUE(fit$interaction)) {
    return(data.frame(
      model = "reporter_lm",
      n_observations = sum(unique(fit$reporter_effects[c("reporter", "n_observations")])$n_observations),
      n_coefficients = sum(unique(fit$reporter_effects[c("reporter", "n_coefficients")])$n_coefficients),
      residual_df = sum(unique(fit$reporter_effects[c("reporter", "residual_df")])$residual_df),
      sigma = stats::median(unique(fit$reporter_effects[c("reporter", "sigma")])$sigma, na.rm = TRUE)
    ))
  }
  data.frame(
    model = c("additive", "interaction"),
    n_observations = c(stats::nobs(fit$total_fit), stats::nobs(fit$full_fit)),
    n_coefficients = c(length(stats::coef(fit$total_fit)), length(stats::coef(fit$full_fit))),
    residual_df = c(stats::df.residual(fit$total_fit), stats::df.residual(fit$full_fit)),
    sigma = c(stats::sigma(fit$total_fit), stats::sigma(fit$full_fit))
  )
}
