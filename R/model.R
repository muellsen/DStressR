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
    full = stats::as.formula(paste(".response ~ .promoter * .compound +", tech)),
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
  beta[is.na(beta)] <- 0
  estimate <- sum(contrast * beta)
  vc <- safe_vcov(fit)
  se <- sqrt(drop(t(contrast) %*% vc %*% contrast))
  df <- fit_df_residual(fit)
  statistic <- estimate / se
  pvalue <- 2 * stats::pt(abs(statistic), df = df, lower.tail = FALSE)
  c(estimate = estimate, std_error = se, statistic = statistic, pvalue = pvalue)
}

contrast_estimates <- function(fit, newdata_a, newdata_b) {
  terms_obj <- stats::delete.response(fit_terms(fit))
  x_a <- stats::model.matrix(terms_obj, newdata_a, contrasts.arg = fit_contrasts(fit))
  x_b <- stats::model.matrix(terms_obj, newdata_b, contrasts.arg = fit_contrasts(fit))
  contrast <- x_a - x_b
  beta <- fit_coef(fit)
  beta[is.na(beta)] <- 0
  estimate <- as.numeric(contrast %*% beta)
  vc <- safe_vcov(fit)
  se <- sqrt(rowSums((contrast %*% vc) * contrast))
  df <- fit_df_residual(fit)
  statistic <- estimate / se
  pvalue <- 2 * stats::pt(abs(statistic), df = df, lower.tail = FALSE)
  cbind(estimate = estimate, std_error = se, statistic = statistic, pvalue = pvalue)
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

low_rank_background_effect <- function(table, effect, rank) {
  rank <- validate_background_rank(rank)
  if (rank == 0 || nrow(table) == 0) {
    return(rep(0, nrow(table)))
  }

  required <- c("promoter", "compound", effect)
  missing <- setdiff(required, names(table))
  if (length(missing) > 0) {
    stop("Cannot estimate low-rank background; missing columns: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }

  promoters <- unique(as.character(table$promoter))
  compounds <- unique(as.character(table$compound))
  mat <- matrix(
    NA_real_,
    nrow = length(promoters),
    ncol = length(compounds),
    dimnames = list(promoters, compounds)
  )
  idx <- cbind(
    match(as.character(table$promoter), promoters),
    match(as.character(table$compound), compounds)
  )
  mat[idx] <- as.numeric(table[[effect]])

  observed <- is.finite(mat)
  col_mean <- colMeans(mat, na.rm = TRUE)
  col_mean[!is.finite(col_mean)] <- 0
  centered <- sweep(mat, 2, col_mean, "-")
  centered[!observed] <- 0

  rank <- min(rank, nrow(centered), ncol(centered))
  if (rank == 0 || all(abs(centered[observed]) < sqrt(.Machine$double.eps))) {
    return(rep(0, nrow(table)))
  }

  sv <- svd(centered, nu = rank, nv = rank)
  keep <- seq_len(rank)
  low_rank <- sv$u[, keep, drop = FALSE] %*%
    (diag(sv$d[keep], nrow = rank, ncol = rank) %*% t(sv$v[, keep, drop = FALSE]))

  # Keep the low-rank term orthogonal to the additive compound mean already
  # removed by global_effect, so residual specific effects still average to zero
  # within each compound over observed promoters.
  for (j in seq_len(ncol(low_rank))) {
    obs_j <- observed[, j]
    if (any(obs_j)) {
      low_rank[obs_j, j] <- low_rank[obs_j, j] - mean(low_rank[obs_j, j], na.rm = TRUE)
    }
  }
  low_rank[!observed] <- NA_real_
  as.numeric(low_rank[idx])
}

effect_matrix_from_table <- function(table, effect, promoter = "promoter", compound = "compound") {
  required <- c(promoter, compound, effect)
  missing <- setdiff(required, names(table))
  if (length(missing) > 0) {
    stop("Cannot build effect matrix; missing columns: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }

  promoters <- unique(as.character(table[[promoter]]))
  compounds <- unique(as.character(table[[compound]]))
  mat <- matrix(
    NA_real_,
    nrow = length(promoters),
    ncol = length(compounds),
    dimnames = list(promoters, compounds)
  )
  idx <- cbind(
    match(as.character(table[[promoter]]), promoters),
    match(as.character(table[[compound]]), compounds)
  )
  mat[idx] <- as.numeric(table[[effect]])
  mat
}

#' Diagnose low-rank background structure
#'
#' Computes singular values of a promoter-by-compound effect matrix and compares
#' them with a permutation null. The default permutation shuffles promoter
#' labels within each compound, preserving the compound-wise marginal
#' distribution while breaking shared promoter-loading structure.
#'
#' @param table Data frame with promoter, compound, and effect columns.
#' @param effect Numeric effect column to decompose, usually
#'   `specific_effect`.
#' @param promoter,compound Column names identifying promoters and compounds.
#' @param rank_max Maximum component index to report.
#' @param permutations Number of null permutations. Use `0` to skip the null.
#' @param seed Optional random seed for reproducible permutations.
#' @return A data frame with observed singular values, variance fractions, and
#'   optional permutation summaries.
#' @export
background_rank_diagnostics <- function(table,
                                        effect = "specific_effect",
                                        promoter = "promoter",
                                        compound = "compound",
                                        rank_max = 10,
                                        permutations = 100,
                                        seed = NULL) {
  rank_max <- validate_background_rank(rank_max)
  permutations <- validate_background_rank(permutations)
  mat <- effect_matrix_from_table(
    table,
    effect = effect,
    promoter = promoter,
    compound = compound
  )
  observed <- is.finite(mat)
  col_mean <- colMeans(mat, na.rm = TRUE)
  col_mean[!is.finite(col_mean)] <- 0
  centered <- sweep(mat, 2, col_mean, "-")
  centered[!observed] <- 0

  rank_max <- min(rank_max, nrow(centered), ncol(centered))
  if (rank_max == 0) {
    return(data.frame())
  }

  singular_values <- svd(centered, nu = 0, nv = 0)$d
  total_ss <- sum(singular_values^2)
  component <- seq_len(rank_max)
  observed_sv <- singular_values[component]
  prop_var <- if (total_ss > 0) observed_sv^2 / total_ss else rep(NA_real_, rank_max)

  null_median <- null_q95 <- null_q99 <- rep(NA_real_, rank_max)
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
      permuted <- centered
      for (j in seq_len(ncol(permuted))) {
        obs_j <- observed[, j]
        if (sum(obs_j) > 1) {
          permuted[obs_j, j] <- sample(permuted[obs_j, j])
        }
      }
      sv_b <- svd(permuted, nu = 0, nv = 0)$d
      null_sv[b, ] <- sv_b[component]
    }
    null_median <- apply(null_sv, 2, stats::median, na.rm = TRUE)
    null_q95 <- apply(null_sv, 2, stats::quantile, probs = 0.95, na.rm = TRUE)
    null_q99 <- apply(null_sv, 2, stats::quantile, probs = 0.99, na.rm = TRUE)
  }

  data.frame(
    component = component,
    observed = observed_sv,
    prop_variance = prop_var,
    cumulative_prop_variance = cumsum(prop_var),
    null_median = null_median,
    null_q95 = null_q95,
    null_q99 = null_q99,
    n_promoters = nrow(centered),
    n_compounds = ncol(centered),
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

fit_promoter_effects <- function(assay, technical, control) {
  technical_terms <- technical_formula(technical)
  technical_design_formula <- stats::as.formula(paste("~", technical_terms))
  by_promoter <- split(assay, assay$.promoter)

  rows <- lapply(names(by_promoter), function(promoter) {
    d <- by_promoter[[promoter]]
    d <- d[is.finite(d$.response), , drop = FALSE]
    if (nrow(d) == 0) {
      return(NULL)
    }

    compound <- as.character(d$.compound)
    compounds <- sort(setdiff(unique(compound), control))
    if (length(compounds) == 0 || !any(compound == control)) {
      return(NULL)
    }

    y <- d$.response
    z <- stats::model.matrix(technical_design_formula, data = d)
    z <- independent_columns(z)
    ztz_inv <- solve(crossprod(z))
    mz_y <- as.numeric(stats::lm.fit(z, y)$residuals)
    y_mz_y <- sum(mz_y^2)

    counts_all <- table(factor(compound, levels = c(control, compounds)))
    n_compound <- as.numeric(counts_all[compounds])
    names(n_compound) <- compounds

    z_sum_all <- rowsum(z, factor(compound, levels = c(control, compounds)), reorder = FALSE)
    b <- z_sum_all[compounds, , drop = FALSE]
    u_all <- rowsum(mz_y, factor(compound, levels = c(control, compounds)), reorder = FALSE)
    u <- as.numeric(u_all[compounds, , drop = TRUE])
    names(u) <- compounds

    n_inv <- 1 / n_compound
    s <- solve(solve(ztz_inv) - crossprod(b, b * n_inv))
    w <- n_inv * u
    beta <- w + n_inv * as.numeric(b %*% (s %*% crossprod(b, w)))
    a_inv_diag <- n_inv + n_inv^2 * rowSums((b %*% s) * b)

    sse <- y_mz_y - sum(beta * u)
    if (is.finite(sse) && sse < 0 && abs(sse) < sqrt(.Machine$double.eps)) {
      sse <- 0
    }
    df <- length(y) - ncol(z) - length(compounds)
    sigma2 <- sse / df
    se <- sqrt(sigma2 * a_inv_diag)
    test <- wald_t_test(beta, se, df)

    data.frame(
      promoter = promoter,
      compound = compounds,
      total_effect = beta,
      total_se = se,
      total_statistic = test$statistic,
      total_pvalue = test$pvalue,
      residual_df = df,
      sigma = sqrt(sigma2),
      n_observations = length(y),
      n_coefficients = ncol(z) + length(compounds),
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  if (is.null(out)) {
    out <- data.frame()
  }
  out
}

promoter_effect_results <- function(fit, compounds = NULL, promoters = NULL) {
  out <- fit$promoter_effects
  if (!is.null(compounds)) {
    out <- out[out$compound %in% compounds, , drop = FALSE]
  }
  if (nrow(out) == 0) {
    return(data.frame())
  }

  out$additive_total_effect <- NA_real_
  out$additive_total_se <- NA_real_
  df <- min(out$residual_df, na.rm = TRUE)
  if (isTRUE(fit$empirical_bayes)) {
    sigma <- stats::median(out$sigma, na.rm = TRUE)
    out$total_se <- eb_shrink(out$total_se, sigma, df)
    total_test <- wald_t_test(out$total_effect, out$total_se, df)
    out$total_statistic <- total_test$statistic
    out$total_pvalue <- total_test$pvalue
  }

  out$total_var <- out$total_se^2
  evc <- fit$empty_vector_promoter
  if (!is.null(evc) && nzchar(evc)) {
    evc_rows <- out[out$promoter == evc, c("compound", "total_effect", "total_var"), drop = FALSE]
    if (nrow(evc_rows) == 0) {
      stop("No fitted effects found for `empty_vector_promoter = \"", evc, "\"`.", call. = FALSE)
    }
    names(evc_rows) <- c("compound", "empty_vector_effect", "empty_vector_var")
    out <- merge(out, evc_rows, by = "compound", all.x = TRUE, sort = FALSE)
    out$background_adjusted_effect <- out$total_effect - out$empty_vector_effect
    out$background_adjusted_var <- out$total_var + out$empty_vector_var
    out <- out[out$promoter != evc, , drop = FALSE]
  } else {
    out$empty_vector_effect <- NA_real_
    out$empty_vector_var <- NA_real_
    out$background_adjusted_effect <- out$total_effect
    out$background_adjusted_var <- out$total_var
  }

  if (!is.null(promoters)) {
    out <- out[out$promoter %in% promoters, , drop = FALSE]
  }
  if (nrow(out) == 0) {
    return(data.frame())
  }

  global <- stats::aggregate(background_adjusted_effect ~ compound, out, mean, na.rm = TRUE)
  names(global)[2] <- "global_effect"
  out <- merge(out, global, by = "compound", all.x = TRUE, sort = FALSE)
  out$centering_effect <- out$global_effect
  out$low_rank_effect <- low_rank_background_effect(
    out,
    effect = "background_adjusted_effect",
    rank = fit$background_rank
  )
  out$specific_effect <- out$background_adjusted_effect - out$centering_effect - out$low_rank_effect

  variance_compounds <- unique(out$compound)
  total_variance_split <- split(out$total_var, out$compound)
  evc_variance_split <- split(out$empty_vector_var, out$compound)
  background_variance_split <- split(out$background_adjusted_var, out$compound)
  variance_summary <- data.frame(
    compound = variance_compounds,
    sum_total_var = vapply(total_variance_split[variance_compounds], sum, numeric(1), na.rm = TRUE),
    sum_background_adjusted_var = vapply(background_variance_split[variance_compounds], sum, numeric(1), na.rm = TRUE),
    empty_vector_var_for_compound = vapply(evc_variance_split[variance_compounds], function(x) {
      vals <- unique(x[is.finite(x)])
      if (length(vals) == 0) {
        0
      } else {
        vals[1]
      }
    }, numeric(1)),
    n_promoters_for_compound = vapply(total_variance_split[variance_compounds], function(x) sum(is.finite(x)), numeric(1)),
    stringsAsFactors = FALSE
  )
  out <- merge(out, variance_summary, by = "compound", all.x = TRUE, sort = FALSE)
  m <- out$n_promoters_for_compound
  out$global_se <- sqrt(out$sum_total_var / m^2 + out$empty_vector_var_for_compound)
  global_test <- wald_t_test(out$global_effect, out$global_se, df)
  out$global_statistic <- global_test$statistic
  out$global_pvalue <- global_test$pvalue
  out$specific_var <- ((m - 1) / m)^2 * out$total_var +
    (out$sum_total_var - out$total_var) / m^2
  out$specific_se <- sqrt(out$specific_var)
  out$specific_se[!is.finite(out$specific_se) | m <= 1] <- NA_real_
  if (isTRUE(fit$empirical_bayes)) {
    sigma <- stats::median(out$sigma, na.rm = TRUE)
    out$specific_se <- eb_shrink(out$specific_se, sigma, df)
  }
  specific_test <- wald_t_test(out$specific_effect, out$specific_se, df)
  out$specific_statistic <- specific_test$statistic
  out$specific_pvalue <- specific_test$pvalue

  out$total_padj_global <- adjust_destress_pvalues(out$total_pvalue, out$promoter, "global")
  out$total_padj_by_promoter <- adjust_destress_pvalues(out$total_pvalue, out$promoter, "by_promoter")
  out$specific_padj_global <- adjust_destress_pvalues(out$specific_pvalue, out$promoter, "global")
  out$specific_padj_by_promoter <- adjust_destress_pvalues(out$specific_pvalue, out$promoter, "by_promoter")
  adjustment <- if (is.null(fit$adjustment)) "global" else fit$adjustment
  out$total_padj <- if (adjustment == "by_promoter") out$total_padj_by_promoter else if (adjustment == "none") out$total_pvalue else out$total_padj_global
  out$specific_padj <- if (adjustment == "by_promoter") out$specific_padj_by_promoter else if (adjustment == "none") out$specific_pvalue else out$specific_padj_global

  out <- out[, c(
    "promoter", "compound",
    "total_effect", "total_se", "total_statistic", "total_pvalue",
    "additive_total_effect", "additive_total_se",
    "empty_vector_effect", "background_adjusted_effect",
    "global_effect", "global_se", "global_statistic", "global_pvalue",
    "low_rank_effect",
    "specific_effect", "specific_se", "specific_statistic", "specific_pvalue",
    "total_padj_global", "total_padj_by_promoter", "specific_padj_global",
    "specific_padj_by_promoter", "total_padj", "specific_padj"
  ), drop = FALSE]
  out[order(out$promoter, out$compound), ]
}

promoter_lm_results <- function(fit, compounds = NULL, promoters = NULL) {
  control <- fit$assay_info$control
  all_promoters <- fit$levels$promoter
  all_compounds <- setdiff(fit$levels$compound, control)
  promoters <- if (is.null(promoters)) all_promoters else intersect(promoters, all_promoters)
  compounds <- if (is.null(compounds)) all_compounds else intersect(compounds, all_compounds)

  if (length(promoters) == 0 || length(compounds) == 0) {
    return(data.frame())
  }

  rows <- lapply(promoters, function(promoter) {
    promoter_fit <- fit$promoter_fits[[promoter]]
    promoter_data <- fit_model_frame(promoter_fit)
    compound_levels <- levels(promoter_data$.compound)
    promoter_compounds <- intersect(compounds, setdiff(compound_levels, control))
    if (length(promoter_compounds) == 0) {
      return(NULL)
    }

    base <- data.frame(.compound = factor(control, levels = compound_levels))
    base <- base[rep(1, length(promoter_compounds)), , drop = FALSE]
    comp <- base
    comp$.compound <- factor(promoter_compounds, levels = compound_levels)
    for (col in fit$technical) {
      level <- modal_factor_level(promoter_data[[col]])
      base[[col]] <- factor(level, levels = levels(promoter_data[[col]]))
      comp[[col]] <- base[[col]]
    }

    total <- contrast_estimates(promoter_fit, comp, base)

    data.frame(
      promoter = promoter,
      compound = promoter_compounds,
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

  global <- stats::aggregate(total_effect ~ compound, out, mean, na.rm = TRUE)
  names(global)[2] <- "global_effect"
  out <- merge(out, global, by = "compound", all.x = TRUE, sort = FALSE)
  out$low_rank_effect <- low_rank_background_effect(
    out,
    effect = "total_effect",
    rank = fit$background_rank
  )
  out$specific_effect <- out$total_effect - out$global_effect - out$low_rank_effect

  df <- min(vapply(fit$promoter_fits[promoters], fit_df_residual, numeric(1)), na.rm = TRUE)
  if (isTRUE(fit$empirical_bayes)) {
    sigma <- stats::median(vapply(fit$promoter_fits[promoters], fit_sigma, numeric(1)), na.rm = TRUE)
    out$total_se <- eb_shrink(out$total_se, sigma, df)
    total_test <- wald_t_test(out$total_effect, out$total_se, df)
    out$total_statistic <- total_test$statistic
    out$total_pvalue <- total_test$pvalue
  }

  out$total_var <- out$total_se^2
  variance_summary <- stats::aggregate(
    total_var ~ compound,
    out,
    function(x) c(sum = sum(x, na.rm = TRUE), m = sum(is.finite(x)))
  )
  variance_summary <- do.call(data.frame, variance_summary)
  names(variance_summary) <- c("compound", "sum_total_var", "n_promoters_for_compound")
  out <- merge(out, variance_summary, by = "compound", all.x = TRUE, sort = FALSE)

  m <- out$n_promoters_for_compound
  out$global_se <- sqrt(out$sum_total_var) / m
  global_test <- wald_t_test(out$global_effect, out$global_se, df)
  out$global_statistic <- global_test$statistic
  out$global_pvalue <- global_test$pvalue
  out$specific_var <- ((m - 1) / m)^2 * out$total_var +
    (out$sum_total_var - out$total_var) / m^2
  out$specific_se <- sqrt(out$specific_var)
  out$specific_se[!is.finite(out$specific_se) | m <= 1] <- NA_real_

  if (isTRUE(fit$empirical_bayes)) {
    out$specific_se <- eb_shrink(out$specific_se, sigma, df)
  }
  specific_test <- wald_t_test(out$specific_effect, out$specific_se, df)
  out$specific_statistic <- specific_test$statistic
  out$specific_pvalue <- specific_test$pvalue

  out$total_padj_global <- adjust_destress_pvalues(out$total_pvalue, out$promoter, "global")
  out$total_padj_by_promoter <- adjust_destress_pvalues(out$total_pvalue, out$promoter, "by_promoter")
  out$specific_padj_global <- adjust_destress_pvalues(out$specific_pvalue, out$promoter, "global")
  out$specific_padj_by_promoter <- adjust_destress_pvalues(out$specific_pvalue, out$promoter, "by_promoter")
  adjustment <- if (is.null(fit$adjustment)) "global" else fit$adjustment
  out$total_padj <- if (adjustment == "by_promoter") out$total_padj_by_promoter else if (adjustment == "none") out$total_pvalue else out$total_padj_global
  out$specific_padj <- if (adjustment == "by_promoter") out$specific_padj_by_promoter else if (adjustment == "none") out$specific_pvalue else out$specific_padj_global

  out <- out[, c(
    "promoter", "compound",
    "total_effect", "total_se", "total_statistic", "total_pvalue",
    "additive_total_effect", "additive_total_se",
    "global_effect", "global_se", "global_statistic", "global_pvalue",
    "low_rank_effect",
    "specific_effect", "specific_se", "specific_statistic", "specific_pvalue",
    "total_padj_global", "total_padj_by_promoter", "specific_padj_global",
    "specific_padj_by_promoter", "total_padj", "specific_padj"
  ), drop = FALSE]
  out[order(out$promoter, out$compound), ]
}

observed_mean_results <- function(fit, compounds = NULL, promoters = NULL) {
  assay <- fit$assay_data
  control <- fit$assay_info$control
  all_promoters <- fit$levels$promoter
  promoters <- if (is.null(promoters)) all_promoters else intersect(promoters, all_promoters)

  d <- assay[assay$.promoter %in% promoters, , drop = FALSE]
  if (!is.null(compounds)) {
    d <- d[d$.compound %in% c(control, compounds), , drop = FALSE]
  }
  if (nrow(d) == 0) {
    return(data.frame())
  }

  d$.adjusted_response <- technical_adjusted_response(fit$total_fit, d, fit$technical)
  cell_mean <- stats::aggregate(
    .adjusted_response ~ .promoter + .compound,
    d,
    function(x) c(mean = mean(x, na.rm = TRUE), n = sum(is.finite(x)))
  )
  cell_mean <- do.call(data.frame, cell_mean)
  names(cell_mean) <- c("promoter", "compound", "mean_response", "n")
  cell_mean$promoter <- as.character(cell_mean$promoter)
  cell_mean$compound <- as.character(cell_mean$compound)
  cell_mean$n <- as.numeric(cell_mean$n)

  cell_key <- paste(cell_mean$promoter, cell_mean$compound, sep = "\r")
  mean_by_cell <- stats::setNames(cell_mean$mean_response, cell_key)
  d_key <- paste(as.character(d$.promoter), as.character(d$.compound), sep = "\r")
  within_residual <- d$.adjusted_response - mean_by_cell[d_key]
  residual_df <- sum(is.finite(within_residual)) - nrow(cell_mean)
  sigma <- sqrt(sum(within_residual^2, na.rm = TRUE) / residual_df)
  if (!is.finite(sigma) || residual_df <= 0) {
    sigma <- fit_sigma(fit$total_fit)
    residual_df <- fit_df_residual(fit$total_fit)
  }

  control_mean <- cell_mean[cell_mean$compound == control, c("promoter", "mean_response", "n"), drop = FALSE]
  names(control_mean) <- c("promoter", "control_mean_response", "control_n")

  out <- merge(
    cell_mean[cell_mean$compound != control, , drop = FALSE],
    control_mean,
    by = "promoter",
    all.x = TRUE,
    sort = FALSE
  )
  out <- out[is.finite(out$mean_response) & is.finite(out$control_mean_response), , drop = FALSE]
  if (!is.null(compounds)) {
    out <- out[out$compound %in% compounds, , drop = FALSE]
  }
  if (nrow(out) == 0) {
    return(data.frame())
  }

  out$total_effect <- out$mean_response - out$control_mean_response
  out$total_var <- sigma^2 * (1 / out$n + 1 / out$control_n)
  global <- stats::aggregate(total_effect ~ compound, out, mean, na.rm = TRUE)
  names(global)[2] <- "global_effect"
  out <- merge(out, global, by = "compound", all.x = TRUE, sort = FALSE)
  out$low_rank_effect <- low_rank_background_effect(
    out,
    effect = "total_effect",
    rank = fit$background_rank
  )
  out$specific_effect <- out$total_effect - out$global_effect - out$low_rank_effect

  variance_summary <- stats::aggregate(
    total_var ~ compound,
    out,
    function(x) c(sum = sum(x, na.rm = TRUE), m = sum(is.finite(x)))
  )
  variance_summary <- do.call(data.frame, variance_summary)
  names(variance_summary) <- c("compound", "sum_total_var", "n_promoters_for_compound")
  out <- merge(out, variance_summary, by = "compound", all.x = TRUE, sort = FALSE)

  out$total_se <- sqrt(out$total_var)
  m <- out$n_promoters_for_compound
  out$specific_var <- ((m - 1) / m)^2 * out$total_var +
    (out$sum_total_var - out$total_var) / m^2
  out$specific_se <- sqrt(out$specific_var)
  out$specific_se[!is.finite(out$specific_se) | m <= 1] <- NA_real_
  if (isTRUE(fit$empirical_bayes)) {
    out$specific_se <- eb_shrink(out$specific_se, sigma, residual_df)
  }
  out$specific_statistic <- out$specific_effect / out$specific_se
  out$specific_pvalue <- 2 * stats::pt(
    abs(out$specific_statistic),
    df = residual_df,
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

  out$total_padj_global <- adjust_destress_pvalues(out$total_pvalue, out$promoter, "global")
  out$total_padj_by_promoter <- adjust_destress_pvalues(out$total_pvalue, out$promoter, "by_promoter")
  out$specific_padj_global <- adjust_destress_pvalues(out$specific_pvalue, out$promoter, "global")
  out$specific_padj_by_promoter <- adjust_destress_pvalues(out$specific_pvalue, out$promoter, "by_promoter")
  adjustment <- if (is.null(fit$adjustment)) "global" else fit$adjustment
  out$total_padj <- if (adjustment == "by_promoter") out$total_padj_by_promoter else if (adjustment == "none") out$total_pvalue else out$total_padj_global
  out$specific_padj <- if (adjustment == "by_promoter") out$specific_padj_by_promoter else if (adjustment == "none") out$specific_pvalue else out$specific_padj_global

  out <- out[, c(
    "promoter", "compound",
    "total_effect", "total_se", "total_statistic", "total_pvalue",
    "additive_total_effect", "additive_total_se",
    "global_effect", "low_rank_effect",
    "specific_effect", "specific_se", "specific_statistic",
    "specific_pvalue",
    "total_padj_global", "total_padj_by_promoter", "specific_padj_global",
    "specific_padj_by_promoter", "total_padj", "specific_padj"
  ), drop = FALSE]
  out[order(out$promoter, out$compound), ]
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
#' @param empty_vector_promoter Optional promoter/control strain used as an
#'   empty-vector reporter in the model-based path. When supplied, its
#'   reference-relative compound effect is subtracted from every promoter's
#'   reference-relative compound effect before promoter-library centering.
#' @param background_rank Non-negative integer. The default `0` removes only
#'   the additive compound-wide mean. Values `1` or `2` additionally subtract a
#'   low-rank background term from the promoter-by-compound effect matrix before
#'   testing promoter-specific residual effects.
#' @param normalization One of `"linear_model"`, `"median_polish"`, or
#'   `"empty_vector"`. `"model"` and `"evc"` are accepted aliases.
#' @param testing One of `"student_t"`, `"moderated_t"`, or `"gaussian_z"`.
#' @param aggregation One of `"none"` or `"max_p"`.
#' @param adjustment One of `"global"`, `"by_promoter"`, or `"none"`.
#' @param interaction If `FALSE`, fit one Gaussian linear
#'   model per promoter with the control compound as reference and the supplied
#'   technical covariates as design terms. The latter is the scalable path for
#'   promoter-specific compound effects. If `TRUE`, fit the historical full
#'   promoter-by-compound interaction model.
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
                         empty_vector_promoter = NULL,
                         background_rank = 0,
                         normalization = NULL,
                         testing = NULL,
                         aggregation = NULL,
                         adjustment = NULL,
                         interaction = FALSE,
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
  background_rank <- validate_background_rank(background_rank)
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
  interaction <- isTRUE(interaction)
  if (!is.null(empty_vector_promoter)) {
    empty_vector_promoter <- as.character(empty_vector_promoter)
    if (length(empty_vector_promoter) != 1 || is.na(empty_vector_promoter) || !nzchar(empty_vector_promoter)) {
      stop("`empty_vector_promoter` must be one promoter label.", call. = FALSE)
    }
    if (!empty_vector_promoter %in% levels(assay$.promoter)) {
      stop("Empty-vector promoter '", empty_vector_promoter, "' was not found in the assay.", call. = FALSE)
    }
  }
  formulas <- make_formulas(technical)
  total_fit <- if (interaction) {
    stats::lm(formulas$total, data = assay, na.action = stats::na.exclude)
  } else {
    NULL
  }
  full_fit <- if (interaction) {
    stats::lm(formulas$full, data = assay, na.action = stats::na.exclude)
  } else {
    NULL
  }
  assay_data <- if (interaction) {
    NULL
  } else {
    assay
  }
  promoter_formula <- stats::as.formula(paste(".response ~ .compound +", technical_formula(technical)))
  promoter_fits <- if (interaction) {
    NULL
  } else {
    NULL
  }
  promoter_effects <- if (interaction) {
    NULL
  } else {
    fit_promoter_effects(assay, technical, attr(assay, "destress")$control)
  }

  structure(
    list(
      total_fit = total_fit,
      full_fit = full_fit,
      interaction = interaction,
      assay_data = assay_data,
      promoter_fits = promoter_fits,
      promoter_effects = promoter_effects,
      growth_exponents = attr(assay, "destress")$growth_exponent_fit,
      assay_info = attr(assay, "destress"),
      levels = list(
        promoter = levels(assay$.promoter),
        compound = levels(assay$.compound)
      ),
      technical = technical,
      empirical_bayes = empirical_bayes,
      empty_vector_promoter = empty_vector_promoter,
      background_rank = background_rank,
      stages = stages,
      preset = if (is.null(preset)) "model" else preset,
      adjustment = stages$adjustment
    ),
    class = "destress_fit"
  )
}

#' Extract estimated model parameters
#'
#' @param fit A `destress_fit` object.
#' @return A named list of estimated parameter tables available for the fitted
#'   model. The scalable model path includes promoter-specific growth
#'   normalization estimates and promoter-compound effect estimates.
#' @export
model_parameters <- function(fit) {
  if (!inherits(fit, "destress_fit")) {
    stop("`fit` must be a destress_fit.", call. = FALSE)
  }

  out <- list(
    background = data.frame(
      background_rank = validate_background_rank(fit$background_rank)
    ),
    growth_exponents = fit$growth_exponents,
    promoter_effects = fit$promoter_effects
  )

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
  if (!isTRUE(fit$interaction)) {
    return(promoter_effect_results(fit, compounds = compounds, promoters = promoters))
  }
  if (!is.null(fit$empty_vector_promoter)) {
    stop("Model-based empty-vector adjustment is currently implemented for the scalable promoter-specific path.", call. = FALSE)
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
  out$low_rank_effect <- low_rank_background_effect(
    out,
    effect = "total_effect",
    rank = fit$background_rank
  )
  out$specific_effect <- out$total_effect - out$global_effect - out$low_rank_effect
  out$specific_se <- out$total_se
  if (isTRUE(fit$empirical_bayes)) {
    out$specific_se <- eb_shrink(out$specific_se, fit_sigma(fit$full_fit), fit_df_residual(fit$full_fit))
  }
  out$specific_statistic <- out$specific_effect / out$specific_se
  out$specific_pvalue <- 2 * stats::pt(abs(out$specific_statistic),
                                       df = fit_df_residual(fit$full_fit),
                                       lower.tail = FALSE)
  out$total_padj_global <- adjust_destress_pvalues(out$total_pvalue, out$promoter, "global")
  out$total_padj_by_promoter <- adjust_destress_pvalues(out$total_pvalue, out$promoter, "by_promoter")
  out$specific_padj_global <- adjust_destress_pvalues(out$specific_pvalue, out$promoter, "global")
  out$specific_padj_by_promoter <- adjust_destress_pvalues(out$specific_pvalue, out$promoter, "by_promoter")
  adjustment <- if (is.null(fit$adjustment)) "global" else fit$adjustment
  out$total_padj <- if (adjustment == "by_promoter") out$total_padj_by_promoter else if (adjustment == "none") out$total_pvalue else out$total_padj_global
  out$specific_padj <- if (adjustment == "by_promoter") out$specific_padj_by_promoter else if (adjustment == "none") out$specific_pvalue else out$specific_padj_global
  out <- out[, c(
    "promoter", "compound",
    "total_effect", "total_se", "total_statistic", "total_pvalue",
    "additive_total_effect", "additive_total_se",
    "global_effect", "low_rank_effect",
    "specific_effect", "specific_se", "specific_statistic", "specific_pvalue",
    "total_padj_global", "total_padj_by_promoter", "specific_padj_global",
    "specific_padj_by_promoter", "total_padj", "specific_padj"
  ), drop = FALSE]
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
  if (!isTRUE(fit$interaction)) {
    return(data.frame(
      model = "promoter_glm",
      n_observations = sum(unique(fit$promoter_effects[c("promoter", "n_observations")])$n_observations),
      n_coefficients = sum(unique(fit$promoter_effects[c("promoter", "n_coefficients")])$n_coefficients),
      residual_df = sum(unique(fit$promoter_effects[c("promoter", "residual_df")])$residual_df),
      sigma = stats::median(unique(fit$promoter_effects[c("promoter", "sigma")])$sigma, na.rm = TRUE)
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
