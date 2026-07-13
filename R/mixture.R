scaled_t_density <- function(x, location, scale, df) {
  stats::dt((x - location) / scale, df = df) / scale
}

fit_three_part_t <- function(x,
                             df = 4,
                             max_iter = 500,
                             tol = 1e-7,
                             min_scale = 1e-4,
                             min_prior = 1e-4) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < 10) {
    stop("Need at least 10 finite effects to fit the three-part mixture.", call. = FALSE)
  }
  global_scale <- stats::mad(x, constant = 1.4826, na.rm = TRUE)
  if (!is.finite(global_scale) || global_scale <= 0) {
    global_scale <- stats::sd(x, na.rm = TRUE)
  }
  if (!is.finite(global_scale) || global_scale <= 0) {
    global_scale <- 1
  }

  mu <- as.numeric(stats::quantile(x, c(0.15, 0.5, 0.85), na.rm = TRUE, names = FALSE))
  if (length(unique(mu)) < 3) {
    mu <- stats::median(x, na.rm = TRUE) + c(-0.5, 0, 0.5) * global_scale
  }
  scale <- rep(max(global_scale, min_scale), 3)
  prior <- c(0.15, 0.70, 0.15)
  loglik_old <- -Inf

  for (iter in seq_len(max_iter)) {
    dens <- vapply(seq_len(3), function(k) {
      prior[k] * scaled_t_density(x, mu[k], scale[k], df)
    }, numeric(n))
    denom <- rowSums(dens)
    denom[!is.finite(denom) | denom <= 0] <- .Machine$double.xmin
    z <- dens / denom

    latent_precision <- vapply(seq_len(3), function(k) {
      residual2 <- ((x - mu[k]) / scale[k])^2
      (df + 1) / (df + residual2)
    }, numeric(n))

    nk <- colSums(z)
    prior <- pmax(nk / n, min_prior)
    prior <- prior / sum(prior)
    for (k in seq_len(3)) {
      w <- z[, k] * latent_precision[, k]
      if (sum(w) > 0) {
        mu[k] <- sum(w * x) / sum(w)
        scale[k] <- sqrt(sum(w * (x - mu[k])^2) / max(nk[k], .Machine$double.eps))
        scale[k] <- max(scale[k], min_scale)
      }
    }

    ord <- order(mu)
    mu <- mu[ord]
    scale <- scale[ord]
    prior <- prior[ord]

    dens <- vapply(seq_len(3), function(k) {
      prior[k] * scaled_t_density(x, mu[k], scale[k], df)
    }, numeric(n))
    loglik <- sum(log(pmax(rowSums(dens), .Machine$double.xmin)))
    if (is.finite(loglik_old) && abs(loglik - loglik_old) < tol * (abs(loglik_old) + tol)) {
      break
    }
    loglik_old <- loglik
  }

  list(
    prior = prior,
    location = mu,
    scale = scale,
    df = df,
    logLik = loglik_old,
    iterations = iter,
    converged = iter < max_iter
  )
}

predict_three_part_t <- function(x, fit) {
  dens <- vapply(seq_len(3), function(k) {
    fit$prior[k] * scaled_t_density(x, fit$location[k], fit$scale[k], fit$df)
  }, numeric(length(x)))
  denom <- rowSums(dens)
  denom[!is.finite(denom) | denom <= 0] <- .Machine$double.xmin
  post <- dens / denom
  colnames(post) <- c("prob_repressed", "prob_null", "prob_activated")
  null_z <- (x - fit$location[2]) / fit$scale[2]
  null_p <- 2 * stats::pt(abs(null_z), df = fit$df, lower.tail = FALSE)
  class_id <- max.col(post, ties.method = "first")
  posterior_class <- c("repressed", "null", "activated")[class_id]
  data.frame(
    prob_repressed = post[, 1],
    prob_null = post[, 2],
    prob_activated = post[, 3],
    local_fdr = post[, 2],
    empirical_null_z = null_z,
    empirical_null_pvalue = null_p,
    posterior_class = posterior_class,
    stringsAsFactors = FALSE
  )
}

#' Fit a three-part empirical-null mixture to promoter-compound effects
#'
#' Fits, separately for each promoter, a three-component Student-t mixture to
#' adjusted promoter-compound effects. The ordered components are interpreted as
#' repressed, null, and activated effects. This second-stage model is intended
#' for empirical-null calibration after the first-stage DStressR model has
#' already adjusted growth, technical factors, compound-wide effects, and
#' promoter-specific variance.
#'
#' @param table A data frame with one row per promoter-compound pair.
#' @param value Numeric effect column, usually `specific_effect` or a centered
#'   DStressR EB effect column.
#' @param promoter Column identifying promoters.
#' @param df Degrees of freedom for each Student-t component. Smaller values
#'   give heavier tails.
#' @param max_iter Maximum EM iterations per promoter.
#' @param tol Relative log-likelihood convergence tolerance.
#' @param min_scale Lower bound for component scale.
#' @param min_prior Lower bound for component mixing proportions.
#' @param padj_method Multiple-testing correction method passed to
#'   [stats::p.adjust()], applied within promoter to empirical-null p-values.
#' @return The input table with posterior probabilities, local FDR,
#'   empirical-null p-values, promoter-wise adjusted p-values, and posterior
#'   class appended. The promoter-level fitted parameters are available as
#'   `attr(result, "mixture_summary")`.
#' @export
fit_effect_mixture <- function(table,
                               value = "specific_effect",
                               promoter = "promoter",
                               df = 4,
                               max_iter = 2000,
                               tol = 1e-6,
                               min_scale = 1e-4,
                               min_prior = 1e-4,
                               padj_method = "BH") {
  stopifnot(is.data.frame(table))
  required <- c(value, promoter)
  missing_cols <- setdiff(required, names(table))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  out <- table
  out$.effect_mixture_value <- as.numeric(out[[value]])
  out$.effect_mixture_promoter <- as.character(out[[promoter]])
  out$prob_repressed <- NA_real_
  out$prob_null <- NA_real_
  out$prob_activated <- NA_real_
  out$local_fdr <- NA_real_
  out$posterior_nonnull <- NA_real_
  out$local_fdr_qvalue_by_promoter <- NA_real_
  out$empirical_null_z <- NA_real_
  out$empirical_null_pvalue <- NA_real_
  out$empirical_null_padj_by_promoter <- NA_real_
  out$posterior_class <- NA_character_

  split_idx <- split(seq_len(nrow(out)), out$.effect_mixture_promoter)
  summary_rows <- list()
  for (prom in names(split_idx)) {
    idx <- split_idx[[prom]]
    finite_idx <- idx[is.finite(out$.effect_mixture_value[idx])]
    if (length(finite_idx) < 10) {
      next
    }
    x <- out$.effect_mixture_value[finite_idx]
    fit <- fit_three_part_t(
      x,
      df = df,
      max_iter = max_iter,
      tol = tol,
      min_scale = min_scale,
      min_prior = min_prior
    )
    pred <- predict_three_part_t(x, fit)
    out[finite_idx, c(
      "prob_repressed",
      "prob_null",
      "prob_activated",
      "local_fdr",
      "empirical_null_z",
      "empirical_null_pvalue",
      "posterior_class"
    )] <- pred
    out$posterior_nonnull[finite_idx] <- 1 - out$local_fdr[finite_idx]
    out$empirical_null_padj_by_promoter[finite_idx] <- stats::p.adjust(
      out$empirical_null_pvalue[finite_idx],
      method = padj_method
    )
    lfdr_order <- finite_idx[order(out$local_fdr[finite_idx])]
    cumulative_fdr <- cumsum(out$local_fdr[lfdr_order]) / seq_along(lfdr_order)
    local_q <- rev(cummin(rev(cumulative_fdr)))
    out$local_fdr_qvalue_by_promoter[lfdr_order] <- local_q
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      promoter = prom,
      n = length(finite_idx),
      df = fit$df,
      logLik = fit$logLik,
      iterations = fit$iterations,
      converged = fit$converged,
      prior_repressed = fit$prior[1],
      prior_null = fit$prior[2],
      prior_activated = fit$prior[3],
      location_repressed = fit$location[1],
      location_null = fit$location[2],
      location_activated = fit$location[3],
      scale_repressed = fit$scale[1],
      scale_null = fit$scale[2],
      scale_activated = fit$scale[3],
      stringsAsFactors = FALSE
    )
  }

  out$.effect_mixture_value <- NULL
  out$.effect_mixture_promoter <- NULL
  summary <- if (length(summary_rows) > 0) {
    do.call(rbind, summary_rows)
  } else {
    data.frame()
  }
  attr(out, "mixture_summary") <- summary
  out
}
