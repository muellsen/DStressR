# Avoid R CMD check notes for data-frame columns used inside ggplot2 aesthetics.
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    ".color_effect", ".density", ".effect_variance", ".log10_variance", ".rank",
    ".signed_mean_effect", "curve", "density", "variance_trend", "x", "xend",
    "y", "yend"
  ))
}

#' Perturbation-level response diagnostics
#'
#' Summarizes a reporter-by-perturbation effect table into perturbation-level
#' diagnostics. The primary use is to compare the absolute mean response of a
#' perturbation with the variance of its reporter-specific responses, analogous
#' in spirit to mean-variance diagnostic plots used for count-data workflows.
#'
#' The function is deliberately generic: columns are referred to as reporters
#' and perturbations, although DStressR result tables usually contain promoters
#' and compounds. Large absolute mean effects with low variance indicate
#' coherent, broad responses, whereas unusually large variance conditional on
#' the absolute mean effect indicates heterogeneous reporter behavior.
#'
#' @param table A data frame with one row per reporter-perturbation pair.
#' @param mean_effect Numeric column used to compute the perturbation-level mean
#'   effect. For DStressR results, `total_effect` or
#'   `rank_adjusted_total_effect` are typical choices.
#' @param variance_effect Numeric column whose cross-reporter variance is
#'   computed for each perturbation. Defaults to `mean_effect`.
#' @param reporter,perturbation Columns identifying reporters and perturbations.
#' @param perturbation_label Optional human-readable perturbation label column.
#'   Defaults to `perturbation`.
#' @param min_reporters Minimum number of finite reporter-level effects required
#'   for a perturbation.
#' @param trend_span Span used for the loess trend of variance over
#'   `rank(abs(mean effect))`.
#' @return A data frame with one row per perturbation.
#' @export
perturbation_diagnostics <- function(table,
                                     mean_effect = "total_effect",
                                     variance_effect = mean_effect,
                                     reporter = "promoter",
                                     perturbation = "compound",
                                     perturbation_label = perturbation,
                                     min_reporters = 2,
                                     trend_span = 0.45) {
  stopifnot(is.data.frame(table))
  required <- c(mean_effect, variance_effect, reporter, perturbation, perturbation_label)
  missing_cols <- setdiff(required, names(table))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  if (!is.numeric(min_reporters) || length(min_reporters) != 1 || min_reporters < 2) {
    stop("`min_reporters` must be a number of at least 2.", call. = FALSE)
  }

  d <- table[, required, drop = FALSE]
  names(d) <- c(".mean_effect", ".variance_effect", ".reporter", ".perturbation", ".perturbation_label")
  d$.mean_effect <- as.numeric(d$.mean_effect)
  d$.variance_effect <- as.numeric(d$.variance_effect)
  d$.reporter <- as.character(d$.reporter)
  d$.perturbation <- as.character(d$.perturbation)
  d$.perturbation_label <- as.character(d$.perturbation_label)
  missing_label <- is.na(d$.perturbation_label) | !nzchar(d$.perturbation_label)
  d$.perturbation_label[missing_label] <- d$.perturbation[missing_label]

  split_idx <- split(seq_len(nrow(d)), d$.perturbation)
  out <- do.call(rbind, lapply(split_idx, function(idx) {
    sub <- d[idx, , drop = FALSE]
    ok_mean <- is.finite(sub$.mean_effect)
    ok_var <- is.finite(sub$.variance_effect)
    n_mean <- sum(ok_mean)
    n_variance <- sum(ok_var)
    if (n_mean < min_reporters || n_variance < min_reporters) {
      return(NULL)
    }
    data.frame(
      perturbation = sub$.perturbation[1],
      perturbation_label = sub$.perturbation_label[1],
      n_reporters = length(unique(sub$.reporter[ok_mean | ok_var])),
      n_mean_effect = n_mean,
      n_variance_effect = n_variance,
      mean_effect = mean(sub$.mean_effect[ok_mean]),
      abs_mean_effect = abs(mean(sub$.mean_effect[ok_mean])),
      effect_variance = stats::var(sub$.variance_effect[ok_var]),
      stringsAsFactors = FALSE
    )
  }))
  if (is.null(out) || nrow(out) == 0) {
    stop("No perturbations have enough finite reporter effects.", call. = FALSE)
  }

  out <- out[is.finite(out$effect_variance), , drop = FALSE]
  out <- out[order(out$abs_mean_effect, out$perturbation_label, out$perturbation), , drop = FALSE]
  out$rank_abs_mean_effect <- seq_len(nrow(out))
  out$rank_effect_variance <- rank(-out$effect_variance, ties.method = "min")

  out$variance_trend <- NA_real_
  if (nrow(out) >= 10 && length(unique(out$rank_abs_mean_effect)) >= 10) {
    trend_fit <- try(suppressWarnings(
      stats::loess(
        effect_variance ~ rank_abs_mean_effect,
        data = out,
        span = trend_span,
        degree = 1,
        control = stats::loess.control(surface = "direct")
      )),
      silent = TRUE
    )
    if (!inherits(trend_fit, "try-error")) {
      out$variance_trend <- as.numeric(suppressWarnings(stats::predict(trend_fit, newdata = out)))
    }
  }
  out$variance_residual <- out$effect_variance - out$variance_trend
  out$rank_variance_residual <- rank(-out$variance_residual, ties.method = "min", na.last = "keep")
  out$rank_abs_mean_effect_desc <- rank(-out$abs_mean_effect, ties.method = "min")

  rownames(out) <- NULL
  out
}

#' Plot perturbation-level mean-variance diagnostics
#'
#' Draws a DESeq-style diagnostic plot in which perturbations are ordered by
#' their absolute mean response and the y-axis shows the variance of
#' reporter-level effects. A non-parametric trend can be overlaid, and the most
#' heterogeneous perturbations can be labelled.
#'
#' @inheritParams perturbation_diagnostics
#' @param diagnostics Optional output from [perturbation_diagnostics()]. If
#'   supplied, `table` is ignored.
#' @param add_trend If `TRUE`, overlay the loess trend.
#' @param label_by Which perturbations to label. `"residual"` labels
#'   perturbations with the largest variance residual above the trend,
#'   `"variance"` labels the largest variances, `"abs_mean"` labels the largest
#'   absolute mean effects, and `"none"` suppresses labels.
#' @param top_n Number of perturbations to label.
#' @param point_size Point size.
#' @param color_limits Optional numeric vector of length two for the color
#'   scale. If omitted, symmetric robust limits are computed from
#'   `color_quantile`.
#' @param color_quantile Quantile of the absolute mean effect used for robust
#'   color limits when `color_limits` is omitted.
#' @param title,subtitle,xlab,ylab Plot labels.
#' @param legend_position,legend_justification Numeric vectors passed to
#'   `ggplot2::theme()` to place the legend inside the plotting region.
#' @return A `ggplot` object. The diagnostic table is available as
#'   `attr(plot, "diagnostics")`.
#' @export
plot_mean_variance_diagnostic <- function(table = NULL,
                                          diagnostics = NULL,
                                          mean_effect = "total_effect",
                                          variance_effect = mean_effect,
                                          reporter = "promoter",
                                          perturbation = "compound",
                                          perturbation_label = perturbation,
                                          min_reporters = 2,
                                          trend_span = 0.45,
                                          add_trend = TRUE,
                                          label_by = c("residual", "variance", "abs_mean", "none"),
                                          top_n = 8,
                                          point_size = 1.8,
                                          color_limits = NULL,
                                          color_quantile = 0.92,
                                          title = NULL,
                                          subtitle = NULL,
                                          xlab = NULL,
                                          ylab = NULL,
                                          legend_position = c(0.035, 0.965),
                                          legend_justification = c(0, 1)) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package `ggplot2` is required for plot_mean_variance_diagnostic().", call. = FALSE)
  }
  label_by <- match.arg(label_by)
  if (is.null(diagnostics)) {
    if (is.null(table)) {
      stop("Provide either `table` or `diagnostics`.", call. = FALSE)
    }
    diagnostics <- perturbation_diagnostics(
      table,
      mean_effect = mean_effect,
      variance_effect = variance_effect,
      reporter = reporter,
      perturbation = perturbation,
      perturbation_label = perturbation_label,
      min_reporters = min_reporters,
      trend_span = trend_span
    )
  }
  stopifnot(is.data.frame(diagnostics))
  required <- c("rank_abs_mean_effect", "effect_variance", "mean_effect", "abs_mean_effect", "perturbation_label")
  missing_cols <- setdiff(required, names(diagnostics))
  if (length(missing_cols) > 0) {
    stop("`diagnostics` is missing columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  d <- diagnostics
  d$.signed_mean_effect <- as.numeric(d$mean_effect)
  d$.effect_variance <- as.numeric(d$effect_variance)
  d$.rank <- as.numeric(d$rank_abs_mean_effect)
  d <- d[is.finite(d$.rank) & is.finite(d$.effect_variance), , drop = FALSE]
  if (nrow(d) == 0) {
    stop("No finite diagnostic rows available for plotting.", call. = FALSE)
  }

  if (is.null(color_limits)) {
    finite_color <- abs(d$.signed_mean_effect[is.finite(d$.signed_mean_effect)])
    color_max <- if (length(finite_color) > 0) {
      stats::quantile(finite_color, probs = color_quantile, names = FALSE, na.rm = TRUE)
    } else {
      NA_real_
    }
    if (!is.finite(color_max) || color_max <= 0) {
      color_max <- max(finite_color, na.rm = TRUE)
    }
    if (!is.finite(color_max) || color_max <= 0) {
      color_max <- 1
    }
    color_limits <- c(-color_max, color_max)
  } else {
    color_limits <- as.numeric(color_limits)
    if (length(color_limits) != 2 || any(!is.finite(color_limits)) || color_limits[1] >= color_limits[2]) {
      stop("`color_limits` must be a finite increasing numeric vector of length two.", call. = FALSE)
    }
  }
  d$.color_effect <- pmax(pmin(d$.signed_mean_effect, color_limits[2]), color_limits[1])

  label_df <- d[FALSE, , drop = FALSE]
  if (label_by != "none" && top_n > 0) {
    order_col <- switch(
      label_by,
      residual = if ("variance_residual" %in% names(d)) d$variance_residual else d$effect_variance,
      variance = d$effect_variance,
      abs_mean = d$abs_mean_effect
    )
    label_df <- d[order(-order_col), , drop = FALSE]
    label_df <- label_df[is.finite(order_col[order(-order_col)]), , drop = FALSE]
    label_df <- utils::head(label_df, top_n)
  }

  if (is.null(subtitle)) {
    subtitle <- paste0(
      "Mean: ", mean_effect,
      "; variance: ", variance_effect
    )
  }
  if (is.null(xlab)) {
    xlab <- destress_effect_axis_label(mean_effect, type = "rank_abs_mean")
  }
  if (is.null(ylab)) {
    ylab <- destress_effect_axis_label(variance_effect, type = "variance")
  }
  color_label <- destress_effect_axis_label(mean_effect, type = "mean")

  p <- ggplot2::ggplot(d, ggplot2::aes(.rank, .effect_variance)) +
    ggplot2::geom_point(
      ggplot2::aes(color = .color_effect),
      size = point_size,
      alpha = 0.82
    ) +
    ggplot2::scale_color_gradient2(
      low = "#1D4ED8",
      mid = "#9CA3AF",
      high = "#B91C1C",
      midpoint = 0,
      limits = color_limits,
      name = color_label
    ) +
    destress_diagnostic_theme(base_size = 10, legend_position = legend_position,
                              legend_justification = legend_justification) +
    ggplot2::guides(
      color = ggplot2::guide_colorbar(
        direction = "vertical",
        title.position = "top",
        title.hjust = 0.5,
        barwidth = grid::unit(0.08, "in"),
        barheight = grid::unit(0.58, "in")
      )
    ) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = xlab,
      y = ylab
    )

  if (isTRUE(add_trend) && "variance_trend" %in% names(d) && any(is.finite(d$variance_trend))) {
    trend_df <- d[is.finite(d$variance_trend), , drop = FALSE]
    p <- p +
      ggplot2::geom_line(
        data = trend_df,
        ggplot2::aes(.rank, variance_trend),
        inherit.aes = FALSE,
        color = "#008A5B",
        linewidth = 0.65
      )
  }

  if (nrow(label_df) > 0) {
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p +
        ggrepel::geom_text_repel(
          data = label_df,
          ggplot2::aes(label = perturbation_label),
          size = 2.6,
          min.segment.length = 0,
          box.padding = 0.25,
          point.padding = 0.12,
          max.overlaps = Inf,
          seed = 1,
          show.legend = FALSE
        )
    } else {
      p <- p +
        ggplot2::geom_text(
          data = label_df,
          ggplot2::aes(label = perturbation_label),
          size = 2.6,
          vjust = -0.65,
          check_overlap = TRUE,
          show.legend = FALSE
        )
    }
  }

  attr(p, "diagnostics") <- diagnostics
  p
}

#' Fit diagnostic distributions to perturbation-level variances
#'
#' Fits simple positive-support or log-scale distributions to the variance
#' column produced by [perturbation_diagnostics()]. This is intended for
#' diagnostic use: it can help identify broad distributional structure in
#' perturbation-level heterogeneity, but it does not change DStressR inference.
#'
#' The beta-prime and inverse-gamma distributions are fitted on the raw
#' variance scale and evaluated on the log10 scale. The log-normal and log-t
#' distributions are fitted directly to log10 variances. A two-Gaussian mixture
#' on log10 variances is available as an explicitly requested exploratory
#' option, but is not used by default.
#'
#' @param diagnostics A data frame, usually from [perturbation_diagnostics()].
#' @param variance Numeric column containing positive variance estimates.
#' @param distributions Character vector of distributions to fit. Supported
#'   values are `"beta_prime"`, `"inverse_gamma"`, `"log_normal"`, `"log_t"`,
#'   and `"two_gaussian"`.
#' @param min_n Minimum number of finite positive variances required.
#' @param gmm_starts Number of random starts for the two-Gaussian mixture.
#' @param gmm_min_weight Minimum component weight used during mixture fitting.
#' @param gmm_sigma_floor Optional lower bound for component standard deviations
#'   on the log10 scale. If `NULL`, a data-adaptive floor is used.
#' @param seed Random seed used for mixture starting values.
#' @return A data frame with one row per fitted distribution. The fitted
#'   variance column is stored as the `"variance"` attribute.
#' @export
fit_variance_distribution <- function(diagnostics,
                                      variance = "effect_variance",
                                      distributions = c("beta_prime", "inverse_gamma", "log_normal", "log_t"),
                                      min_n = 10,
                                      gmm_starts = 50,
                                      gmm_min_weight = 0.08,
                                      gmm_sigma_floor = NULL,
                                      seed = 1) {
  stopifnot(is.data.frame(diagnostics))
  if (!variance %in% names(diagnostics)) {
    stop("`diagnostics` is missing variance column `", variance, "`.", call. = FALSE)
  }
  distributions[distributions == "scaled_inverse_chisq"] <- "inverse_gamma"
  supported_distributions <- c("beta_prime", "inverse_gamma", "log_normal", "log_t", "two_gaussian")
  distributions <- match.arg(distributions, supported_distributions, several.ok = TRUE)
  x <- as.numeric(diagnostics[[variance]])
  x <- x[is.finite(x) & x > 0]
  if (length(x) < min_n) {
    stop("At least ", min_n, " finite positive variances are required.", call. = FALSE)
  }
  z <- log10(x)

  rows <- list()
  if ("beta_prime" %in% distributions) {
    fit <- fit_beta_prime_variance(x)
    if (!is.null(fit)) {
      log_density <- beta_prime_log10_density(
        z,
        shape1 = fit$shape1,
        shape2 = fit$shape2,
        scale = fit$scale,
        log = TRUE
      )
      if (all(is.finite(log_density))) {
        loglik <- sum(log_density)
        rows[[length(rows) + 1]] <- variance_distribution_row(
          distribution = "beta-prime",
          n = length(x),
          logLik = loglik,
          n_parameters = 3,
          shape1 = fit$shape1,
          shape2 = fit$shape2,
          scale = fit$scale,
          raw_logLik = fit$logLik
        )
      }
    }
  }
  if ("inverse_gamma" %in% distributions) {
    fit <- fit_inverse_gamma_variance(x)
    if (!is.null(fit)) {
      log_density <- inverse_gamma_log10_density(
        z,
        shape = fit$shape,
        scale = fit$scale,
        log = TRUE
      )
      if (all(is.finite(log_density))) {
        rows[[length(rows) + 1]] <- variance_distribution_row(
          distribution = "inverse-gamma",
          n = length(x),
          logLik = sum(log_density),
          n_parameters = 2,
          raw_logLik = fit$logLik,
          inv_gamma_shape = fit$shape,
          inv_gamma_scale = fit$scale
        )
      }
    }
  }
  if ("log_normal" %in% distributions) {
    fit <- fit_log_normal_variance(z)
    if (!is.null(fit)) {
      rows[[length(rows) + 1]] <- variance_distribution_row(
        distribution = "log-normal",
        n = length(x),
        logLik = fit$logLik,
        n_parameters = 2,
        log_location = fit$mu,
        log_scale = fit$sigma
      )
    }
  }
  if ("log_t" %in% distributions) {
    fit <- fit_log_t_variance(z)
    if (!is.null(fit)) {
      rows[[length(rows) + 1]] <- variance_distribution_row(
        distribution = "log-t",
        n = length(x),
        logLik = fit$logLik,
        n_parameters = 3,
        log_location = fit$mu,
        log_scale = fit$sigma,
        log_df = fit$df
      )
    }
  }
  if ("two_gaussian" %in% distributions) {
    fit <- fit_two_gaussian_variance(
      z,
      n_random = gmm_starts,
      min_weight = gmm_min_weight,
      sigma_floor = gmm_sigma_floor,
      seed = seed
    )
    if (!is.null(fit)) {
      rows[[length(rows) + 1]] <- variance_distribution_row(
        distribution = "two-Gaussian mixture",
        n = length(x),
        logLik = fit$logLik,
        n_parameters = 5,
        weight_low = fit$weight[1],
        weight_high = fit$weight[2],
        mu_low_log10 = fit$mu[1],
        mu_high_log10 = fit$mu[2],
        sigma_low_log10 = fit$sigma[1],
        sigma_high_log10 = fit$sigma[2],
        min_weight = fit$min_weight,
        sigma_floor = fit$sigma_floor,
        n_starts = fit$n_starts,
        em_iterations = fit$iterations
      )
    }
  }
  if (length(rows) == 0) {
    stop("No variance distribution could be fitted.", call. = FALSE)
  }

  out <- do.call(rbind, rows)
  out <- out[order(out$AIC), , drop = FALSE]
  rownames(out) <- NULL
  class(out) <- c("destress_variance_distribution_fit", class(out))
  attr(out, "variance") <- variance
  out
}

#' Plot diagnostic distributions for perturbation-level variances
#'
#' Shows the empirical density of log10 variances together with fitted
#' diagnostic distributions from [fit_variance_distribution()].
#'
#' @param diagnostics A data frame, usually from [perturbation_diagnostics()].
#' @param fit Optional output from [fit_variance_distribution()]. If omitted,
#'   distributions are fitted before plotting.
#' @param variance Numeric column containing positive variance estimates.
#' @param distributions Character vector passed to [fit_variance_distribution()]
#'   when `fit` is omitted.
#' @param adjust Bandwidth adjustment for the empirical density.
#' @param grid_n Number of grid points for fitted curves.
#' @param title,subtitle,xlab,ylab Plot labels.
#' @return A `ggplot` object. The fit table is available as `attr(plot, "fit")`.
#' @export
plot_variance_distribution <- function(diagnostics,
                                       fit = NULL,
                                       variance = "effect_variance",
                                       distributions = c("beta_prime", "inverse_gamma", "log_normal", "log_t"),
                                       adjust = 1.15,
                                       grid_n = 512,
                                       title = NULL,
                                       subtitle = NULL,
                                       xlab = NULL,
                                       ylab = "Density") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package `ggplot2` is required for plot_variance_distribution().", call. = FALSE)
  }
  stopifnot(is.data.frame(diagnostics))
  if (!variance %in% names(diagnostics)) {
    stop("`diagnostics` is missing variance column `", variance, "`.", call. = FALSE)
  }
  x <- as.numeric(diagnostics[[variance]])
  x <- x[is.finite(x) & x > 0]
  if (length(x) < 2) {
    stop("At least two finite positive variances are required for plotting.", call. = FALSE)
  }
  z <- log10(x)

  if (is.null(fit)) {
    fit <- fit_variance_distribution(
      diagnostics,
      variance = variance,
      distributions = distributions
    )
  }
  if (!is.data.frame(fit)) {
    stop("`fit` must be a data frame returned by fit_variance_distribution().", call. = FALSE)
  }
  if (is.null(xlab)) {
    xlab <- destress_effect_axis_label(variance, type = "log10_variance")
  }

  density_data <- data.frame(
    .log10_variance = z,
    stringsAsFactors = FALSE
  )
  z_grid <- seq(min(z, na.rm = TRUE), max(z, na.rm = TRUE), length.out = grid_n)
  fit_curve <- variance_distribution_curves(fit, z_grid)
  fit_curve <- fit_curve[is.finite(fit_curve$.density), , drop = FALSE]
  empirical_density <- stats::density(z, adjust = adjust)
  histogram_density <- graphics::hist(z, breaks = 24, plot = FALSE)$density
  y_max <- max(c(empirical_density$y, histogram_density), na.rm = TRUE)
  x_range <- range(z_grid, na.rm = TRUE)
  curve_order <- c(
    "Empirical KDE",
    "Beta-prime fit",
    "Inverse-gamma fit",
    "Log-normal fit",
    "Log-t fit",
    "Two-Gaussian mixture"
  )
  curve_labels <- c("Empirical KDE", intersect(curve_order[-1], unique(fit_curve$curve)))
  legend_data <- data.frame(
    curve = curve_labels,
    x = x_range[1] + 0.70 * diff(x_range),
    xend = x_range[1] + 0.78 * diff(x_range),
    y = NA_real_,
    yend = NA_real_,
    stringsAsFactors = FALSE
  )
  legend_data$y <- y_max * seq(1.01, by = -0.07, length.out = nrow(legend_data))
  legend_data$yend <- legend_data$y
  legend_text <- data.frame(
    curve = legend_data$curve,
    x = x_range[1] + 0.80 * diff(x_range),
    y = legend_data$y,
    stringsAsFactors = FALSE
  )

  p <- ggplot2::ggplot(density_data, ggplot2::aes(.log10_variance)) +
    ggplot2::geom_histogram(
      ggplot2::aes(y = ggplot2::after_stat(density)),
      bins = 24,
      fill = "#E5E7EB",
      color = "white",
      linewidth = 0.25,
      alpha = 0.85
    ) +
    ggplot2::geom_density(
      ggplot2::aes(color = "Empirical KDE", linetype = "Empirical KDE"),
      linewidth = 0.8,
      adjust = adjust,
      show.legend = FALSE
    ) +
    destress_diagnostic_theme(base_size = 10, legend_position = c(0.985, 0.985),
                              legend_justification = c(1, 1)) +
    ggplot2::scale_color_manual(
      values = c(
        "Empirical KDE" = "#111827",
        "Beta-prime fit" = "#0072B2",
        "Inverse-gamma fit" = "#009E73",
        "Log-normal fit" = "#CC79A7",
        "Log-t fit" = "#D55E00",
        "Two-Gaussian mixture" = "#7C3AED"
      ),
      name = NULL
    ) +
    ggplot2::scale_linetype_manual(
      values = c(
        "Empirical KDE" = "solid",
        "Beta-prime fit" = "longdash",
        "Inverse-gamma fit" = "dotdash",
        "Log-normal fit" = "dashed",
        "Log-t fit" = "twodash",
        "Two-Gaussian mixture" = "22"
      ),
      name = NULL
    ) +
    ggplot2::coord_cartesian(ylim = c(0, y_max * 1.08), clip = "off") +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = xlab,
      y = ylab
    )

  if (nrow(fit_curve) > 0) {
    p <- p +
      ggplot2::geom_line(
        data = fit_curve,
        ggplot2::aes(.log10_variance, .density, color = curve, linetype = curve),
        inherit.aes = FALSE,
        linewidth = 0.85,
        show.legend = FALSE
      )
  }

  p <- p +
    ggplot2::geom_segment(
      data = legend_data,
      ggplot2::aes(x = x, xend = xend, y = y, yend = yend, color = curve, linetype = curve),
      inherit.aes = FALSE,
      linewidth = 0.8,
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      data = legend_text,
      ggplot2::aes(x = x, y = y, label = curve),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 0.5,
      size = 3.0,
      color = "#111827"
    )

  attr(p, "fit") <- fit
  p
}

destress_diagnostic_theme <- function(base_size = 10,
                                      legend_position = c(0.985, 0.985),
                                      legend_justification = c(1, 1)) {
  ggplot2::theme_light(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(linewidth = 0.25, color = "#E5E7EB"),
      axis.title = ggplot2::element_text(color = "#111827"),
      axis.text = ggplot2::element_text(color = "#4B5563"),
      legend.position = legend_position,
      legend.justification = legend_justification,
      legend.background = ggplot2::element_rect(
        fill = grDevices::adjustcolor("white", alpha.f = 0.88),
        color = "#D1D5DB",
        linewidth = 0.25
      ),
      legend.key = ggplot2::element_blank(),
      legend.margin = ggplot2::margin(3, 5, 3, 5),
      legend.key.width = grid::unit(1.35, "lines"),
      legend.key.height = grid::unit(0.38, "lines"),
      plot.title = ggplot2::element_text(face = "bold", color = "#111827"),
      plot.subtitle = ggplot2::element_text(color = "#4B5563")
    )
}

destress_effect_axis_label <- function(effect, type = c("mean", "rank_abs_mean", "variance", "log10_variance")) {
  type <- match.arg(type)
  effect <- as.character(effect)
  if (type == "mean") {
    return(switch(
      effect,
      total_effect = ,
      mean_response = expression(hat(Delta * bar(y))[j]^"tot"),
      rank_adjusted_total_effect = expression(hat(Delta * bar(y))[j]^"tot,k"),
      specific_effect = ,
      effect = expression(hat(Delta * bar(y))[j]^"spec"),
      global_effect = ,
      rank_adjusted_global_effect = expression(hat(Delta * bar(y))[j]^"tot"),
      "Mean effect"
    ))
  }
  if (type == "rank_abs_mean") {
    return(switch(
      effect,
      total_effect = ,
      mean_response = expression("Rank of " * abs(hat(Delta * bar(y))[j]^"tot")),
      rank_adjusted_total_effect = expression("Rank of " * abs(hat(Delta * bar(y))[j]^"tot,k")),
      specific_effect = ,
      effect = expression("Rank of " * abs(hat(Delta * bar(y))[j]^"spec")),
      global_effect = ,
      rank_adjusted_global_effect = expression("Rank of " * abs(hat(Delta * bar(y))[j]^"tot")),
      "Rank of absolute mean effect"
    ))
  }
  if (type == "variance") {
    return(switch(
      effect,
      total_effect = ,
      mean_response = ,
      total_effect_variance = expression("Variance across promoters of " * hat(Delta * y)["\u00b7" * j]^"tot"),
      rank_adjusted_total_effect = ,
      rank_adjusted_total_effect_variance = expression("Variance across promoters of " * hat(Delta * y)["\u00b7" * j]^"tot,k"),
      specific_effect = ,
      effect = ,
      specific_effect_variance = expression("Variance across promoters of " * hat(Delta * y)["\u00b7" * j]^"spec"),
      global_effect = ,
      rank_adjusted_global_effect = expression("Variance across promoters of " * bar(Delta * y)[j]^"tot"),
      expression("Variance across reporters")
    ))
  }
  switch(
    effect,
    total_effect = ,
    mean_response = ,
    total_effect_variance = expression(log[10] * " variance across promoters of " * hat(Delta * y)["\u00b7" * j]^"tot"),
    rank_adjusted_total_effect = ,
    rank_adjusted_total_effect_variance = expression(log[10] * " variance across promoters of " * hat(Delta * y)["\u00b7" * j]^"tot,k"),
    specific_effect = ,
    effect = ,
    specific_effect_variance = expression(log[10] * " variance across promoters of " * hat(Delta * y)["\u00b7" * j]^"spec"),
    expression(log[10] * " variance across reporters")
  )
}

variance_distribution_row <- function(distribution,
                                      n,
                                      logLik,
                                      n_parameters,
                                      shape1 = NA_real_,
                                      shape2 = NA_real_,
                                      scale = NA_real_,
                                      raw_logLik = NA_real_,
                                      inv_gamma_shape = NA_real_,
                                      inv_gamma_scale = NA_real_,
                                      log_location = NA_real_,
                                      log_scale = NA_real_,
                                      log_df = NA_real_,
                                      weight_low = NA_real_,
                                      weight_high = NA_real_,
                                      mu_low_log10 = NA_real_,
                                      mu_high_log10 = NA_real_,
                                      sigma_low_log10 = NA_real_,
                                      sigma_high_log10 = NA_real_,
                                      min_weight = NA_real_,
                                      sigma_floor = NA_real_,
                                      n_starts = NA_integer_,
                                      em_iterations = NA_integer_) {
  data.frame(
    distribution = distribution,
    n = n,
    logLik = logLik,
    AIC = 2 * n_parameters - 2 * logLik,
    n_parameters = n_parameters,
    beta_prime_shape1 = shape1,
    beta_prime_shape2 = shape2,
    beta_prime_scale = scale,
    beta_prime_raw_logLik = raw_logLik,
    inverse_gamma_shape = inv_gamma_shape,
    inverse_gamma_scale = inv_gamma_scale,
    log_location_log10 = log_location,
    log_scale_log10 = log_scale,
    log_t_df = log_df,
    mixture_weight_low = weight_low,
    mixture_weight_high = weight_high,
    mixture_mu_low_log10 = mu_low_log10,
    mixture_mu_high_log10 = mu_high_log10,
    mixture_sigma_low_log10 = sigma_low_log10,
    mixture_sigma_high_log10 = sigma_high_log10,
    mixture_min_weight = min_weight,
    mixture_sigma_floor = sigma_floor,
    mixture_n_starts = n_starts,
    mixture_em_iterations = em_iterations,
    stringsAsFactors = FALSE
  )
}

fit_beta_prime_variance <- function(x) {
  x <- x[is.finite(x) & x > 0]
  if (length(x) < 10) {
    return(NULL)
  }
  lower <- log(c(1e-3, 1e-3, min(x) / 100))
  upper <- log(c(1e3, 1e3, max(x) * 100))
  neg_loglik <- function(theta) {
    shape1 <- exp(theta[1])
    shape2 <- exp(theta[2])
    scale <- exp(theta[3])
    y <- x / scale
    -sum(
      (shape1 - 1) * log(y) -
        (shape1 + shape2) * log1p(y) -
        lbeta(shape1, shape2) -
        log(scale)
    )
  }
  starts <- rbind(
    log(c(2, 2, stats::median(x))),
    log(c(1, 3, stats::median(x))),
    log(c(3, 3, mean(x))),
    log(c(1.5, 1.5, stats::quantile(x, 0.75, names = FALSE)))
  )
  fits <- lapply(seq_len(nrow(starts)), function(i) {
    try(
      stats::optim(
        par = starts[i, ],
        fn = neg_loglik,
        method = "L-BFGS-B",
        lower = lower,
        upper = upper,
        control = list(maxit = 1000)
      ),
      silent = TRUE
    )
  })
  ok <- vapply(fits, function(fit) {
    !inherits(fit, "try-error") && isTRUE(fit$convergence == 0) && is.finite(fit$value)
  }, logical(1))
  if (!any(ok)) {
    return(NULL)
  }
  ok_idx <- which(ok)
  fit <- fits[[ok_idx[which.min(vapply(fits[ok], `[[`, numeric(1), "value"))]]]
  shape1 <- exp(fit$par[1])
  shape2 <- exp(fit$par[2])
  scale <- exp(fit$par[3])
  if (any(!is.finite(c(shape1, shape2, scale)))) {
    return(NULL)
  }
  list(shape1 = shape1, shape2 = shape2, scale = scale, logLik = -fit$value)
}

beta_prime_log10_density <- function(z, shape1, shape2, scale, log = FALSE) {
  x <- 10^z
  y <- x / scale
  out <- (shape1 - 1) * log(y) -
    (shape1 + shape2) * log1p(y) -
    lbeta(shape1, shape2) -
    log(scale) +
    log(log(10)) +
    z * log(10)
  if (log) out else exp(out)
}

fit_inverse_gamma_variance <- function(x) {
  x <- x[is.finite(x) & x > 0]
  if (length(x) < 10) {
    return(NULL)
  }
  lower <- log(c(1e-3, min(x) / 100))
  upper <- log(c(1e3, max(x) * 100))
  neg_loglik <- function(theta) {
    shape <- exp(theta[1])
    scale <- exp(theta[2])
    -sum(
      shape * log(scale) -
        lgamma(shape) -
        (shape + 1) * log(x) -
        scale / x
    )
  }
  m <- mean(x)
  v <- stats::var(x)
  start_shape <- if (is.finite(v) && v > 0) {
    max(2.1, min(100, 2 + m^2 / v))
  } else {
    4
  }
  start_scale <- max(.Machine$double.eps, m * (start_shape - 1))
  starts <- rbind(
    log(c(start_shape, start_scale)),
    log(c(2.5, stats::median(x) * 1.5)),
    log(c(5, mean(x) * 4))
  )
  fits <- lapply(seq_len(nrow(starts)), function(i) {
    try(
      stats::optim(
        par = starts[i, ],
        fn = neg_loglik,
        method = "L-BFGS-B",
        lower = lower,
        upper = upper,
        control = list(maxit = 1000)
      ),
      silent = TRUE
    )
  })
  ok <- vapply(fits, function(fit) {
    !inherits(fit, "try-error") && isTRUE(fit$convergence == 0) && is.finite(fit$value)
  }, logical(1))
  if (!any(ok)) {
    return(NULL)
  }
  ok_idx <- which(ok)
  fit <- fits[[ok_idx[which.min(vapply(fits[ok], `[[`, numeric(1), "value"))]]]
  shape <- exp(fit$par[1])
  scale <- exp(fit$par[2])
  if (any(!is.finite(c(shape, scale)))) {
    return(NULL)
  }
  list(shape = shape, scale = scale, logLik = -fit$value)
}

inverse_gamma_log10_density <- function(z, shape, scale, log = FALSE) {
  x <- 10^z
  out <- shape * log(scale) -
    lgamma(shape) -
    (shape + 1) * log(x) -
    scale / x +
    log(log(10)) +
    z * log(10)
  if (log) out else exp(out)
}

fit_log_normal_variance <- function(z) {
  z <- z[is.finite(z)]
  if (length(z) < 2 || stats::sd(z) == 0) {
    return(NULL)
  }
  mu <- mean(z)
  sigma <- sqrt(mean((z - mu)^2))
  if (!is.finite(sigma) || sigma <= 0) {
    return(NULL)
  }
  list(
    mu = mu,
    sigma = sigma,
    logLik = sum(stats::dnorm(z, mean = mu, sd = sigma, log = TRUE))
  )
}

log_normal_log10_density <- function(z, mu, sigma, log = FALSE) {
  stats::dnorm(z, mean = mu, sd = sigma, log = log)
}

fit_log_t_variance <- function(z) {
  z <- z[is.finite(z)]
  if (length(z) < 10 || stats::sd(z) == 0) {
    return(NULL)
  }
  robust_sd <- stats::IQR(z) / 1.349
  if (!is.finite(robust_sd) || robust_sd <= 0) {
    robust_sd <- stats::sd(z)
  }
  neg_loglik <- function(theta) {
    mu <- theta[1]
    sigma <- exp(theta[2])
    df <- 2 + exp(theta[3])
    -sum(stats::dt((z - mu) / sigma, df = df, log = TRUE) - log(sigma))
  }
  starts <- rbind(
    c(stats::median(z), log(stats::sd(z)), log(6)),
    c(mean(z), log(stats::sd(z)), log(20)),
    c(stats::median(z), log(robust_sd), log(3))
  )
  fits <- lapply(seq_len(nrow(starts)), function(i) {
    try(
      stats::optim(
        par = starts[i, ],
        fn = neg_loglik,
        method = "BFGS",
        control = list(maxit = 1000)
      ),
      silent = TRUE
    )
  })
  ok <- vapply(fits, function(fit) {
    !inherits(fit, "try-error") && isTRUE(fit$convergence == 0) && is.finite(fit$value)
  }, logical(1))
  if (!any(ok)) {
    return(NULL)
  }
  ok_idx <- which(ok)
  fit <- fits[[ok_idx[which.min(vapply(fits[ok], `[[`, numeric(1), "value"))]]]
  list(
    mu = fit$par[1],
    sigma = exp(fit$par[2]),
    df = 2 + exp(fit$par[3]),
    logLik = -fit$value
  )
}

log_t_log10_density <- function(z, mu, sigma, df, log = FALSE) {
  out <- stats::dt((z - mu) / sigma, df = df, log = TRUE) - log(sigma)
  if (log) out else exp(out)
}

fit_two_gaussian_variance <- function(z,
                                      max_iter = 2000,
                                      tol = 1e-9,
                                      n_random = 50,
                                      min_weight = 0.08,
                                      sigma_floor = NULL,
                                      seed = 1) {
  z <- z[is.finite(z)]
  if (length(z) < 10 || stats::sd(z) == 0) {
    return(NULL)
  }
  k <- 2
  n <- length(z)
  sd_z <- stats::sd(z)
  if (is.null(sigma_floor)) {
    sigma_floor <- max(sd_z * 0.25, 0.08)
  }
  log_sum_exp <- function(mat) {
    row_max <- do.call(pmax, as.data.frame(mat))
    row_max + log(rowSums(exp(mat - row_max)))
  }
  run_em <- function(mu, sigma, weight) {
    weight <- pmax(weight, 1e-6)
    weight <- weight / sum(weight)
    sigma <- pmax(sigma, sigma_floor)
    loglik_old <- -Inf
    loglik <- -Inf
    for (iter in seq_len(max_iter)) {
      log_dens <- vapply(
        seq_len(k),
        function(component) {
          log(weight[component]) +
            stats::dnorm(z, mu[component], sigma[component], log = TRUE)
        },
        numeric(n)
      )
      log_denom <- log_sum_exp(log_dens)
      resp <- exp(log_dens - log_denom)
      nk <- pmax(colSums(resp), 1e-8)
      weight <- pmax(nk / n, min_weight)
      weight <- weight / sum(weight)
      mu <- colSums(resp * z) / nk
      centered <- sweep(matrix(z, nrow = n, ncol = k), 2, mu, "-")
      sigma <- sqrt(colSums(resp * centered^2) / nk)
      sigma <- pmax(sigma, sigma_floor)
      loglik <- sum(log_denom)
      if (is.finite(loglik_old) && abs(loglik - loglik_old) < tol) {
        break
      }
      loglik_old <- loglik
    }
    ord <- order(mu)
    list(
      weight = weight[ord],
      mu = mu[ord],
      sigma = sigma[ord],
      logLik = loglik,
      iterations = iter
    )
  }

  starts <- list(
    list(
      mu = as.numeric(stats::quantile(z, c(0.35, 0.85), names = FALSE)),
      sigma = rep(max(sd_z / 2, sigma_floor), k),
      weight = c(0.75, 0.25)
    ),
    list(
      mu = as.numeric(stats::quantile(z, c(0.25, 0.75), names = FALSE)),
      sigma = rep(max(sd_z / 2.5, sigma_floor), k),
      weight = c(0.5, 0.5)
    )
  )
  kmeans_fit <- try(stats::kmeans(z, centers = k, nstart = 50), silent = TRUE)
  if (!inherits(kmeans_fit, "try-error")) {
    clusters <- split(z, kmeans_fit$cluster)
    starts[[length(starts) + 1]] <- list(
      mu = vapply(clusters, mean, numeric(1)),
      sigma = pmax(vapply(clusters, stats::sd, numeric(1)), sigma_floor),
      weight = lengths(clusters) / n
    )
  }
  set.seed(seed)
  for (start_id in seq_len(n_random)) {
    starts[[length(starts) + 1]] <- list(
      mu = sort(sample(z, k)),
      sigma = rep(stats::runif(1, sd_z / 6, sd_z), k),
      weight = as.numeric(stats::rgamma(k, shape = 1))
    )
  }
  fits <- lapply(starts, function(start) run_em(start$mu, start$sigma, start$weight))
  logliks <- vapply(fits, `[[`, numeric(1), "logLik")
  best <- fits[[which.max(logliks)]]
  best$n_starts <- length(starts)
  best$min_weight <- min_weight
  best$sigma_floor <- sigma_floor
  best
}

two_gaussian_log10_density <- function(z, weight_low, weight_high,
                                       mu_low, mu_high, sigma_low, sigma_high) {
  weight_low * stats::dnorm(z, mu_low, sigma_low) +
    weight_high * stats::dnorm(z, mu_high, sigma_high)
}

variance_distribution_curves <- function(fit, z_grid) {
  rows <- list()
  for (i in seq_len(nrow(fit))) {
    distribution <- fit$distribution[i]
    if (identical(distribution, "beta-prime")) {
      rows[[length(rows) + 1]] <- data.frame(
        .log10_variance = z_grid,
        .density = beta_prime_log10_density(
          z_grid,
          shape1 = fit$beta_prime_shape1[i],
          shape2 = fit$beta_prime_shape2[i],
          scale = fit$beta_prime_scale[i]
        ),
        curve = "Beta-prime fit",
        stringsAsFactors = FALSE
      )
    }
    if (identical(distribution, "inverse-gamma")) {
      rows[[length(rows) + 1]] <- data.frame(
        .log10_variance = z_grid,
        .density = inverse_gamma_log10_density(
          z_grid,
          shape = fit$inverse_gamma_shape[i],
          scale = fit$inverse_gamma_scale[i]
        ),
        curve = "Inverse-gamma fit",
        stringsAsFactors = FALSE
      )
    }
    if (identical(distribution, "log-normal")) {
      rows[[length(rows) + 1]] <- data.frame(
        .log10_variance = z_grid,
        .density = log_normal_log10_density(
          z_grid,
          mu = fit$log_location_log10[i],
          sigma = fit$log_scale_log10[i]
        ),
        curve = "Log-normal fit",
        stringsAsFactors = FALSE
      )
    }
    if (identical(distribution, "log-t")) {
      rows[[length(rows) + 1]] <- data.frame(
        .log10_variance = z_grid,
        .density = log_t_log10_density(
          z_grid,
          mu = fit$log_location_log10[i],
          sigma = fit$log_scale_log10[i],
          df = fit$log_t_df[i]
        ),
        curve = "Log-t fit",
        stringsAsFactors = FALSE
      )
    }
    if (identical(distribution, "two-Gaussian mixture")) {
      rows[[length(rows) + 1]] <- data.frame(
        .log10_variance = z_grid,
        .density = two_gaussian_log10_density(
          z_grid,
          weight_low = fit$mixture_weight_low[i],
          weight_high = fit$mixture_weight_high[i],
          mu_low = fit$mixture_mu_low_log10[i],
          mu_high = fit$mixture_mu_high_log10[i],
          sigma_low = fit$mixture_sigma_low_log10[i],
          sigma_high = fit$mixture_sigma_high_log10[i]
        ),
        curve = "Two-Gaussian mixture",
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(rows) == 0) {
    return(data.frame())
  }
  do.call(rbind, rows)
}
