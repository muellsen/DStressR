# Avoid R CMD check notes for data-frame columns used inside ggplot2 aesthetics.
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    ".color_effect", ".component", ".density", ".effect", ".effect_variance",
    ".group_label", ".log10_variance", ".matrix", ".null_reference", ".rank", ".tail_region",
    ".signed_mean_effect", "curve", "density", "effect", "null_median",
    "null_q99", "observed", "perturbation_label", "reporter_order",
    "tail_probability", "variance_trend", "x", "xend", "y", "yend"
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
#' and perturbations, although DStressR result tables usually contain reporters
#' and perturbations. Large absolute mean effects with low variance indicate
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
                                     reporter = "reporter",
                                     perturbation = "perturbation",
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
                                          reporter = "reporter",
                                          perturbation = "perturbation",
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

split_table_by_group <- function(table, group = NULL) {
  if (is.null(group) || length(group) == 0) {
    out <- list(table)
    names(out) <- ""
    return(out)
  }
  missing_group <- setdiff(group, names(table))
  if (length(missing_group) > 0) {
    stop("Missing grouping columns: ", paste(missing_group, collapse = ", "), call. = FALSE)
  }
  key <- do.call(paste, c(table[group], sep = "\r"))
  split(table, key, drop = TRUE)
}

group_label_from_table <- function(table, group = NULL) {
  if (is.null(group) || length(group) == 0) {
    return("")
  }
  vals <- vapply(group, function(col) as.character(table[[col]][1]), character(1))
  paste(paste(group, vals, sep = " = "), collapse = ", ")
}

effect_matrix_from_table_mean <- function(table, effect, reporter, perturbation) {
  d <- table[, c(reporter, perturbation, effect), drop = FALSE]
  names(d) <- c(".reporter", ".perturbation", ".effect")
  d$.effect <- as.numeric(d$.effect)
  d <- d[is.finite(d$.effect), , drop = FALSE]
  if (nrow(d) == 0) {
    stop("No finite effects available for matrix decomposition.", call. = FALSE)
  }
  d <- stats::aggregate(.effect ~ .reporter + .perturbation, d, mean, na.rm = TRUE)
  reporters <- unique(as.character(d$.reporter))
  perturbations <- unique(as.character(d$.perturbation))
  mat <- matrix(
    NA_real_,
    nrow = length(reporters),
    ncol = length(perturbations),
    dimnames = list(reporters, perturbations)
  )
  idx <- cbind(match(as.character(d$.reporter), reporters), match(as.character(d$.perturbation), perturbations))
  mat[idx] <- d$.effect
  mat
}

#' Decompose a reporter-by-perturbation effect matrix
#'
#' Estimates a low-rank background component from an effect matrix and returns
#' the observed effect, low-rank component, and rank-adjusted residual in long
#' format. If grouping columns are supplied, the decomposition is performed
#' independently within each group.
#'
#' @param table A data frame with one row per reporter-perturbation pair.
#' @param effect Numeric effect column to decompose.
#' @param reporter,perturbation Columns identifying reporters and perturbations.
#' @param group Optional grouping columns. A separate low-rank decomposition is
#'   fitted within each group.
#' @param rank Non-negative rank of the background component.
#' @param impute Method used to fill missing entries before singular-value
#'   decomposition. The default `"column_mean"` replaces missing entries by the
#'   observed mean of the corresponding perturbation column.
#' @param reporter_label,perturbation_label Optional display-label columns.
#' @return A data frame containing the original effect, the low-rank effect, the
#'   rank-adjusted effect, and rank-1 reporter/perturbation scores when
#'   `rank >= 1`.
#' @export
low_rank_effect_decomposition <- function(table,
                                          effect = "total_effect",
                                          reporter = "reporter",
                                          perturbation = "perturbation",
                                          group = NULL,
                                          rank = 1,
                                          impute = c("column_mean", "global_mean", "zero"),
                                          reporter_label = reporter,
                                          perturbation_label = perturbation) {
  stopifnot(is.data.frame(table))
  rank <- validate_background_rank(rank)
  impute <- validate_background_impute(impute)
  required <- unique(c(effect, reporter, perturbation, group, reporter_label, perturbation_label))
  missing_cols <- setdiff(required, names(table))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  groups <- split_table_by_group(table, group = group)
  rows <- lapply(groups, function(dg) {
    mat <- effect_matrix_from_table_mean(dg, effect, reporter, perturbation)
    observed <- is.finite(mat)
    x <- impute_effect_matrix(mat, method = impute)
    rank_use <- min(rank, nrow(x), ncol(x))

    low_rank <- matrix(0, nrow = nrow(x), ncol = ncol(x), dimnames = dimnames(x))
    reporter_score_rank1 <- rep(NA_real_, nrow(x))
    perturbation_score_rank1 <- rep(NA_real_, ncol(x))
    names(reporter_score_rank1) <- rownames(x)
    names(perturbation_score_rank1) <- colnames(x)
    if (rank_use > 0 && any(abs(x[observed]) > sqrt(.Machine$double.eps))) {
      sv <- svd(x, nu = rank_use, nv = rank_use)
      keep <- seq_len(rank_use)
      low_rank <- sv$u[, keep, drop = FALSE] %*%
        (diag(sv$d[keep], nrow = rank_use, ncol = rank_use) %*% t(sv$v[, keep, drop = FALSE]))
      dimnames(low_rank) <- dimnames(x)
      reporter_score_rank1 <- sv$u[, 1] * sv$d[1]
      perturbation_score_rank1 <- sv$v[, 1] * sv$d[1]
      names(reporter_score_rank1) <- rownames(x)
      names(perturbation_score_rank1) <- colnames(x)
    }
    low_rank[!observed] <- NA_real_

    out <- expand.grid(
      reporter = rownames(mat),
      perturbation = colnames(mat),
      stringsAsFactors = FALSE
    )
    idx <- cbind(match(out$reporter, rownames(mat)), match(out$perturbation, colnames(mat)))
    out$effect <- mat[idx]
    out$low_rank_effect <- low_rank[idx]
    out$rank_adjusted_effect <- out$effect - out$low_rank_effect
    out$reporter_score_rank1 <- reporter_score_rank1[out$reporter]
    out$perturbation_score_rank1 <- perturbation_score_rank1[out$perturbation]
    out$rank <- rank
    out$group_label <- group_label_from_table(dg, group)
    if (!is.null(group) && length(group) > 0) {
      for (col in rev(group)) {
        out[[col]] <- dg[[col]][1]
        out <- out[, c(col, setdiff(names(out), col)), drop = FALSE]
      }
    }

    reporter_labels <- unique(dg[, c(reporter, reporter_label), drop = FALSE])
    names(reporter_labels) <- c("reporter", "reporter_label")
    perturbation_labels <- unique(dg[, c(perturbation, perturbation_label), drop = FALSE])
    names(perturbation_labels) <- c("perturbation", "perturbation_label")
    out <- merge(out, reporter_labels, by = "reporter", all.x = TRUE, sort = FALSE)
    out <- merge(out, perturbation_labels, by = "perturbation", all.x = TRUE, sort = FALSE)
    out
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Plot low-rank background diagnostics
#'
#' Draws observed singular values together with permutation-null summaries from
#' [background_rank_diagnostics()]. This is a diagnostic for choosing whether a
#' low-rank background component is visible in the effect matrix.
#'
#' @param table Optional effect table. Ignored when `diagnostics` is supplied.
#' @param diagnostics Optional data frame returned by
#'   [background_rank_diagnostics()].
#' @inheritParams background_rank_diagnostics
#' @param group Optional grouping columns. A separate diagnostic is computed and
#'   faceted for each group.
#' @param threshold Permutation reference quantile.
#' @param impute Method used to fill missing entries before singular-value
#'   decomposition.
#' @param title,subtitle,xlab,ylab Plot labels.
#' @return A `ggplot` object. The diagnostic table is available as
#'   `attr(plot, "diagnostics")`.
#' @export
plot_background_rank_diagnostics <- function(table = NULL,
                                             diagnostics = NULL,
                                             effect = "total_effect",
                                             reporter = "reporter",
                                             perturbation = "perturbation",
                                             group = NULL,
                                             rank_max = 10,
                                             permutations = 100,
                                             threshold = 0.99,
                                             impute = c("column_mean", "global_mean", "zero"),
                                             seed = NULL,
                                             title = NULL,
                                             subtitle = NULL,
                                             xlab = "Component",
                                             ylab = "Singular value") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package `ggplot2` is required for plot_background_rank_diagnostics().", call. = FALSE)
  }
  impute <- validate_background_impute(impute)
  if (is.null(diagnostics)) {
    if (is.null(table)) {
      stop("Provide either `table` or `diagnostics`.", call. = FALSE)
    }
    split_groups <- split_table_by_group(table, group = group)
    diagnostics <- do.call(rbind, lapply(seq_along(split_groups), function(i) {
      dg <- split_groups[[i]]
      seed_i <- if (is.null(seed)) NULL else seed + i - 1L
      out <- background_rank_diagnostics(
        dg,
        effect = effect,
        reporter = reporter,
        perturbation = perturbation,
        rank_max = rank_max,
        permutations = permutations,
        threshold = threshold,
        impute = impute,
        seed = seed_i
      )
      out$group_label <- group_label_from_table(dg, group)
      out
    }))
  }
  required <- c("component", "observed", "null_median", "null_q99")
  missing_cols <- setdiff(required, names(diagnostics))
  if (length(missing_cols) > 0) {
    stop("`diagnostics` is missing columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  d <- diagnostics
  if (!"group_label" %in% names(d)) {
    d$group_label <- ""
  }
  d$.null_reference <- if ("null_threshold" %in% names(d)) d$null_threshold else d$null_q99
  d$.group_label <- d$group_label
  d$.component <- as.numeric(d$component)

  p <- ggplot2::ggplot(d, ggplot2::aes(.component, observed)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = null_median, ymax = .null_reference),
      fill = "#BDBDBD",
      alpha = 0.35
    ) +
    ggplot2::geom_line(ggplot2::aes(y = null_median), color = "#737373", linetype = "dashed", linewidth = 0.35) +
    ggplot2::geom_line(color = "#1F78B4", linewidth = 0.45) +
    ggplot2::geom_point(color = "#1F78B4", size = 1.8) +
    ggplot2::scale_x_continuous(breaks = sort(unique(d$.component))) +
    destress_diagnostic_theme(base_size = 9, legend_position = "none") +
    ggplot2::labs(title = title, subtitle = subtitle, x = xlab, y = ylab)
  if (length(unique(d$.group_label)) > 1 || nzchar(d$.group_label[1])) {
    p <- p + ggplot2::facet_wrap(~.group_label)
  }
  attr(p, "diagnostics") <- diagnostics
  p
}

#' Plot a low-rank effect decomposition
#'
#' Shows the observed effect matrix, the estimated low-rank component, and the
#' rank-adjusted residual matrix. Rows are ordered by the first reporter score
#' within each group, making broad background axes visible even when reporter
#' labels are too dense to display.
#'
#' @param decomposition Output from [low_rank_effect_decomposition()].
#' @param matrices Which matrices to display.
#' @param group Optional grouping columns used for faceting.
#' @param show_reporter_labels If `TRUE`, draw reporter labels.
#' @param perturbation_order Optional perturbation order.
#' @param clip_quantile Quantile of absolute effects used for color clipping.
#' @param color_limit Optional positive color limit.
#' @param title,subtitle,xlab,ylab,legend_title Plot labels.
#' @return A `ggplot` object. The plotted data are available as
#'   `attr(plot, "plot_data")`.
#' @export
plot_low_rank_effect_heatmap <- function(decomposition,
                                         matrices = c("effect", "low_rank_effect", "rank_adjusted_effect"),
                                         group = NULL,
                                         show_reporter_labels = FALSE,
                                         perturbation_order = NULL,
                                         clip_quantile = 0.985,
                                         color_limit = NULL,
                                         title = NULL,
                                         subtitle = NULL,
                                         xlab = "Perturbation",
                                         ylab = "Reporters ordered by rank-1 score",
                                         legend_title = "Effect") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package `ggplot2` is required for plot_low_rank_effect_heatmap().", call. = FALSE)
  }
  stopifnot(is.data.frame(decomposition))
  matrices <- match.arg(matrices, c("effect", "low_rank_effect", "rank_adjusted_effect"), several.ok = TRUE)
  required <- c("reporter", "perturbation", "reporter_score_rank1", matrices)
  missing_cols <- setdiff(required, names(decomposition))
  if (length(missing_cols) > 0) {
    stop("`decomposition` is missing columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  d <- decomposition
  if (!"group_label" %in% names(d)) {
    d$group_label <- ""
  }
  if (!"reporter_label" %in% names(d)) {
    d$reporter_label <- d$reporter
  }
  if (!"perturbation_label" %in% names(d)) {
    d$perturbation_label <- d$perturbation
  }
  long <- do.call(rbind, lapply(matrices, function(mat_name) {
    data.frame(
      group_label = d$group_label,
      reporter = d$reporter,
      reporter_label = d$reporter_label,
      perturbation = d$perturbation,
      perturbation_label = d$perturbation_label,
      reporter_score_rank1 = d$reporter_score_rank1,
      matrix = mat_name,
      effect = as.numeric(d[[mat_name]]),
      stringsAsFactors = FALSE
    )
  }))
  matrix_labels <- c(
    effect = "Observed effect",
    low_rank_effect = "Low-rank component",
    rank_adjusted_effect = "Rank-adjusted residual"
  )
  long$.matrix <- factor(matrix_labels[long$matrix], levels = matrix_labels[matrices])
  if (is.null(perturbation_order)) {
    perturbation_order <- unique(long$perturbation_label)
  }
  long$perturbation_label <- factor(long$perturbation_label, levels = perturbation_order)

  order_key <- unique(long[, c("group_label", "reporter", "reporter_label", "reporter_score_rank1")])
  order_key <- order_key[order(order_key$group_label, order_key$reporter_score_rank1, order_key$reporter), , drop = FALSE]
  order_key$reporter_order <- stats::ave(seq_len(nrow(order_key)), order_key$group_label, FUN = seq_along)
  long$reporter_order <- order_key$reporter_order[match(paste(long$group_label, long$reporter), paste(order_key$group_label, order_key$reporter))]

  finite_effect <- abs(long$effect[is.finite(long$effect)])
  if (is.null(color_limit)) {
    limit <- stats::quantile(finite_effect, probs = clip_quantile, names = FALSE, na.rm = TRUE)
    if (!is.finite(limit) || limit <= 0) {
      limit <- max(finite_effect, na.rm = TRUE)
    }
  } else {
    limit <- as.numeric(color_limit)
  }
  if (!is.finite(limit) || limit <= 0) {
    limit <- 1
  }
  long$.effect <- pmax(pmin(long$effect, limit), -limit)

  p <- ggplot2::ggplot(long, ggplot2::aes(perturbation_label, reporter_order, fill = .effect)) +
    ggplot2::geom_raster() +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-limit, limit),
      breaks = c(-limit, 0, limit),
      labels = signif(c(-limit, 0, limit), 2),
      name = legend_title
    ) +
    ggplot2::facet_grid(group_label ~ .matrix, scales = "free_y", space = "free_y") +
    ggplot2::theme_light(base_size = 8) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
      axis.text.y = if (isTRUE(show_reporter_labels)) ggplot2::element_text(size = 5) else ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      legend.position = "bottom",
      plot.title.position = "plot",
      strip.background = ggplot2::element_rect(fill = "#F5F5F5", color = "#D0D0D0"),
      strip.text = ggplot2::element_text(face = "bold", color = "#333333")
    ) +
    ggplot2::labs(title = title, subtitle = subtitle, x = xlab, y = ylab)

  attr(p, "plot_data") <- long
  p
}

#' Robust diagnostic tail scores for effect distributions
#'
#' Computes robust Gaussian tail probabilities for effect values, optionally
#' within strata. These values are diagnostics for unreplicated or exploratory
#' effect matrices and should not be interpreted as formal p-values.
#'
#' @param table A data frame containing an effect column.
#' @param effect Numeric effect column.
#' @param group Optional grouping columns. Centers and scales are estimated
#'   separately within each group.
#' @param min_n Minimum number of finite effects required within a group.
#' @param alternative Tail alternative: `"two_sided"`, `"greater"`, or
#'   `"less"`.
#' @return The input table with diagnostic center, scale, z-score, tail
#'   probability, and negative log10 tail score columns appended.
#' @export
effect_tail_scores <- function(table,
                               effect = "effect",
                               group = NULL,
                               min_n = 5,
                               alternative = c("two_sided", "greater", "less")) {
  stopifnot(is.data.frame(table))
  alternative <- match.arg(alternative)
  required <- c(effect, group)
  missing_cols <- setdiff(required, names(table))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  d <- table
  d$diagnostic_center <- NA_real_
  d$diagnostic_scale <- NA_real_
  d$diagnostic_z <- NA_real_
  d$tail_probability <- NA_real_
  d$tail_score <- NA_real_
  d$.destress_row_id <- seq_len(nrow(d))
  group_list <- split_table_by_group(d, group = group)
  for (idx in group_list) {
    row_idx <- idx$.destress_row_id
    x <- as.numeric(idx[[effect]])
    finite <- is.finite(x)
    if (sum(finite) < min_n) {
      next
    }
    center <- stats::median(x[finite], na.rm = TRUE)
    scale <- stats::mad(x[finite], center = center, constant = 1.4826, na.rm = TRUE)
    if (!is.finite(scale) || scale <= sqrt(.Machine$double.eps)) {
      scale <- stats::sd(x[finite], na.rm = TRUE)
    }
    if (!is.finite(scale) || scale <= sqrt(.Machine$double.eps)) {
      next
    }
    z <- (x - center) / scale
    tail_probability <- switch(
      alternative,
      two_sided = 2 * stats::pnorm(abs(z), lower.tail = FALSE),
      greater = stats::pnorm(z, lower.tail = FALSE),
      less = stats::pnorm(z, lower.tail = TRUE)
    )
    d$diagnostic_center[row_idx] <- center
    d$diagnostic_scale[row_idx] <- scale
    d$diagnostic_z[row_idx] <- z
    d$tail_probability[row_idx] <- tail_probability
    d$tail_score[row_idx] <- -log10(pmax(tail_probability, .Machine$double.xmin))
  }
  d$.destress_row_id <- NULL
  d
}

#' Plot diagnostic effect-tail histograms
#'
#' Shows effect distributions and highlights observations in the diagnostic
#' tails. The plot is intended for exploratory effect matrices, especially when
#' formal replicate-based inference is unavailable.
#'
#' @param table A data frame, usually returned by [effect_tail_scores()].
#' @param effect Numeric effect column to plot.
#' @param tail_probability Diagnostic tail-probability column.
#' @param facet Optional columns used for faceting.
#' @param tail_threshold Tail-probability threshold used for highlighting.
#' @param bins Number of histogram bins.
#' @param title,subtitle,xlab,ylab Plot labels.
#' @return A `ggplot` object.
#' @export
plot_effect_tail_histogram <- function(table,
                                       effect = "effect",
                                       tail_probability = "tail_probability",
                                       facet = NULL,
                                       tail_threshold = 0.05,
                                       bins = 40,
                                       title = NULL,
                                       subtitle = NULL,
                                       xlab = NULL,
                                       ylab = "Effects") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package `ggplot2` is required for plot_effect_tail_histogram().", call. = FALSE)
  }
  stopifnot(is.data.frame(table))
  required <- c(effect, tail_probability, facet)
  missing_cols <- setdiff(required, names(table))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  d <- table
  d$.effect <- as.numeric(d[[effect]])
  d$.tail_probability <- as.numeric(d[[tail_probability]])
  d <- d[is.finite(d$.effect), , drop = FALSE]
  d$.tail_region <- ifelse(is.finite(d$.tail_probability) & d$.tail_probability <= tail_threshold, "Diagnostic tail", "Central")
  if (is.null(xlab)) {
    xlab <- effect
  }
  p <- ggplot2::ggplot(d, ggplot2::aes(.effect)) +
    ggplot2::geom_histogram(
      ggplot2::aes(fill = .tail_region),
      bins = bins,
      color = "white",
      linewidth = 0.15
    ) +
    ggplot2::scale_fill_manual(values = c("Central" = "#D1D5DB", "Diagnostic tail" = "#D55E00"), name = NULL) +
    destress_diagnostic_theme(base_size = 8, legend_position = c(0.985, 0.985), legend_justification = c(1, 1)) +
    ggplot2::labs(title = title, subtitle = subtitle, x = xlab, y = ylab)
  if (!is.null(facet) && length(facet) > 0) {
    p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", paste(facet, collapse = "+"))))
  }
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
      total_effect_variance = expression("Variance across reporters of " * hat(Delta * y)["\u00b7" * j]^"tot"),
      rank_adjusted_total_effect = ,
      rank_adjusted_total_effect_variance = expression("Variance across reporters of " * hat(Delta * y)["\u00b7" * j]^"tot,k"),
      specific_effect = ,
      effect = ,
      specific_effect_variance = expression("Variance across reporters of " * hat(Delta * y)["\u00b7" * j]^"spec"),
      global_effect = ,
      rank_adjusted_global_effect = expression("Variance across reporters of " * bar(Delta * y)[j]^"tot"),
      expression("Variance across reporters")
    ))
  }
  switch(
    effect,
    total_effect = ,
    mean_response = ,
    total_effect_variance = expression(log[10] * " variance across reporters of " * hat(Delta * y)["\u00b7" * j]^"tot"),
    rank_adjusted_total_effect = ,
    rank_adjusted_total_effect_variance = expression(log[10] * " variance across reporters of " * hat(Delta * y)["\u00b7" * j]^"tot,k"),
    specific_effect = ,
    effect = ,
    specific_effect_variance = expression(log[10] * " variance across reporters of " * hat(Delta * y)["\u00b7" * j]^"spec"),
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
