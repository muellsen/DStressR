#' Volcano plot for DStressR reporter-perturbation hits
#'
#' Creates a standard volcano plot from a DStressR result table. The x-axis is a
#' reporter-perturbation effect size and the y-axis is the negative log10 adjusted
#' p-value. Significant hits are emphasized, top reporter groups can be colored,
#' and the most significant reporter-perturbation pairs are annotated.
#'
#' The defaults work with [results()] followed by [adjust_pvalues()] or
#' [call_hits()]. For workflow comparison tables, pass the corresponding column
#' names, for example `effect = "destress_eb_effect_centered"`,
#' `padj = "estimated_alpha_eb_padj_by_reporter"`, `perturbation = "srn_code"`,
#' and `perturbation_label = "ProductName"`.
#'
#' @param table A data frame with one row per reporter-perturbation pair.
#' @param effect Effect-size column to plot on the x-axis.
#' @param padj Adjusted p-value column to plot on the y-axis.
#' @param pvalue Optional raw p-value column used only if `padj = NULL`.
#' @param reporter,perturbation Columns identifying the reporter and perturbation.
#' @param perturbation_label Optional column with human-readable perturbation names used
#'   for annotations. Defaults to `perturbation`.
#' @param fdr FDR threshold for hit highlighting.
#' @param lfc Minimum absolute effect size for hit highlighting.
#' @param top_n Number of significant pairs to annotate.
#' @param top_reporters Number of reporter groups to color. Remaining reporters
#'   are shown in grey.
#' @param title,subtitle Plot title and subtitle.
#' @param xlab,ylab Axis labels. Defaults to readable labels based on the
#'   selected columns.
#' @param label_by Label style for annotated points. The default, `"pair"`,
#'   labels top hits as reporter-perturbation pairs.
#' @param max_label_chars Maximum characters per annotation label. Longer
#'   labels are truncated with `...`. Use `Inf` to keep full labels.
#' @param repel_labels If `TRUE` and the optional `ggrepel` package is
#'   installed, use repelled labels for the annotated top hits.
#' @param point_alpha Point transparency.
#' @return A `ggplot` object.
#' @export
plot_volcano <- function(table,
                         effect = "specific_effect",
                         padj = "specific_padj",
                         pvalue = NULL,
                         reporter = "reporter",
                         perturbation = "perturbation",
                         perturbation_label = perturbation,
                         fdr = 0.05,
                         lfc = 0,
                         top_n = 12,
                         top_reporters = 6,
                         title = "DStressR volcano plot",
                         subtitle = NULL,
                         xlab = NULL,
                         ylab = NULL,
                         label_by = c("pair", "reporter", "perturbation"),
                         max_label_chars = 46,
                         repel_labels = TRUE,
                         point_alpha = 0.65) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package `ggplot2` is required for plot_volcano().", call. = FALSE)
  }
  stopifnot(is.data.frame(table))
  label_by <- match.arg(label_by)
  y_col <- if (!is.null(padj)) padj else pvalue
  if (is.null(y_col)) {
    stop("Provide either `padj` or `pvalue`.", call. = FALSE)
  }
  required <- c(effect, y_col, reporter, perturbation, perturbation_label)
  missing_cols <- setdiff(required, names(table))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  d <- table
  d$.effect <- as.numeric(d[[effect]])
  d$.p_for_plot <- as.numeric(d[[y_col]])
  d$.reporter <- as.character(d[[reporter]])
  d$.perturbation <- as.character(d[[perturbation]])
  d$.perturbation_label <- as.character(d[[perturbation_label]])
  missing_label <- is.na(d$.perturbation_label) | !nzchar(d$.perturbation_label)
  d$.perturbation_label[missing_label] <- d$.perturbation[missing_label]
  d <- d[is.finite(d$.effect) & is.finite(d$.p_for_plot) & d$.p_for_plot > 0, , drop = FALSE]
  if (nrow(d) == 0) {
    stop("No finite effect and p-value rows available for plotting.", call. = FALSE)
  }

  d$.neg_log10_p <- -log10(pmax(d$.p_for_plot, .Machine$double.xmin))
  d$.hit <- d$.p_for_plot < fdr & abs(d$.effect) >= lfc
  d$.direction <- ifelse(
    d$.hit & d$.effect > 0,
    "Up",
    ifelse(d$.hit & d$.effect < 0, "Down", "Not significant")
  )

  hit_counts <- stats::aggregate(
    d$.hit,
    by = list(reporter = d$.reporter),
    FUN = sum
  )
  names(hit_counts)[2] <- "hit_n"
  total_counts <- as.data.frame(table(d$.reporter), stringsAsFactors = FALSE)
  names(total_counts) <- c("reporter", "total_n")
  reporter_counts <- merge(hit_counts, total_counts, by = "reporter", all = TRUE)
  reporter_counts$hit_n[is.na(reporter_counts$hit_n)] <- 0
  reporter_counts <- reporter_counts[order(-reporter_counts$hit_n, -reporter_counts$total_n, reporter_counts$reporter), ]
  colored_reporters <- utils::head(reporter_counts$reporter, top_reporters)
  d$.reporter_group <- ifelse(d$.reporter %in% colored_reporters, d$.reporter, "Other reporters")
  d$.reporter_group <- factor(d$.reporter_group, levels = c(colored_reporters, "Other reporters"))

  label_df <- d[d$.hit, , drop = FALSE]
  label_df <- label_df[order(label_df$.p_for_plot, -abs(label_df$.effect)), , drop = FALSE]
  label_df <- utils::head(label_df, top_n)
  if (label_by == "pair") {
    label_df$.label <- paste(label_df$.reporter, label_df$.perturbation_label, sep = " / ")
  } else if (label_by == "reporter") {
    label_df$.label <- label_df$.reporter
  } else {
    label_df$.label <- label_df$.perturbation_label
  }
  if (is.finite(max_label_chars)) {
    too_long <- nchar(label_df$.label) > max_label_chars
    label_df$.label[too_long] <- paste0(
      substr(label_df$.label[too_long], 1, max_label_chars - 3),
      "..."
    )
  }

  palette <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9")
  values <- stats::setNames(rep(palette, length.out = length(colored_reporters)), colored_reporters)
  values <- c(values, "Other reporters" = "#C5C5C5")
  if (is.null(subtitle)) {
    subtitle <- paste0(
      "Hits: adjusted p < ", fdr,
      if (lfc > 0) paste0(" and |effect| >= ", lfc) else ""
    )
  }
  if (is.null(xlab)) {
    xlab <- paste0("Effect size: ", effect)
  }
  if (is.null(ylab)) {
    ylab <- paste0("-log10 adjusted p-value: ", y_col)
  }

  p <- ggplot2::ggplot(d, ggplot2::aes(x = .effect, y = .neg_log10_p)) +
    ggplot2::geom_point(
      ggplot2::aes(color = .reporter_group, shape = .direction, alpha = .direction),
      size = 1.8,
      stroke = 0.25
    ) +
    ggplot2::geom_hline(
      yintercept = -log10(fdr),
      linetype = "longdash",
      color = "#505050",
      linewidth = 0.35
    ) +
    ggplot2::geom_vline(
      xintercept = c(-lfc, lfc),
      linetype = if (lfc > 0) "longdash" else "blank",
      color = "#505050",
      linewidth = 0.35
    ) +
    ggplot2::scale_color_manual(values = values, drop = FALSE) +
    ggplot2::scale_shape_manual(values = c("Down" = 25, "Not significant" = 16, "Up" = 24)) +
    ggplot2::scale_alpha_manual(values = c("Down" = 0.95, "Not significant" = point_alpha, "Up" = 0.95), guide = "none") +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme_light() +
    ggplot2::theme(
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank(),
      plot.title.position = "plot",
      plot.margin = ggplot2::margin(8, 34, 8, 8)
    ) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = xlab,
      y = ylab,
      color = "Reporter",
      shape = "Hit"
    )

  if (nrow(label_df) > 0) {
    if (isTRUE(repel_labels) && requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p +
        ggrepel::geom_text_repel(
          data = label_df,
          ggplot2::aes(label = .label),
          size = 2.7,
          min.segment.length = 0,
          box.padding = 0.3,
          point.padding = 0.15,
          max.overlaps = Inf,
          seed = 1,
          show.legend = FALSE
        )
    } else {
      p <- p +
        ggplot2::geom_text(
          data = label_df,
          ggplot2::aes(label = .label),
          size = 2.7,
          vjust = -0.65,
          check_overlap = TRUE,
          show.legend = FALSE
        )
    }
  }

  p
}

make_response_matrix <- function(table, value, reporter, perturbation_display) {
  d <- table[, c(reporter, perturbation_display, value), drop = FALSE]
  names(d) <- c(".reporter", ".perturbation_display", ".value")
  d$.value <- as.numeric(d$.value)
  d <- d[is.finite(d$.value), , drop = FALSE]
  if (nrow(d) == 0) {
    stop("No finite response values available for plotting.", call. = FALSE)
  }
  d <- stats::aggregate(
    .value ~ .reporter + .perturbation_display,
    d,
    mean,
    na.rm = TRUE
  )
  wide <- stats::reshape(
    d,
    idvar = ".reporter",
    timevar = ".perturbation_display",
    direction = "wide"
  )
  names(wide) <- sub("^\\.value\\.", "", names(wide))
  rownames(wide) <- wide$.reporter
  wide$.reporter <- NULL
  as.matrix(wide)
}

destress_known_reporter_order <- c(
  "acrABp", "marRABp", "micFp", "ompFp", "robp", "soxSp", "tolCp"
)

destress_reporter_order <- function(reporters, reporter_order = NULL) {
  reporters <- unique(as.character(reporters))
  reporters <- reporters[!is.na(reporters) & nzchar(reporters)]
  if (length(reporters) == 0) {
    return(character())
  }

  if (is.null(reporter_order)) {
    reporter_order <- getOption("DStressR.reporter_order", NULL)
  }
  if (is.null(reporter_order)) {
    reporter_order <- destress_known_reporter_order
  }
  reporter_order <- unique(as.character(reporter_order))
  reporter_order <- reporter_order[!is.na(reporter_order) & nzchar(reporter_order)]

  c(
    reporter_order[reporter_order %in% reporters],
    sort(setdiff(reporters, reporter_order))
  )
}

destress_perturbation_axis_labels <- function(perturbations,
                                          show_perturbation_labels,
                                          top_n = 40,
                                          score = NULL,
                                          min_gap = NULL) {
  perturbations <- as.character(perturbations)
  if (isTRUE(show_perturbation_labels)) {
    return(perturbations)
  }
  if (isFALSE(show_perturbation_labels)) {
    return(rep("", length(perturbations)))
  }
  if (length(perturbations) <= top_n) {
    return(perturbations)
  }
  if (is.null(score)) {
    score <- seq_along(perturbations)
  }
  score <- as.numeric(score)
  if (length(score) != length(perturbations)) {
    stop("`score` must have one value per perturbation label.", call. = FALSE)
  }
  score[!is.finite(score)] <- -Inf
  labels <- rep("", length(perturbations))
  if (is.null(min_gap)) {
    min_gap <- max(1, ceiling(length(perturbations) / (top_n * 1.5)))
  }
  min_gap <- as.integer(min_gap)
  if (length(min_gap) != 1 || is.na(min_gap) || min_gap < 1) {
    stop("`min_gap` must be a positive integer.", call. = FALSE)
  }
  candidates <- order(-score, perturbations)
  keep <- integer()
  for (idx in candidates) {
    if (length(keep) >= top_n) {
      break
    }
    if (length(keep) == 0 || min(abs(idx - keep)) >= min_gap) {
      keep <- c(keep, idx)
    }
  }
  labels[keep] <- perturbations[keep]
  labels
}

matrix_cluster_order <- function(mat, margin) {
  if (margin == 1) {
    keep <- rowSums(is.finite(mat)) > 1
    if (sum(keep) < 2) return(seq_len(nrow(mat)))
    ord <- seq_len(nrow(mat))
    ord[keep] <- which(keep)[stats::hclust(stats::dist(mat[keep, , drop = FALSE]))$order]
    ord
  } else {
    keep <- colSums(is.finite(mat)) > 1
    if (sum(keep) < 2) return(seq_len(ncol(mat)))
    ord <- seq_len(ncol(mat))
    ord[keep] <- which(keep)[stats::hclust(stats::dist(t(mat[, keep, drop = FALSE])))$order]
    ord
  }
}

#' Heatmap of a DStressR reporter-by-perturbation response matrix
#'
#' Creates a standard heatmap for normalized reporter-perturbation responses. The
#' default `value` is `specific_effect`, matching [results()], but workflow
#' tables can use columns such as `destress_eb_effect_centered`.
#'
#' @param table A data frame with one row per reporter-perturbation pair.
#' @param value Numeric response/effect column to show in the heatmap.
#' @param reporter,perturbation Columns identifying reporters and perturbations.
#' @param perturbation_label Optional human-readable perturbation-name column. Defaults
#'   to `perturbation`.
#' @param show_perturbation_ids If `TRUE`, append perturbation IDs in square brackets
#'   to perturbation labels.
#' @param top_n_perturbations If finite, show only the top perturbations by mean
#'   absolute response. Use `Inf` to show all perturbations.
#' @param reporter_order Optional global reporter order used when
#'   `cluster_rows = FALSE`. If omitted, the option `DStressR.reporter_order` is
#'   used when set; otherwise known DStressR paper reporters are shown in their
#'   manuscript order and remaining reporters are sorted alphabetically.
#' @param cluster_rows,cluster_cols If `TRUE`, hierarchically cluster reporters
#'   and/or perturbations.
#' @param clip_quantile Quantile of absolute response values used to clip the
#'   color scale. Set to `1` to use the observed maximum.
#' @param color_limit Optional positive color-scale limit. If supplied, values
#'   are clipped to `[-color_limit, color_limit]`; otherwise the limit is
#'   computed from `clip_quantile`.
#' @param show_perturbation_labels If `TRUE`, draw all x-axis perturbation labels. If
#'   `FALSE`, suppress x-axis perturbation labels. The default labels the
#'   `top_perturbation_labels` perturbations with largest absolute column sums, or all
#'   perturbations when fewer are plotted.
#' @param top_perturbation_labels Number of highest-signal perturbations to label when
#'   `show_perturbation_labels = NULL`.
#' @param perturbation_label_score Optional numeric score used to choose the
#'   top-labelled perturbations when `show_perturbation_labels = NULL`. If named, values
#'   are matched to perturbation labels; otherwise the order must match the displayed
#'   matrix columns.
#' @param perturbation_label_min_gap Minimum number of matrix columns between
#'   automatically selected labels. The default chooses a gap from the displayed
#'   matrix size and `top_perturbation_labels`.
#' @param perturbation_label_angle Angle used for visible perturbation labels.
#' @param title,subtitle,xlab,ylab Plot labels.
#' @param legend_title Colorbar title. Defaults to the selected `value` column.
#' @param low,mid,high Colors for negative, zero, and positive responses.
#' @return A `ggplot` object. The plotted matrix is available as
#'   `attr(plot, "response_matrix")`.
#' @export
plot_response_heatmap <- function(table,
                                  value = "specific_effect",
                                  reporter = "reporter",
                                  perturbation = "perturbation",
                                  perturbation_label = perturbation,
                                  show_perturbation_ids = TRUE,
                                  top_n_perturbations = 160,
                                  reporter_order = NULL,
                                  cluster_rows = FALSE,
                                  cluster_cols = TRUE,
                                  clip_quantile = 0.98,
                                  color_limit = NULL,
                                  show_perturbation_labels = NULL,
                                  top_perturbation_labels = 40,
                                  perturbation_label_score = NULL,
                                  perturbation_label_min_gap = NULL,
                                  perturbation_label_angle = 45,
                                  title = "DStressR reporter-by-perturbation matrix",
                                  subtitle = NULL,
                                  xlab = "Perturbations",
                                  ylab = "Reporters",
                                  legend_title = value,
                                  low = "#2166AC",
                                  mid = "white",
                                  high = "#B2182B") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package `ggplot2` is required for plot_response_heatmap().", call. = FALSE)
  }
  stopifnot(is.data.frame(table))
  required <- c(value, reporter, perturbation, perturbation_label)
  missing_cols <- setdiff(required, names(table))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  d <- table
  d$.reporter <- as.character(d[[reporter]])
  d$.perturbation <- as.character(d[[perturbation]])
  d$.perturbation_label <- as.character(d[[perturbation_label]])
  missing_label <- is.na(d$.perturbation_label) | !nzchar(d$.perturbation_label)
  d$.perturbation_label[missing_label] <- d$.perturbation[missing_label]
  d$.perturbation_display <- d$.perturbation_label
  if (isTRUE(show_perturbation_ids)) {
    d$.perturbation_display <- paste0(d$.perturbation_label, " [", d$.perturbation, "]")
  }
  d$.value <- as.numeric(d[[value]])
  d <- d[is.finite(d$.value), , drop = FALSE]
  if (nrow(d) == 0) {
    stop("No finite response values available for plotting.", call. = FALSE)
  }

  if (is.finite(top_n_perturbations)) {
    perturbation_summary <- stats::aggregate(
      abs(.value) ~ .perturbation_display,
      d,
      mean,
      na.rm = TRUE
    )
    names(perturbation_summary)[2] <- ".mean_abs_value"
    perturbation_summary <- perturbation_summary[order(-perturbation_summary$.mean_abs_value), , drop = FALSE]
    keep_perturbations <- utils::head(perturbation_summary$.perturbation_display, top_n_perturbations)
    d <- d[d$.perturbation_display %in% keep_perturbations, , drop = FALSE]
  }

  mat <- make_response_matrix(d, ".value", ".reporter", ".perturbation_display")
  if (isTRUE(cluster_rows)) {
    mat <- mat[matrix_cluster_order(mat, 1), , drop = FALSE]
  } else {
    row_order <- destress_reporter_order(rownames(mat), reporter_order = reporter_order)
    mat <- mat[row_order, , drop = FALSE]
  }
  if (isTRUE(cluster_cols)) {
    mat <- mat[, matrix_cluster_order(mat, 2), drop = FALSE]
  }

  plot_df <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
  names(plot_df) <- c(".reporter", ".perturbation_display", ".value")
  plot_df$.reporter <- factor(plot_df$.reporter, levels = rev(rownames(mat)))
  plot_df$.perturbation_display <- factor(plot_df$.perturbation_display, levels = colnames(mat))

  if (is.null(color_limit)) {
    limit <- stats::quantile(abs(plot_df$.value), clip_quantile, na.rm = TRUE)
    if (!is.finite(limit) || limit <= 0) {
      limit <- max(abs(plot_df$.value), na.rm = TRUE)
    }
  } else {
    limit <- as.numeric(color_limit)
    if (length(limit) != 1 || !is.finite(limit) || limit <= 0) {
      stop("`color_limit` must be a positive finite number.", call. = FALSE)
    }
  }
  plot_df$.plot_value <- pmax(pmin(plot_df$.value, limit), -limit)
  if (is.null(subtitle)) {
    subtitle <- if (is.finite(top_n_perturbations)) {
      paste0("Top ", top_n_perturbations, " perturbations by mean absolute ", value)
    } else {
      "All perturbations"
    }
  }
  label_score <- colSums(abs(mat), na.rm = TRUE)
  if (!is.null(perturbation_label_score)) {
    perturbation_label_score <- as.numeric(perturbation_label_score)
    if (!is.null(names(perturbation_label_score))) {
      label_score <- perturbation_label_score[colnames(mat)]
    } else {
      if (length(perturbation_label_score) != ncol(mat)) {
        stop("Unnamed `perturbation_label_score` must have one value per displayed perturbation.", call. = FALSE)
      }
      label_score <- perturbation_label_score
    }
  }
  perturbation_axis_labels <- destress_perturbation_axis_labels(
    colnames(mat),
    show_perturbation_labels,
    top_n = top_perturbation_labels,
    score = label_score,
    min_gap = perturbation_label_min_gap
  )

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(.perturbation_display, .reporter, fill = .plot_value)
  ) +
    ggplot2::geom_raster() +
    ggplot2::scale_fill_gradient2(
      low = low,
      mid = mid,
      high = high,
      midpoint = 0,
      limits = c(-limit, limit),
      breaks = c(-limit, 0, limit),
      labels = signif(c(-limit, 0, limit), 2),
      na.value = "#e5e7eb",
      name = legend_title
    ) +
    ggplot2::scale_x_discrete(
      labels = stats::setNames(perturbation_axis_labels, colnames(mat))
    ) +
    ggplot2::theme_light(base_size = 8) +
    ggplot2::theme(
      axis.ticks.x = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      legend.position = "bottom",
      plot.title.position = "plot"
    ) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = xlab,
      y = ylab
    )
  p <- p + ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = perturbation_label_angle,
      hjust = 1,
      vjust = 1
    )
  )

  attr(p, "response_matrix") <- mat
  attr(p, "color_limit") <- limit
  p
}

hit_table_summary <- function(d, group_col) {
  groups <- split(d, d[[group_col]])
  out <- lapply(names(groups), function(g) {
    x <- groups[[g]]
    hit_x <- x[x$.hit, , drop = FALSE]
    data.frame(
      group = g,
      n_pairs = nrow(x),
      n_hits = sum(x$.hit, na.rm = TRUE),
      n_positive_hits = sum(x$.hit & x$.effect > 0, na.rm = TRUE),
      n_negative_hits = sum(x$.hit & x$.effect < 0, na.rm = TRUE),
      mean_abs_effect = mean(abs(x$.effect), na.rm = TRUE),
      mean_abs_hit_effect = if (nrow(hit_x) > 0) mean(abs(hit_x$.effect), na.rm = TRUE) else NA_real_,
      max_abs_effect = max(abs(x$.effect), na.rm = TRUE),
      min_padj = min(x$.padj, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, out)
  out <- out[order(-out$n_hits, -out$max_abs_effect, out$group), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Summarize significant reporter-perturbation hits
#'
#' Builds compact pair-, reporter-, and perturbation-level summaries from a DStressR
#' result table. Hits are defined by an adjusted p-value threshold and, optionally,
#' a minimum absolute effect size.
#'
#' @param table A data frame with one row per reporter-perturbation pair.
#' @param effect Effect-size column used for hit direction and effect summaries.
#' @param padj Adjusted p-value column used for hit calls.
#' @param reporter,perturbation Columns identifying reporters and perturbations.
#' @param perturbation_label Optional human-readable perturbation-name column. Defaults
#'   to `perturbation`.
#' @param fdr FDR threshold for hit calls.
#' @param lfc Minimum absolute effect size for hit calls.
#' @return A list of class `destress_hit_summary` with pair-level hits and
#'   reporter- and perturbation-level summaries.
#' @export
summarize_hits <- function(table,
                           effect = "specific_effect",
                           padj = "specific_padj_by_reporter",
                           reporter = "reporter",
                           perturbation = "perturbation",
                           perturbation_label = perturbation,
                           fdr = 0.05,
                           lfc = 0) {
  stopifnot(is.data.frame(table))
  required <- c(effect, padj, reporter, perturbation, perturbation_label)
  missing_cols <- setdiff(required, names(table))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  d <- table
  d$.effect <- as.numeric(d[[effect]])
  d$.padj <- as.numeric(d[[padj]])
  d$.reporter <- as.character(d[[reporter]])
  d$.perturbation <- as.character(d[[perturbation]])
  d$.perturbation_label <- as.character(d[[perturbation_label]])
  missing_label <- is.na(d$.perturbation_label) | !nzchar(d$.perturbation_label)
  d$.perturbation_label[missing_label] <- d$.perturbation[missing_label]
  d <- d[is.finite(d$.effect) & is.finite(d$.padj), , drop = FALSE]
  if (nrow(d) == 0) {
    stop("No finite effect and adjusted p-value rows available.", call. = FALSE)
  }

  d$.hit <- d$.padj <= fdr & abs(d$.effect) >= lfc
  d$.direction <- ifelse(
    d$.hit & d$.effect > 0,
    "positive",
    ifelse(d$.hit & d$.effect < 0, "negative", "not_significant")
  )

  pairs <- data.frame(
    reporter = d$.reporter,
    perturbation = d$.perturbation,
    perturbation_label = d$.perturbation_label,
    effect = d$.effect,
    padj = d$.padj,
    hit = d$.hit,
    direction = d$.direction,
    stringsAsFactors = FALSE
  )
  pairs <- pairs[order(pairs$padj, -abs(pairs$effect), pairs$reporter, pairs$perturbation_label), , drop = FALSE]
  rownames(pairs) <- NULL

  reporter_summary <- hit_table_summary(d, ".reporter")
  names(reporter_summary)[names(reporter_summary) == "group"] <- "reporter"
  perturbation_summary <- hit_table_summary(d, ".perturbation_label")
  names(perturbation_summary)[names(perturbation_summary) == "group"] <- "perturbation_label"

  out <- list(
    pairs = pairs,
    reporters = reporter_summary,
    perturbations = perturbation_summary,
    thresholds = list(fdr = fdr, lfc = lfc, effect = effect, padj = padj)
  )
  class(out) <- "destress_hit_summary"
  out
}

#' Heatmap of significant DStressR hits
#'
#' Shows a reporter-by-perturbation effect matrix with significant pairs highlighted
#' by color and non-significant pairs shown as a light background. This plot is a
#' compact companion to [plot_response_heatmap()] for inspecting the discovered
#' hit structure.
#'
#' @param table A data frame with one row per reporter-perturbation pair.
#' @param effect Effect-size column shown by color for significant hits.
#' @param padj Adjusted p-value column used for hit calls.
#' @param reporter,perturbation Columns identifying reporters and perturbations.
#' @param perturbation_label Optional human-readable perturbation-name column. Defaults
#'   to `perturbation`.
#' @param show_perturbation_ids If `TRUE`, append perturbation IDs in square brackets
#'   to perturbation labels.
#' @param top_n_perturbations If finite, show only perturbations with the strongest hit
#'   evidence, ranked by hit count and effect size. Use `Inf` to show all
#'   perturbations.
#' @param fdr FDR threshold for hit highlighting.
#' @param lfc Minimum absolute effect size for hit highlighting.
#' @param drop_empty_perturbations If `TRUE`, remove perturbations with no significant
#'   hits from the displayed matrix.
#' @param reporter_order Optional global reporter order used when
#'   `order_rows = "global"`. If omitted, the option `DStressR.reporter_order`
#'   is used when set; otherwise known DStressR paper reporters are shown in
#'   their manuscript order and remaining reporters are sorted alphabetically.
#' @param order_rows,order_cols Ordering strategy for reporters and perturbations.
#'   Use `"global"` for the package-wide reporter order, `"input"` to preserve
#'   the input/factor order, `"frequency"` to order by number of hits, or
#'   `"cluster"` for hierarchical clustering of the hit matrix.
#' @param clip_quantile Quantile of absolute significant effects used to clip the
#'   color scale. Set to `1` to use the observed maximum.
#' @param color_limit Optional positive color-scale limit. If supplied, values
#'   are clipped to `[-color_limit, color_limit]`; otherwise the limit is
#'   computed from `clip_quantile`.
#' @param show_perturbation_labels If `TRUE`, draw all x-axis perturbation labels. If
#'   `FALSE`, suppress x-axis perturbation labels. The default labels the
#'   `top_perturbation_labels` perturbations with largest absolute column sums, or all
#'   perturbations when fewer are plotted.
#' @param top_perturbation_labels Number of highest-signal perturbations to label when
#'   `show_perturbation_labels = NULL`.
#' @param perturbation_label_score Optional numeric score used to choose the
#'   top-labelled perturbations when `show_perturbation_labels = NULL`. If named, values
#'   are matched to perturbation labels; otherwise the order must match the displayed
#'   matrix columns.
#' @param perturbation_label_min_gap Minimum number of matrix columns between
#'   automatically selected labels. The default chooses a gap from the displayed
#'   matrix size and `top_perturbation_labels`.
#' @param perturbation_label_angle Angle used for visible perturbation labels.
#' @param title,subtitle,xlab,ylab Plot labels.
#' @param legend_title Colorbar title. Defaults to the selected `effect` column.
#' @param low,mid,high Colors for negative, zero, and positive hit effects.
#' @return A `ggplot` object with attributes `hit_matrix`, `hit_summary`,
#'   `plotted_pairs`, and `color_limit`.
#' @export
plot_hit_heatmap <- function(table,
                             effect = "specific_effect",
                             padj = "specific_padj_by_reporter",
                             reporter = "reporter",
                             perturbation = "perturbation",
                             perturbation_label = perturbation,
                             show_perturbation_ids = TRUE,
                             top_n_perturbations = 160,
                             fdr = 0.05,
                             lfc = 0,
                             drop_empty_perturbations = TRUE,
                             reporter_order = NULL,
                             order_rows = c("global", "cluster", "input", "frequency"),
                             order_cols = c("frequency", "cluster", "input"),
                             clip_quantile = 0.98,
                             color_limit = NULL,
                             show_perturbation_labels = NULL,
                             top_perturbation_labels = 40,
                             perturbation_label_score = NULL,
                             perturbation_label_min_gap = NULL,
                             perturbation_label_angle = 45,
                             title = "DStressR significant-hit matrix",
                             subtitle = NULL,
                             xlab = "Perturbations",
                             ylab = "Reporters",
                             legend_title = effect,
                             low = "#2166AC",
                             mid = "white",
                             high = "#B2182B") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package `ggplot2` is required for plot_hit_heatmap().", call. = FALSE)
  }
  stopifnot(is.data.frame(table))
  order_rows <- match.arg(order_rows)
  order_cols <- match.arg(order_cols)
  required <- c(effect, padj, reporter, perturbation, perturbation_label)
  missing_cols <- setdiff(required, names(table))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  d <- table
  reporter_levels <- if (is.factor(d[[reporter]])) levels(d[[reporter]]) else unique(as.character(d[[reporter]]))
  perturbation_levels <- if (is.factor(d[[perturbation_label]])) levels(d[[perturbation_label]]) else unique(as.character(d[[perturbation_label]]))
  d$.reporter <- as.character(d[[reporter]])
  d$.perturbation <- as.character(d[[perturbation]])
  d$.perturbation_label <- as.character(d[[perturbation_label]])
  missing_label <- is.na(d$.perturbation_label) | !nzchar(d$.perturbation_label)
  d$.perturbation_label[missing_label] <- d$.perturbation[missing_label]
  d$.perturbation_display <- d$.perturbation_label
  if (isTRUE(show_perturbation_ids)) {
    d$.perturbation_display <- paste0(d$.perturbation_label, " [", d$.perturbation, "]")
    perturbation_levels <- unique(d$.perturbation_display)
  }
  d$.effect <- as.numeric(d[[effect]])
  d$.padj <- as.numeric(d[[padj]])
  d <- d[is.finite(d$.effect) & is.finite(d$.padj), , drop = FALSE]
  if (nrow(d) == 0) {
    stop("No finite effect and adjusted p-value rows available for plotting.", call. = FALSE)
  }
  d$.hit <- d$.padj <= fdr & abs(d$.effect) >= lfc
  d$.cluster_value <- ifelse(d$.hit, d$.effect, 0)
  if (isTRUE(drop_empty_perturbations)) {
    hit_perturbations <- unique(d$.perturbation_display[d$.hit])
    d <- d[d$.perturbation_display %in% hit_perturbations, , drop = FALSE]
    if (nrow(d) == 0) {
      stop("No perturbations have significant hits under the selected thresholds.", call. = FALSE)
    }
  }

  if (is.finite(top_n_perturbations)) {
    perturbation_summary <- hit_table_summary(d, ".perturbation_display")
    names(perturbation_summary)[names(perturbation_summary) == "group"] <- ".perturbation_display"
    perturbation_summary <- perturbation_summary[
      order(
        -perturbation_summary$n_hits,
        -perturbation_summary$max_abs_effect,
        perturbation_summary$.perturbation_display
      ),
      ,
      drop = FALSE
    ]
    keep_perturbations <- utils::head(perturbation_summary$.perturbation_display, top_n_perturbations)
    d <- d[d$.perturbation_display %in% keep_perturbations, , drop = FALSE]
  }

  mat <- make_response_matrix(d, ".cluster_value", ".reporter", ".perturbation_display")
  input_reporters <- reporter_levels[reporter_levels %in% rownames(mat)]
  if (length(input_reporters) == 0) {
    input_reporters <- rownames(mat)
  }
  input_perturbations <- perturbation_levels[perturbation_levels %in% colnames(mat)]
  if (length(input_perturbations) == 0) {
    input_perturbations <- colnames(mat)
  }

  if (order_rows == "cluster") {
    row_order <- rownames(mat)[matrix_cluster_order(mat, 1)]
  } else if (order_rows == "frequency") {
    row_hits <- rowSums(mat != 0, na.rm = TRUE)
    row_order <- rownames(mat)[order(-row_hits, rownames(mat))]
  } else if (order_rows == "global") {
    row_order <- destress_reporter_order(rownames(mat), reporter_order = reporter_order)
  } else {
    row_order <- input_reporters
  }
  if (order_cols == "cluster") {
    col_order <- colnames(mat)[matrix_cluster_order(mat, 2)]
  } else if (order_cols == "frequency") {
    col_hits <- colSums(mat != 0, na.rm = TRUE)
    col_strength <- colMeans(abs(mat), na.rm = TRUE)
    col_order <- colnames(mat)[order(-col_hits, -col_strength, colnames(mat))]
  } else {
    col_order <- input_perturbations
  }
  mat <- mat[row_order, col_order, drop = FALSE]

  plot_df <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
  names(plot_df) <- c(".reporter", ".perturbation_display", ".plot_value")
  plot_df$.hit <- plot_df$.plot_value != 0
  plot_df$.reporter <- factor(plot_df$.reporter, levels = rev(rownames(mat)))
  plot_df$.perturbation_display <- factor(plot_df$.perturbation_display, levels = colnames(mat))
  plot_df$.fill_value <- ifelse(plot_df$.hit, plot_df$.plot_value, NA_real_)

  value_for_limit <- abs(plot_df$.fill_value[is.finite(plot_df$.fill_value)])
  if (length(value_for_limit) == 0) {
    value_for_limit <- abs(d$.effect)
  }
  if (is.null(color_limit)) {
    limit <- stats::quantile(value_for_limit, clip_quantile, na.rm = TRUE)
    if (!is.finite(limit) || limit <= 0) {
      limit <- max(value_for_limit, na.rm = TRUE)
    }
  } else {
    limit <- as.numeric(color_limit)
    if (length(limit) != 1 || !is.finite(limit) || limit <= 0) {
      stop("`color_limit` must be a positive finite number.", call. = FALSE)
    }
  }
  plot_df$.fill_value <- pmax(pmin(plot_df$.fill_value, limit), -limit)

  if (is.null(subtitle)) {
    subtitle <- paste0(
      "Hits: adjusted p <= ", fdr,
      if (lfc > 0) paste0(" and |effect| >= ", lfc) else ""
    )
  }
  label_score <- colSums(abs(mat), na.rm = TRUE)
  if (!is.null(perturbation_label_score)) {
    perturbation_label_score <- as.numeric(perturbation_label_score)
    if (!is.null(names(perturbation_label_score))) {
      label_score <- perturbation_label_score[colnames(mat)]
    } else {
      if (length(perturbation_label_score) != ncol(mat)) {
        stop("Unnamed `perturbation_label_score` must have one value per displayed perturbation.", call. = FALSE)
      }
      label_score <- perturbation_label_score
    }
  }
  perturbation_axis_labels <- destress_perturbation_axis_labels(
    colnames(mat),
    show_perturbation_labels,
    top_n = top_perturbation_labels,
    score = label_score,
    min_gap = perturbation_label_min_gap
  )

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(.perturbation_display, .reporter)) +
    ggplot2::geom_tile(fill = "#F1F3F5", color = "white", linewidth = 0.18) +
    ggplot2::geom_tile(ggplot2::aes(fill = .fill_value), color = "white", linewidth = 0.18) +
    ggplot2::scale_fill_gradient2(
      low = low,
      mid = mid,
      high = high,
      midpoint = 0,
      limits = c(-limit, limit),
      breaks = c(-limit, 0, limit),
      labels = signif(c(-limit, 0, limit), 2),
      na.value = "#F1F3F5",
      name = legend_title
    ) +
    ggplot2::scale_x_discrete(
      labels = stats::setNames(perturbation_axis_labels, colnames(mat))
    ) +
    ggplot2::theme_light(base_size = 8) +
    ggplot2::theme(
      axis.ticks.x = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      legend.position = "bottom",
      plot.title.position = "plot"
    ) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = xlab,
      y = ylab
    )
  p <- p + ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = perturbation_label_angle,
      hjust = 1,
      vjust = 1
    )
  )

  attr(p, "hit_matrix") <- mat
  attr(p, "hit_summary") <- summarize_hits(
    table,
    effect = effect,
    padj = padj,
    reporter = reporter,
    perturbation = perturbation,
    perturbation_label = perturbation_label,
    fdr = fdr,
    lfc = lfc
  )
  attr(p, "plotted_pairs") <- d
  attr(p, "color_limit") <- limit
  p
}

#' Histogram of DStressR reporter-perturbation effects
#'
#' Shows the empirical distribution of normalized reporter-perturbation effects,
#' either over all matrix entries or faceted by reporter.
#'
#' @param table A data frame with one row per reporter-perturbation pair.
#' @param value Numeric effect column to plot.
#' @param reporter Column identifying reporters, used when `by = "reporter"`.
#' @param by Plot one pooled histogram (`"all"`) or reporter-faceted
#'   histograms (`"reporter"`).
#' @param bins Number of histogram bins.
#' @param xlim Optional two-element x-axis limit.
#' @param scales Facet scale behavior for `by = "reporter"`.
#' @param title,subtitle,xlab,ylab Plot labels.
#' @param fill,border Histogram fill and border colors.
#' @return A `ggplot` object.
#' @export
plot_effect_histogram <- function(table,
                                  value = "specific_effect",
                                  reporter = "reporter",
                                  by = c("all", "reporter"),
                                  bins = 80,
                                  xlim = NULL,
                                  scales = "fixed",
                                  title = NULL,
                                  subtitle = NULL,
                                  xlab = NULL,
                                  ylab = "Count",
                                  fill = "#4E79A7",
                                  border = "white") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package `ggplot2` is required for plot_effect_histogram().", call. = FALSE)
  }
  stopifnot(is.data.frame(table))
  by <- match.arg(by)
  required <- value
  if (by == "reporter") {
    required <- c(required, reporter)
  }
  missing_cols <- setdiff(required, names(table))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  d <- table
  d$.effect_hist <- as.numeric(d[[value]])
  d <- d[is.finite(d$.effect_hist), , drop = FALSE]
  if (nrow(d) == 0) {
    stop("No finite effect values available for plotting.", call. = FALSE)
  }
  if (by == "reporter") {
    d$.reporter <- as.character(d[[reporter]])
  }
  if (is.null(xlab)) {
    xlab <- value
  }
  if (is.null(title)) {
    title <- if (by == "reporter") {
      "Effect distributions by reporter"
    } else {
      "Effect distribution over all reporter-perturbation entries"
    }
  }
  if (is.null(subtitle)) {
    subtitle <- paste0(
      "n = ", nrow(d),
      "; median = ", signif(stats::median(d$.effect_hist), 3),
      "; MAD = ", signif(stats::mad(d$.effect_hist), 3)
    )
  }

  p <- ggplot2::ggplot(d, ggplot2::aes(.effect_hist)) +
    ggplot2::geom_histogram(bins = bins, fill = fill, color = border, linewidth = 0.15) +
    ggplot2::geom_vline(xintercept = 0, color = "#303030", linewidth = 0.35) +
    ggplot2::theme_light(base_size = 10) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title.position = "plot"
    ) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = xlab,
      y = ylab
    )
  if (!is.null(xlim)) {
    p <- p + ggplot2::coord_cartesian(xlim = xlim)
  }
  if (by == "reporter") {
    p <- p +
      ggplot2::facet_wrap(ggplot2::vars(.reporter), scales = scales, ncol = 5) +
      ggplot2::theme(strip.text = ggplot2::element_text(size = 8))
  }
  p
}

#' Clustered block map of a DStressR reporter-by-perturbation response matrix
#'
#' Hierarchically clusters reporters and perturbations, cuts the dendrograms into
#' interpretable groups, and plots the mean response for each reporter-cluster by
#' perturbation-cluster block. This is useful as a compact overview when the full
#' perturbation library is too large for individual perturbation labels.
#'
#' @param table A data frame with one row per reporter-perturbation pair.
#' @param value Numeric response/effect column to summarize.
#' @param reporter,perturbation Columns identifying reporters and perturbations.
#' @param perturbation_label Optional human-readable perturbation-name column. Defaults
#'   to `perturbation`.
#' @param show_perturbation_ids If `TRUE`, append perturbation IDs in square brackets
#'   to perturbation labels before clustering.
#' @param n_reporter_clusters,n_perturbation_clusters Number of dendrogram clusters
#'   to use for reporters and perturbations.
#' @param missing_value Value used only for clustering missing matrix entries.
#'   Block summaries are still computed from observed finite values.
#' @param clip_quantile Quantile of absolute block means used to clip the color
#'   scale. Set to `1` to use the observed maximum.
#' @param show_counts If `TRUE`, annotate each tile with the number of perturbations
#'   in that perturbation cluster.
#' @param title,subtitle,xlab,ylab Plot labels.
#' @param low,mid,high Colors for negative, zero, and positive responses.
#' @return A `ggplot` object with attributes `response_matrix`,
#'   `reporter_clusters`, `perturbation_clusters`, `block_summary`, `row_hclust`,
#'   and `col_hclust`.
#' @export
plot_response_cluster_blocks <- function(table,
                                         value = "specific_effect",
                                         reporter = "reporter",
                                         perturbation = "perturbation",
                                         perturbation_label = perturbation,
                                         show_perturbation_ids = TRUE,
                                         n_reporter_clusters = 6,
                                         n_perturbation_clusters = 14,
                                         missing_value = 0,
                                         clip_quantile = 0.98,
                                         show_counts = TRUE,
                                         title = "DStressR clustered response map",
                                         subtitle = NULL,
                                         xlab = "Perturbation clusters",
                                         ylab = "Reporter clusters",
                                         low = "#2166AC",
                                         mid = "white",
                                         high = "#B2182B") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package `ggplot2` is required for plot_response_cluster_blocks().", call. = FALSE)
  }
  stopifnot(is.data.frame(table))
  required <- c(value, reporter, perturbation, perturbation_label)
  missing_cols <- setdiff(required, names(table))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  d <- table
  d$.reporter <- as.character(d[[reporter]])
  d$.perturbation <- as.character(d[[perturbation]])
  d$.perturbation_label <- as.character(d[[perturbation_label]])
  missing_label <- is.na(d$.perturbation_label) | !nzchar(d$.perturbation_label)
  d$.perturbation_label[missing_label] <- d$.perturbation[missing_label]
  d$.perturbation_display <- d$.perturbation_label
  if (isTRUE(show_perturbation_ids)) {
    d$.perturbation_display <- paste0(d$.perturbation_label, " [", d$.perturbation, "]")
  }
  d$.value <- as.numeric(d[[value]])
  d <- d[is.finite(d$.value), , drop = FALSE]
  if (nrow(d) == 0) {
    stop("No finite response values available for plotting.", call. = FALSE)
  }

  mat <- make_response_matrix(d, ".value", ".reporter", ".perturbation_display")
  if (nrow(mat) < 2 || ncol(mat) < 2) {
    stop("Clustered block plots require at least two reporters and two perturbations.", call. = FALSE)
  }
  if (n_reporter_clusters < 1 || n_reporter_clusters > nrow(mat)) {
    stop("`n_reporter_clusters` must be between 1 and the number of reporters.", call. = FALSE)
  }
  if (n_perturbation_clusters < 1 || n_perturbation_clusters > ncol(mat)) {
    stop("`n_perturbation_clusters` must be between 1 and the number of perturbations.", call. = FALSE)
  }

  cluster_mat <- mat
  cluster_mat[!is.finite(cluster_mat)] <- missing_value
  row_hc <- stats::hclust(stats::dist(cluster_mat))
  col_hc <- stats::hclust(stats::dist(t(cluster_mat)))
  row_cluster <- stats::cutree(row_hc, k = n_reporter_clusters)[rownames(mat)]
  col_cluster <- stats::cutree(col_hc, k = n_perturbation_clusters)[colnames(mat)]

  reporter_assignments <- data.frame(
    reporter = rownames(mat),
    reporter_cluster = paste0("P", row_cluster),
    dendrogram_order = match(seq_along(rownames(mat)), row_hc$order),
    stringsAsFactors = FALSE
  )
  reporter_assignments <- reporter_assignments[
    order(row_cluster, reporter_assignments$dendrogram_order),
    ,
    drop = FALSE
  ]

  perturbation_assignments <- data.frame(
    perturbation_display = colnames(mat),
    perturbation_cluster = paste0("C", col_cluster),
    dendrogram_order = match(seq_along(colnames(mat)), col_hc$order),
    mean_abs_effect = colMeans(abs(mat), na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  perturbation_assignments <- perturbation_assignments[
    order(col_cluster, perturbation_assignments$dendrogram_order),
    ,
    drop = FALSE
  ]

  block_rows <- list()
  row_levels <- paste0("P", sort(unique(row_cluster)))
  col_levels <- paste0("C", sort(unique(col_cluster)))
  for (pc in sort(unique(row_cluster))) {
    for (cc in sort(unique(col_cluster))) {
      sub <- mat[row_cluster == pc, col_cluster == cc, drop = FALSE]
      block_rows[[length(block_rows) + 1]] <- data.frame(
        .reporter_cluster = paste0("P", pc),
        .perturbation_cluster = paste0("C", cc),
        n_reporters = sum(row_cluster == pc),
        n_perturbations = sum(col_cluster == cc),
        mean_effect = mean(sub, na.rm = TRUE),
        median_effect = stats::median(sub, na.rm = TRUE),
        mean_abs_effect = mean(abs(sub), na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
  }
  block_df <- do.call(rbind, block_rows)
  block_df$.reporter_cluster <- factor(block_df$.reporter_cluster, levels = row_levels)
  block_df$.perturbation_cluster <- factor(block_df$.perturbation_cluster, levels = col_levels)
  block_df$.count_label <- paste0("n=", block_df$n_perturbations)

  limit <- stats::quantile(abs(block_df$mean_effect), clip_quantile, na.rm = TRUE)
  if (!is.finite(limit) || limit <= 0) {
    limit <- max(abs(block_df$mean_effect), na.rm = TRUE)
  }
  block_df$.plot_value <- pmax(pmin(block_df$mean_effect, limit), -limit)
  if (is.null(subtitle)) {
    subtitle <- paste0(
      n_reporter_clusters, " reporter clusters x ", n_perturbation_clusters,
      " perturbation clusters"
    )
  }

  p <- ggplot2::ggplot(
    block_df,
    ggplot2::aes(.perturbation_cluster, .reporter_cluster, fill = .plot_value)
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.4) +
    ggplot2::scale_fill_gradient2(
      low = low,
      mid = mid,
      high = high,
      midpoint = 0,
      limits = c(-limit, limit),
      name = paste0("Mean\n", value)
    ) +
    ggplot2::theme_light(base_size = 9) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      legend.position = "bottom",
      plot.title.position = "plot"
    ) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = xlab,
      y = ylab
    )

  if (isTRUE(show_counts)) {
    p <- p +
      ggplot2::geom_text(
        label = block_df$.count_label,
        size = 2.6,
        color = "#334155"
      )
  }

  attr(p, "response_matrix") <- mat
  attr(p, "reporter_clusters") <- reporter_assignments
  attr(p, "perturbation_clusters") <- perturbation_assignments
  attr(p, "block_summary") <- block_df
  attr(p, "row_hclust") <- row_hc
  attr(p, "col_hclust") <- col_hc
  attr(p, "color_limit") <- limit
  p
}

draw_clustered_heatmap_base <- function(mat,
                                        row_hc,
                                        col_hc,
                                        row_clusters,
                                        col_clusters,
                                        color_limit,
                                        low,
                                        mid,
                                        high,
                                        title,
                                        subtitle,
                                        legend_title,
                                        show_rownames,
                                        show_colnames) {
  ordered_mat <- mat[row_hc$order, col_hc$order, drop = FALSE]
  plot_mat <- ordered_mat
  plot_mat[!is.finite(plot_mat)] <- NA_real_
  plot_mat <- pmax(pmin(plot_mat, color_limit), -color_limit)

  heat_cols <- grDevices::colorRampPalette(c(low, mid, high))(101)
  cluster_cols <- grDevices::colorRampPalette(c(
    "#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00",
    "#56B4E9", "#999999", "#F0E442", "#332288", "#88CCEE",
    "#44AA99", "#117733", "#882255", "#AA4499", "#DDCC77"
  ))(max(length(unique(row_clusters)), length(unique(col_clusters)), 3))

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::layout(
    matrix(c(0, 1, 0, 2, 3, 4, 0, 5, 0), nrow = 3, byrow = TRUE),
    widths = c(1.5, 8.5, 0.85),
    heights = c(1.55, 6.4, 0.35)
  )

  heat_left_margin <- if (isTRUE(show_rownames)) 5.2 else 0.4
  heat_bottom_margin <- if (isTRUE(show_colnames)) 6.5 else 0.4
  heat_right_margin <- 0.2

  graphics::par(mar = c(0, heat_left_margin, 2.6, heat_right_margin), xaxs = "i", yaxs = "i")
  graphics::plot(stats::as.dendrogram(col_hc), axes = FALSE, leaflab = "none")
  graphics::title(main = title, sub = subtitle, cex.main = 1.1, cex.sub = 0.78, line = 1.2)

  graphics::par(mar = c(0.4, 0.2, 0, 0), xaxs = "i", yaxs = "i")
  graphics::plot(stats::as.dendrogram(row_hc), horiz = TRUE, axes = FALSE, leaflab = "none")

  graphics::par(
    mar = c(heat_bottom_margin, heat_left_margin, 0, heat_right_margin),
    xaxs = "i",
    yaxs = "i"
  )
  graphics::image(
    x = seq_len(ncol(plot_mat)),
    y = seq_len(nrow(plot_mat)),
    z = t(plot_mat[nrow(plot_mat):1, , drop = FALSE]),
    col = heat_cols,
    zlim = c(-color_limit, color_limit),
    axes = FALSE,
    xlab = "",
    ylab = ""
  )
  if (isTRUE(show_rownames)) {
    graphics::axis(
      2,
      at = seq_len(nrow(plot_mat)),
      labels = rev(rownames(plot_mat)),
      las = 2,
      cex.axis = 0.62,
      tick = FALSE,
      line = -0.25
    )
  }
  if (isTRUE(show_colnames)) {
    graphics::axis(
      1,
      at = seq_len(ncol(plot_mat)),
      labels = colnames(plot_mat),
      las = 2,
      cex.axis = 0.42,
      tick = FALSE,
      line = -0.25
    )
  }
  graphics::box(col = "#334155", lwd = 0.6)

  graphics::par(mar = c(1.2, 0.4, 1.2, 3.5), xaxs = "i", yaxs = "i")
  legend_values <- seq(-color_limit, color_limit, length.out = length(heat_cols))
  graphics::image(
    x = 1,
    y = legend_values,
    z = matrix(legend_values, nrow = 1),
    col = heat_cols,
    zlim = c(-color_limit, color_limit),
    axes = FALSE,
    xlab = "",
    ylab = ""
  )
  graphics::axis(4, las = 1, cex.axis = 0.65)
  graphics::mtext(legend_title, side = 4, line = 2.4, cex = 0.7)

  graphics::par(mar = c(0, heat_left_margin, 0.08, heat_right_margin), xaxs = "i", yaxs = "i")
  ordered_col_clusters <- as.integer(factor(col_clusters[colnames(ordered_mat)]))
  graphics::image(
    x = seq_len(ncol(ordered_mat)),
    y = 1,
    z = matrix(ordered_col_clusters, nrow = ncol(ordered_mat), ncol = 1),
    col = cluster_cols,
    axes = FALSE,
    xlab = "",
    ylab = ""
  )
}

#' Clustered heatmap with reporter and perturbation dendrograms
#'
#' Draws a clustered reporter-by-perturbation response heatmap with hierarchical
#' trees on both axes. Unlike [plot_response_cluster_blocks()], this keeps the
#' individual matrix cells visible and uses the dendrograms to reveal structure
#' without collapsing the data into coarse blocks.
#'
#' @param table A data frame with one row per reporter-perturbation pair.
#' @param value Numeric response/effect column to show in the heatmap.
#' @param reporter,perturbation Columns identifying reporters and perturbations.
#' @param perturbation_label Optional human-readable perturbation-name column. Defaults
#'   to `perturbation`.
#' @param show_perturbation_ids If `TRUE`, append perturbation IDs in square brackets
#'   to perturbation labels before clustering.
#' @param top_n_perturbations If finite, show only the top perturbations by mean
#'   absolute response. Use `Inf` to show all perturbations.
#' @param n_reporter_clusters,n_perturbation_clusters Number of dendrogram clusters
#'   returned in the cluster assignment tables.
#' @param missing_value Value used only for clustering missing matrix entries.
#'   Heatmap cells with missing values are left missing.
#' @param clip_quantile Quantile of absolute response values used to clip the
#'   color scale. Set to `1` to use the observed maximum.
#' @param file Optional output file. Supports `.png` and `.pdf`. If `NULL`, the
#'   plot is drawn on the active graphics device.
#' @param width,height Plot size in inches when `file` is supplied.
#' @param res PNG resolution in dots per inch.
#' @param title,subtitle Plot title and subtitle.
#' @param show_rownames,show_colnames Whether to draw row and column labels.
#' @param low,mid,high Colors for negative, zero, and positive responses.
#' @return Invisibly returns a list containing the response matrix, clustering
#'   objects, ordered matrix, cluster assignments, and color limit.
#' @export
plot_response_clustered_heatmap <- function(table,
                                            value = "specific_effect",
                                            reporter = "reporter",
                                            perturbation = "perturbation",
                                            perturbation_label = perturbation,
                                            show_perturbation_ids = TRUE,
                                            top_n_perturbations = 400,
                                            n_reporter_clusters = 6,
                                            n_perturbation_clusters = 14,
                                            missing_value = 0,
                                            clip_quantile = 0.98,
                                            file = NULL,
                                            width = 14,
                                            height = 8,
                                            res = 300,
                                            title = "DStressR clustered response heatmap",
                                            subtitle = NULL,
                                            show_rownames = TRUE,
                                            show_colnames = FALSE,
                                            low = "#2166AC",
                                            mid = "white",
                                            high = "#B2182B") {
  stopifnot(is.data.frame(table))
  required <- c(value, reporter, perturbation, perturbation_label)
  missing_cols <- setdiff(required, names(table))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  d <- table
  d$.reporter <- as.character(d[[reporter]])
  d$.perturbation <- as.character(d[[perturbation]])
  d$.perturbation_label <- as.character(d[[perturbation_label]])
  missing_label <- is.na(d$.perturbation_label) | !nzchar(d$.perturbation_label)
  d$.perturbation_label[missing_label] <- d$.perturbation[missing_label]
  d$.perturbation_display <- d$.perturbation_label
  if (isTRUE(show_perturbation_ids)) {
    d$.perturbation_display <- paste0(d$.perturbation_label, " [", d$.perturbation, "]")
  }
  d$.value <- as.numeric(d[[value]])
  d <- d[is.finite(d$.value), , drop = FALSE]
  if (nrow(d) == 0) {
    stop("No finite response values available for plotting.", call. = FALSE)
  }

  if (is.finite(top_n_perturbations)) {
    perturbation_summary <- stats::aggregate(
      abs(.value) ~ .perturbation_display,
      d,
      mean,
      na.rm = TRUE
    )
    names(perturbation_summary)[2] <- ".mean_abs_value"
    perturbation_summary <- perturbation_summary[order(-perturbation_summary$.mean_abs_value), , drop = FALSE]
    keep_perturbations <- utils::head(perturbation_summary$.perturbation_display, top_n_perturbations)
    d <- d[d$.perturbation_display %in% keep_perturbations, , drop = FALSE]
  }

  mat <- make_response_matrix(d, ".value", ".reporter", ".perturbation_display")
  if (nrow(mat) < 2 || ncol(mat) < 2) {
    stop("Clustered heatmaps require at least two reporters and two perturbations.", call. = FALSE)
  }

  n_reporter_clusters <- min(max(1, n_reporter_clusters), nrow(mat))
  n_perturbation_clusters <- min(max(1, n_perturbation_clusters), ncol(mat))
  cluster_mat <- mat
  cluster_mat[!is.finite(cluster_mat)] <- missing_value
  row_hc <- stats::hclust(stats::dist(cluster_mat))
  col_hc <- stats::hclust(stats::dist(t(cluster_mat)))
  row_clusters <- stats::cutree(row_hc, k = n_reporter_clusters)[rownames(mat)]
  col_clusters <- stats::cutree(col_hc, k = n_perturbation_clusters)[colnames(mat)]
  ordered_mat <- mat[row_hc$order, col_hc$order, drop = FALSE]

  color_limit <- stats::quantile(abs(mat), clip_quantile, na.rm = TRUE)
  if (!is.finite(color_limit) || color_limit <= 0) {
    color_limit <- max(abs(mat), na.rm = TRUE)
  }
  if (is.null(subtitle)) {
    subtitle <- if (is.finite(top_n_perturbations)) {
      paste0("Top ", top_n_perturbations, " perturbations by mean absolute ", value)
    } else {
      "All perturbations"
    }
  }

  reporter_assignments <- data.frame(
    reporter = rownames(mat),
    reporter_cluster = paste0("P", row_clusters),
    dendrogram_order = match(seq_along(rownames(mat)), row_hc$order),
    stringsAsFactors = FALSE
  )
  reporter_assignments <- reporter_assignments[
    order(row_clusters, reporter_assignments$dendrogram_order),
    ,
    drop = FALSE
  ]
  perturbation_assignments <- data.frame(
    perturbation_display = colnames(mat),
    perturbation_cluster = paste0("C", col_clusters),
    dendrogram_order = match(seq_along(colnames(mat)), col_hc$order),
    mean_abs_effect = colMeans(abs(mat), na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  perturbation_assignments <- perturbation_assignments[
    order(col_clusters, perturbation_assignments$dendrogram_order),
    ,
    drop = FALSE
  ]

  if (!is.null(file)) {
    ext <- tolower(tools::file_ext(file))
    if (ext == "pdf") {
      grDevices::pdf(file, width = width, height = height, onefile = FALSE)
    } else if (ext == "png") {
      grDevices::png(file, width = width, height = height, units = "in", res = res, bg = "white")
    } else {
      stop("`file` must end in `.png` or `.pdf`.", call. = FALSE)
    }
    on.exit(grDevices::dev.off(), add = TRUE)
  }

  draw_clustered_heatmap_base(
    mat = mat,
    row_hc = row_hc,
    col_hc = col_hc,
    row_clusters = row_clusters,
    col_clusters = col_clusters,
    color_limit = color_limit,
    low = low,
    mid = mid,
    high = high,
    title = title,
    subtitle = subtitle,
    legend_title = value,
    show_rownames = show_rownames,
    show_colnames = show_colnames
  )

  invisible(list(
    response_matrix = mat,
    ordered_matrix = ordered_mat,
    row_hclust = row_hc,
    col_hclust = col_hc,
    reporter_clusters = reporter_assignments,
    perturbation_clusters = perturbation_assignments,
    color_limit = color_limit
  ))
}

utils::globalVariables(c(
  ".fill_value",
  ".perturbation_cluster",
  ".perturbation_display",
  ".direction",
  ".effect",
  ".effect_hist",
  ".label",
  ".neg_log10_p",
  ".plot_value",
  ".reporter",
  ".reporter_cluster",
  ".reporter_group"
))
