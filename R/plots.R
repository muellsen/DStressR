#' Volcano plot for DStressR promoter-compound hits
#'
#' Creates a standard volcano plot from a DStressR result table. The x-axis is a
#' promoter-compound effect size and the y-axis is the negative log10 adjusted
#' p-value. Significant hits are emphasized, top promoter groups can be colored,
#' and the most significant promoter-compound pairs are annotated.
#'
#' The defaults work with [results()] followed by [adjust_pvalues()] or
#' [call_hits()]. For workflow comparison tables, pass the corresponding column
#' names, for example `effect = "destress_eb_effect_centered"`,
#' `padj = "estimated_alpha_eb_padj_by_promoter"`, `compound = "srn_code"`,
#' and `compound_label = "ProductName"`.
#'
#' @param table A data frame with one row per promoter-compound pair.
#' @param effect Effect-size column to plot on the x-axis.
#' @param padj Adjusted p-value column to plot on the y-axis.
#' @param pvalue Optional raw p-value column used only if `padj = NULL`.
#' @param promoter,compound Columns identifying the promoter and compound.
#' @param compound_label Optional column with human-readable compound names used
#'   for annotations. Defaults to `compound`.
#' @param fdr FDR threshold for hit highlighting.
#' @param lfc Minimum absolute effect size for hit highlighting.
#' @param top_n Number of significant pairs to annotate.
#' @param top_promoters Number of promoter groups to color. Remaining promoters
#'   are shown in grey.
#' @param title,subtitle Plot title and subtitle.
#' @param xlab,ylab Axis labels. Defaults to readable labels based on the
#'   selected columns.
#' @param label_by Label style for annotated points. The default, `"pair"`,
#'   labels top hits as promoter-compound pairs.
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
                         promoter = "promoter",
                         compound = "compound",
                         compound_label = compound,
                         fdr = 0.05,
                         lfc = 0,
                         top_n = 12,
                         top_promoters = 6,
                         title = "DStressR volcano plot",
                         subtitle = NULL,
                         xlab = NULL,
                         ylab = NULL,
                         label_by = c("pair", "promoter", "compound"),
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
  required <- c(effect, y_col, promoter, compound, compound_label)
  missing_cols <- setdiff(required, names(table))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  d <- table
  d$.effect <- as.numeric(d[[effect]])
  d$.p_for_plot <- as.numeric(d[[y_col]])
  d$.promoter <- as.character(d[[promoter]])
  d$.compound <- as.character(d[[compound]])
  d$.compound_label <- as.character(d[[compound_label]])
  missing_label <- is.na(d$.compound_label) | !nzchar(d$.compound_label)
  d$.compound_label[missing_label] <- d$.compound[missing_label]
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
    by = list(promoter = d$.promoter),
    FUN = sum
  )
  names(hit_counts)[2] <- "hit_n"
  total_counts <- as.data.frame(table(d$.promoter), stringsAsFactors = FALSE)
  names(total_counts) <- c("promoter", "total_n")
  promoter_counts <- merge(hit_counts, total_counts, by = "promoter", all = TRUE)
  promoter_counts$hit_n[is.na(promoter_counts$hit_n)] <- 0
  promoter_counts <- promoter_counts[order(-promoter_counts$hit_n, -promoter_counts$total_n, promoter_counts$promoter), ]
  colored_promoters <- utils::head(promoter_counts$promoter, top_promoters)
  d$.promoter_group <- ifelse(d$.promoter %in% colored_promoters, d$.promoter, "Other promoters")
  d$.promoter_group <- factor(d$.promoter_group, levels = c(colored_promoters, "Other promoters"))

  label_df <- d[d$.hit, , drop = FALSE]
  label_df <- label_df[order(label_df$.p_for_plot, -abs(label_df$.effect)), , drop = FALSE]
  label_df <- utils::head(label_df, top_n)
  if (label_by == "pair") {
    label_df$.label <- paste(label_df$.promoter, label_df$.compound_label, sep = " / ")
  } else if (label_by == "promoter") {
    label_df$.label <- label_df$.promoter
  } else {
    label_df$.label <- label_df$.compound_label
  }
  if (is.finite(max_label_chars)) {
    too_long <- nchar(label_df$.label) > max_label_chars
    label_df$.label[too_long] <- paste0(
      substr(label_df$.label[too_long], 1, max_label_chars - 3),
      "..."
    )
  }

  palette <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9")
  values <- stats::setNames(rep(palette, length.out = length(colored_promoters)), colored_promoters)
  values <- c(values, "Other promoters" = "#C5C5C5")
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
      ggplot2::aes(color = .promoter_group, shape = .direction, alpha = .direction),
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
      color = "Promoter",
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

make_response_matrix <- function(table, value, promoter, compound_display) {
  d <- table[, c(promoter, compound_display, value), drop = FALSE]
  names(d) <- c(".promoter", ".compound_display", ".value")
  d$.value <- as.numeric(d$.value)
  d <- d[is.finite(d$.value), , drop = FALSE]
  if (nrow(d) == 0) {
    stop("No finite response values available for plotting.", call. = FALSE)
  }
  d <- stats::aggregate(
    .value ~ .promoter + .compound_display,
    d,
    mean,
    na.rm = TRUE
  )
  wide <- stats::reshape(
    d,
    idvar = ".promoter",
    timevar = ".compound_display",
    direction = "wide"
  )
  names(wide) <- sub("^\\.value\\.", "", names(wide))
  rownames(wide) <- wide$.promoter
  wide$.promoter <- NULL
  as.matrix(wide)
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

#' Heatmap of a DStressR promoter-by-compound response matrix
#'
#' Creates a standard heatmap for normalized promoter-compound responses. The
#' default `value` is `specific_effect`, matching [results()], but workflow
#' tables can use columns such as `destress_eb_effect_centered`.
#'
#' @param table A data frame with one row per promoter-compound pair.
#' @param value Numeric response/effect column to show in the heatmap.
#' @param promoter,compound Columns identifying promoters and compounds.
#' @param compound_label Optional human-readable compound-name column. Defaults
#'   to `compound`.
#' @param show_compound_ids If `TRUE`, append compound IDs in square brackets
#'   to compound labels.
#' @param top_n_compounds If finite, show only the top compounds by mean
#'   absolute response. Use `Inf` to show all compounds.
#' @param cluster_rows,cluster_cols If `TRUE`, hierarchically cluster promoters
#'   and/or compounds.
#' @param clip_quantile Quantile of absolute response values used to clip the
#'   color scale. Set to `1` to use the observed maximum.
#' @param title,subtitle,xlab,ylab Plot labels.
#' @param low,mid,high Colors for negative, zero, and positive responses.
#' @return A `ggplot` object. The plotted matrix is available as
#'   `attr(plot, "response_matrix")`.
#' @export
plot_response_heatmap <- function(table,
                                  value = "specific_effect",
                                  promoter = "promoter",
                                  compound = "compound",
                                  compound_label = compound,
                                  show_compound_ids = TRUE,
                                  top_n_compounds = 160,
                                  cluster_rows = TRUE,
                                  cluster_cols = TRUE,
                                  clip_quantile = 0.98,
                                  title = "DStressR promoter-by-compound matrix",
                                  subtitle = NULL,
                                  xlab = "Compounds",
                                  ylab = "Promoters",
                                  low = "#2166AC",
                                  mid = "white",
                                  high = "#B2182B") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package `ggplot2` is required for plot_response_heatmap().", call. = FALSE)
  }
  stopifnot(is.data.frame(table))
  required <- c(value, promoter, compound, compound_label)
  missing_cols <- setdiff(required, names(table))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  d <- table
  d$.promoter <- as.character(d[[promoter]])
  d$.compound <- as.character(d[[compound]])
  d$.compound_label <- as.character(d[[compound_label]])
  missing_label <- is.na(d$.compound_label) | !nzchar(d$.compound_label)
  d$.compound_label[missing_label] <- d$.compound[missing_label]
  d$.compound_display <- d$.compound_label
  if (isTRUE(show_compound_ids)) {
    d$.compound_display <- paste0(d$.compound_label, " [", d$.compound, "]")
  }
  d$.value <- as.numeric(d[[value]])
  d <- d[is.finite(d$.value), , drop = FALSE]
  if (nrow(d) == 0) {
    stop("No finite response values available for plotting.", call. = FALSE)
  }

  if (is.finite(top_n_compounds)) {
    compound_summary <- stats::aggregate(
      abs(.value) ~ .compound_display,
      d,
      mean,
      na.rm = TRUE
    )
    names(compound_summary)[2] <- ".mean_abs_value"
    compound_summary <- compound_summary[order(-compound_summary$.mean_abs_value), , drop = FALSE]
    keep_compounds <- utils::head(compound_summary$.compound_display, top_n_compounds)
    d <- d[d$.compound_display %in% keep_compounds, , drop = FALSE]
  }

  mat <- make_response_matrix(d, ".value", ".promoter", ".compound_display")
  if (isTRUE(cluster_rows)) {
    mat <- mat[matrix_cluster_order(mat, 1), , drop = FALSE]
  }
  if (isTRUE(cluster_cols)) {
    mat <- mat[, matrix_cluster_order(mat, 2), drop = FALSE]
  }

  plot_df <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
  names(plot_df) <- c(".promoter", ".compound_display", ".value")
  plot_df$.promoter <- factor(plot_df$.promoter, levels = rownames(mat))
  plot_df$.compound_display <- factor(plot_df$.compound_display, levels = colnames(mat))

  limit <- stats::quantile(abs(plot_df$.value), clip_quantile, na.rm = TRUE)
  if (!is.finite(limit) || limit <= 0) {
    limit <- max(abs(plot_df$.value), na.rm = TRUE)
  }
  plot_df$.plot_value <- pmax(pmin(plot_df$.value, limit), -limit)
  if (is.null(subtitle)) {
    subtitle <- if (is.finite(top_n_compounds)) {
      paste0("Top ", top_n_compounds, " compounds by mean absolute ", value)
    } else {
      "All compounds"
    }
  }

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(.compound_display, .promoter, fill = .plot_value)
  ) +
    ggplot2::geom_raster() +
    ggplot2::scale_fill_gradient2(
      low = low,
      mid = mid,
      high = high,
      midpoint = 0,
      limits = c(-limit, limit),
      name = value
    ) +
    ggplot2::theme_light(base_size = 8) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
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

  attr(p, "response_matrix") <- mat
  attr(p, "color_limit") <- limit
  p
}

#' Histogram of DStressR promoter-compound effects
#'
#' Shows the empirical distribution of normalized promoter-compound effects,
#' either over all matrix entries or faceted by promoter.
#'
#' @param table A data frame with one row per promoter-compound pair.
#' @param value Numeric effect column to plot.
#' @param promoter Column identifying promoters, used when `by = "promoter"`.
#' @param by Plot one pooled histogram (`"all"`) or promoter-faceted
#'   histograms (`"promoter"`).
#' @param bins Number of histogram bins.
#' @param xlim Optional two-element x-axis limit.
#' @param scales Facet scale behavior for `by = "promoter"`.
#' @param title,subtitle,xlab,ylab Plot labels.
#' @param fill,border Histogram fill and border colors.
#' @return A `ggplot` object.
#' @export
plot_effect_histogram <- function(table,
                                  value = "specific_effect",
                                  promoter = "promoter",
                                  by = c("all", "promoter"),
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
  if (by == "promoter") {
    required <- c(required, promoter)
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
  if (by == "promoter") {
    d$.promoter <- as.character(d[[promoter]])
  }
  if (is.null(xlab)) {
    xlab <- value
  }
  if (is.null(title)) {
    title <- if (by == "promoter") {
      "Effect distributions by promoter"
    } else {
      "Effect distribution over all promoter-compound entries"
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
  if (by == "promoter") {
    p <- p +
      ggplot2::facet_wrap(ggplot2::vars(.promoter), scales = scales, ncol = 5) +
      ggplot2::theme(strip.text = ggplot2::element_text(size = 8))
  }
  p
}

#' Clustered block map of a DStressR promoter-by-compound response matrix
#'
#' Hierarchically clusters promoters and compounds, cuts the dendrograms into
#' interpretable groups, and plots the mean response for each promoter-cluster by
#' compound-cluster block. This is useful as a compact overview when the full
#' compound library is too large for individual compound labels.
#'
#' @param table A data frame with one row per promoter-compound pair.
#' @param value Numeric response/effect column to summarize.
#' @param promoter,compound Columns identifying promoters and compounds.
#' @param compound_label Optional human-readable compound-name column. Defaults
#'   to `compound`.
#' @param show_compound_ids If `TRUE`, append compound IDs in square brackets
#'   to compound labels before clustering.
#' @param n_promoter_clusters,n_compound_clusters Number of dendrogram clusters
#'   to use for promoters and compounds.
#' @param missing_value Value used only for clustering missing matrix entries.
#'   Block summaries are still computed from observed finite values.
#' @param clip_quantile Quantile of absolute block means used to clip the color
#'   scale. Set to `1` to use the observed maximum.
#' @param show_counts If `TRUE`, annotate each tile with the number of compounds
#'   in that compound cluster.
#' @param title,subtitle,xlab,ylab Plot labels.
#' @param low,mid,high Colors for negative, zero, and positive responses.
#' @return A `ggplot` object with attributes `response_matrix`,
#'   `promoter_clusters`, `compound_clusters`, `block_summary`, `row_hclust`,
#'   and `col_hclust`.
#' @export
plot_response_cluster_blocks <- function(table,
                                         value = "specific_effect",
                                         promoter = "promoter",
                                         compound = "compound",
                                         compound_label = compound,
                                         show_compound_ids = TRUE,
                                         n_promoter_clusters = 6,
                                         n_compound_clusters = 14,
                                         missing_value = 0,
                                         clip_quantile = 0.98,
                                         show_counts = TRUE,
                                         title = "DStressR clustered response map",
                                         subtitle = NULL,
                                         xlab = "Compound clusters",
                                         ylab = "Promoter clusters",
                                         low = "#2166AC",
                                         mid = "white",
                                         high = "#B2182B") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package `ggplot2` is required for plot_response_cluster_blocks().", call. = FALSE)
  }
  stopifnot(is.data.frame(table))
  required <- c(value, promoter, compound, compound_label)
  missing_cols <- setdiff(required, names(table))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  d <- table
  d$.promoter <- as.character(d[[promoter]])
  d$.compound <- as.character(d[[compound]])
  d$.compound_label <- as.character(d[[compound_label]])
  missing_label <- is.na(d$.compound_label) | !nzchar(d$.compound_label)
  d$.compound_label[missing_label] <- d$.compound[missing_label]
  d$.compound_display <- d$.compound_label
  if (isTRUE(show_compound_ids)) {
    d$.compound_display <- paste0(d$.compound_label, " [", d$.compound, "]")
  }
  d$.value <- as.numeric(d[[value]])
  d <- d[is.finite(d$.value), , drop = FALSE]
  if (nrow(d) == 0) {
    stop("No finite response values available for plotting.", call. = FALSE)
  }

  mat <- make_response_matrix(d, ".value", ".promoter", ".compound_display")
  if (nrow(mat) < 2 || ncol(mat) < 2) {
    stop("Clustered block plots require at least two promoters and two compounds.", call. = FALSE)
  }
  if (n_promoter_clusters < 1 || n_promoter_clusters > nrow(mat)) {
    stop("`n_promoter_clusters` must be between 1 and the number of promoters.", call. = FALSE)
  }
  if (n_compound_clusters < 1 || n_compound_clusters > ncol(mat)) {
    stop("`n_compound_clusters` must be between 1 and the number of compounds.", call. = FALSE)
  }

  cluster_mat <- mat
  cluster_mat[!is.finite(cluster_mat)] <- missing_value
  row_hc <- stats::hclust(stats::dist(cluster_mat))
  col_hc <- stats::hclust(stats::dist(t(cluster_mat)))
  row_cluster <- stats::cutree(row_hc, k = n_promoter_clusters)[rownames(mat)]
  col_cluster <- stats::cutree(col_hc, k = n_compound_clusters)[colnames(mat)]

  promoter_assignments <- data.frame(
    promoter = rownames(mat),
    promoter_cluster = paste0("P", row_cluster),
    dendrogram_order = match(seq_along(rownames(mat)), row_hc$order),
    stringsAsFactors = FALSE
  )
  promoter_assignments <- promoter_assignments[
    order(row_cluster, promoter_assignments$dendrogram_order),
    ,
    drop = FALSE
  ]

  compound_assignments <- data.frame(
    compound_display = colnames(mat),
    compound_cluster = paste0("C", col_cluster),
    dendrogram_order = match(seq_along(colnames(mat)), col_hc$order),
    mean_abs_effect = colMeans(abs(mat), na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  compound_assignments <- compound_assignments[
    order(col_cluster, compound_assignments$dendrogram_order),
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
        .promoter_cluster = paste0("P", pc),
        .compound_cluster = paste0("C", cc),
        n_promoters = sum(row_cluster == pc),
        n_compounds = sum(col_cluster == cc),
        mean_effect = mean(sub, na.rm = TRUE),
        median_effect = stats::median(sub, na.rm = TRUE),
        mean_abs_effect = mean(abs(sub), na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
  }
  block_df <- do.call(rbind, block_rows)
  block_df$.promoter_cluster <- factor(block_df$.promoter_cluster, levels = row_levels)
  block_df$.compound_cluster <- factor(block_df$.compound_cluster, levels = col_levels)
  block_df$.count_label <- paste0("n=", block_df$n_compounds)

  limit <- stats::quantile(abs(block_df$mean_effect), clip_quantile, na.rm = TRUE)
  if (!is.finite(limit) || limit <= 0) {
    limit <- max(abs(block_df$mean_effect), na.rm = TRUE)
  }
  block_df$.plot_value <- pmax(pmin(block_df$mean_effect, limit), -limit)
  if (is.null(subtitle)) {
    subtitle <- paste0(
      n_promoter_clusters, " promoter clusters x ", n_compound_clusters,
      " compound clusters"
    )
  }

  p <- ggplot2::ggplot(
    block_df,
    ggplot2::aes(.compound_cluster, .promoter_cluster, fill = .plot_value)
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
  attr(p, "promoter_clusters") <- promoter_assignments
  attr(p, "compound_clusters") <- compound_assignments
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

#' Clustered heatmap with promoter and compound dendrograms
#'
#' Draws a clustered promoter-by-compound response heatmap with hierarchical
#' trees on both axes. Unlike [plot_response_cluster_blocks()], this keeps the
#' individual matrix cells visible and uses the dendrograms to reveal structure
#' without collapsing the data into coarse blocks.
#'
#' @param table A data frame with one row per promoter-compound pair.
#' @param value Numeric response/effect column to show in the heatmap.
#' @param promoter,compound Columns identifying promoters and compounds.
#' @param compound_label Optional human-readable compound-name column. Defaults
#'   to `compound`.
#' @param show_compound_ids If `TRUE`, append compound IDs in square brackets
#'   to compound labels before clustering.
#' @param top_n_compounds If finite, show only the top compounds by mean
#'   absolute response. Use `Inf` to show all compounds.
#' @param n_promoter_clusters,n_compound_clusters Number of dendrogram clusters
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
                                            promoter = "promoter",
                                            compound = "compound",
                                            compound_label = compound,
                                            show_compound_ids = TRUE,
                                            top_n_compounds = 400,
                                            n_promoter_clusters = 6,
                                            n_compound_clusters = 14,
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
  required <- c(value, promoter, compound, compound_label)
  missing_cols <- setdiff(required, names(table))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  d <- table
  d$.promoter <- as.character(d[[promoter]])
  d$.compound <- as.character(d[[compound]])
  d$.compound_label <- as.character(d[[compound_label]])
  missing_label <- is.na(d$.compound_label) | !nzchar(d$.compound_label)
  d$.compound_label[missing_label] <- d$.compound[missing_label]
  d$.compound_display <- d$.compound_label
  if (isTRUE(show_compound_ids)) {
    d$.compound_display <- paste0(d$.compound_label, " [", d$.compound, "]")
  }
  d$.value <- as.numeric(d[[value]])
  d <- d[is.finite(d$.value), , drop = FALSE]
  if (nrow(d) == 0) {
    stop("No finite response values available for plotting.", call. = FALSE)
  }

  if (is.finite(top_n_compounds)) {
    compound_summary <- stats::aggregate(
      abs(.value) ~ .compound_display,
      d,
      mean,
      na.rm = TRUE
    )
    names(compound_summary)[2] <- ".mean_abs_value"
    compound_summary <- compound_summary[order(-compound_summary$.mean_abs_value), , drop = FALSE]
    keep_compounds <- utils::head(compound_summary$.compound_display, top_n_compounds)
    d <- d[d$.compound_display %in% keep_compounds, , drop = FALSE]
  }

  mat <- make_response_matrix(d, ".value", ".promoter", ".compound_display")
  if (nrow(mat) < 2 || ncol(mat) < 2) {
    stop("Clustered heatmaps require at least two promoters and two compounds.", call. = FALSE)
  }

  n_promoter_clusters <- min(max(1, n_promoter_clusters), nrow(mat))
  n_compound_clusters <- min(max(1, n_compound_clusters), ncol(mat))
  cluster_mat <- mat
  cluster_mat[!is.finite(cluster_mat)] <- missing_value
  row_hc <- stats::hclust(stats::dist(cluster_mat))
  col_hc <- stats::hclust(stats::dist(t(cluster_mat)))
  row_clusters <- stats::cutree(row_hc, k = n_promoter_clusters)[rownames(mat)]
  col_clusters <- stats::cutree(col_hc, k = n_compound_clusters)[colnames(mat)]
  ordered_mat <- mat[row_hc$order, col_hc$order, drop = FALSE]

  color_limit <- stats::quantile(abs(mat), clip_quantile, na.rm = TRUE)
  if (!is.finite(color_limit) || color_limit <= 0) {
    color_limit <- max(abs(mat), na.rm = TRUE)
  }
  if (is.null(subtitle)) {
    subtitle <- if (is.finite(top_n_compounds)) {
      paste0("Top ", top_n_compounds, " compounds by mean absolute ", value)
    } else {
      "All compounds"
    }
  }

  promoter_assignments <- data.frame(
    promoter = rownames(mat),
    promoter_cluster = paste0("P", row_clusters),
    dendrogram_order = match(seq_along(rownames(mat)), row_hc$order),
    stringsAsFactors = FALSE
  )
  promoter_assignments <- promoter_assignments[
    order(row_clusters, promoter_assignments$dendrogram_order),
    ,
    drop = FALSE
  ]
  compound_assignments <- data.frame(
    compound_display = colnames(mat),
    compound_cluster = paste0("C", col_clusters),
    dendrogram_order = match(seq_along(colnames(mat)), col_hc$order),
    mean_abs_effect = colMeans(abs(mat), na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  compound_assignments <- compound_assignments[
    order(col_clusters, compound_assignments$dendrogram_order),
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
    promoter_clusters = promoter_assignments,
    compound_clusters = compound_assignments,
    color_limit = color_limit
  ))
}

utils::globalVariables(c(
  ".compound_cluster",
  ".compound_display",
  ".direction",
  ".effect",
  ".effect_hist",
  ".label",
  ".neg_log10_p",
  ".plot_value",
  ".promoter",
  ".promoter_cluster",
  ".promoter_group"
))
