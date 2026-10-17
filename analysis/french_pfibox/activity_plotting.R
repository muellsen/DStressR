# Display-only helper for already reconstructed activity trajectories.
# Column names are configurable; no smoothing, fitting or inference is performed.
plot_activity_trajectories <- function(data, time = "time_h", response = "value",
  condition = "condition", trajectory = "colony", panel = "panel",
  ncol = 2, scales = "free_y", summary = c("median", "mean", "none"),
  reference = "Untreated references",
  colors = c("Untreated references" = "#747474", Treated = "#DF7E32"),
  xlab = "Time (h)", ylab = "Unscaled promoter activity (per h)",
  base_size = 11, strip_size = 10, summary_linewidth = .95) {
  summary <- match.arg(summary)
  stopifnot(all(c(time, response, condition, trajectory, panel) %in% names(data)))
  d <- data.frame(.time = data[[time]], .value = data[[response]],
    .condition = data[[condition]], .panel = data[[panel]])
  d$.trajectory <- interaction(data[[panel]], data[[condition]], data[[trajectory]], drop = TRUE)
  d <- d[is.finite(d$.time) & is.finite(d$.value), ]
  levels <- if (is.factor(d$.condition)) levels(d$.condition) else unique(d$.condition)
  widths <- stats::setNames(rep(.9, length(levels)), levels)
  opacity <- stats::setNames(rep(.8, length(levels)), levels)
  if (reference %in% levels) {
    widths[reference] <- .3
    opacity[reference] <- .16
  }
  p <- ggplot2::ggplot(d, ggplot2::aes(.time, .value, group = .trajectory, color = .condition)) +
    ggplot2::geom_line(ggplot2::aes(linewidth = .condition, alpha = .condition)) +
    ggplot2::facet_wrap(~.panel, ncol = ncol, scales = scales) +
    ggplot2::scale_linewidth_manual(values = widths, guide = "none") +
    ggplot2::scale_alpha_manual(values = opacity, guide = "none") +
    ggplot2::labs(x = xlab, y = ylab, color = NULL) +
    ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = strip_size),
      legend.position = "bottom", panel.spacing = grid::unit(12, "pt"))
  if (!is.null(colors)) p <- p + ggplot2::scale_color_manual(values = colors, name = NULL)
  if (summary != "none") {
    fun <- if (summary == "median") stats::median else mean
    med <- stats::aggregate(.value ~ .panel + .condition + .time, d, fun)
    p <- p + ggplot2::geom_line(data = med, ggplot2::aes(group = .condition),
      linewidth = summary_linewidth)
  }
  p
}
