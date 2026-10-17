# Numerical processing audit only: no effect models or clustering.
out <- file.path("analysis", "outputs", "french_pfibox", "figure3_input_audit")
d <- readRDS(file.path(out, "selected_source_trajectories.rds"))
scale01 <- function(x) {
  if (!all(is.finite(x)) || diff(range(x)) == 0) return(rep(NA_real_, length(x)))
  (x - min(x)) / diff(range(x))
}
smooth <- function(y, span, degree) {
  time <- seq_along(y)
  good <- is.finite(y)
  if (sum(good) < 30) return(rep(NA_real_, length(y)))
  suppressWarnings(tryCatch(as.numeric(predict(loess(y ~ time,
    data = data.frame(y = y[good], time = time[good]), span = span,
    degree = degree), newdata = data.frame(time = time))),
    error = function(e) rep(NA_real_, length(y))))
}
rules <- expand.grid(offset = 0:19, span = 0.5, degree = 2,
  order = "trim_then_smooth", stringsAsFactors = FALSE)
rules <- unique(rbind(rules, expand.grid(offset = 19, span = c(.1, .2, .3, .5, .75, 1),
  degree = c(1, 2), order = c("trim_then_smooth", "smooth_then_trim"),
  stringsAsFactors = FALSE)))
rules$rule <- sprintf("%s_offset%d_span%g_degree%d", rules$order, rules$offset,
  rules$span, rules$degree)
results <- list()
curves <- list()
if (identical(Sys.getenv("PFIBOX_AUDIT_PLOT_ONLY"), "1")) {
  results <- read.delim(file.path(out, "processing_order_candidates.tsv"))
  curves <- readRDS(file.path(out, "trimmed_derivative_activity_candidates.rds"))
} else {
for (control in names(d$traces)) {
  message("Reconstruction audit: ", control)
  x <- d$traces[[control]]
  for (i in seq_len(nrow(d$ids))) {
    pa <- x$promoter.activity[i, ]
    if (!all(is.finite(pa))) next
    dq <- diff(x$constit.normalized[i, ])
    for (k in seq_len(nrow(rules))) {
      z <- rules[k, ]
      idx <- seq_along(pa) + z$offset
      v <- if (z$order == "trim_then_smooth") smooth(dq[idx], z$span, z$degree) else
        smooth(dq, z$span, z$degree)[idx]
      scaled <- scale01(v)
      complete <- all(is.finite(scaled))
      results[[length(results) + 1L]] <- data.frame(control = control,
        reporter = d$ids$Promoter[i], position = d$ids$position[i], rule = z$rule,
        complete = complete, correlation = if (complete) cor(scaled, pa) else NA_real_,
        rmse = if (complete) sqrt(mean((scaled - pa)^2)) else NA_real_,
        max_error = if (complete) max(abs(scaled - pa)) else NA_real_)
      if (z$order == "trim_then_smooth" && z$offset == 19 && z$span == .5 && z$degree == 2)
        curves[[paste(control, d$ids$position[i], sep = "::")]] <- list(
          control = control, reporter = d$ids$Promoter[i], position = d$ids$position[i],
          unscaled_activity = v, scaled_activity = scaled, deposited_activity = pa)
    }
  }
}
results <- do.call(rbind, results)
}
summary <- do.call(rbind, lapply(split(results, results$rule), function(x) {
  good <- x$complete
  data.frame(rule = x$rule[1], n_eligible = nrow(x), n_complete = sum(good),
    median_correlation = median(x$correlation[good]), median_rmse = median(x$rmse[good]),
    maximum_rmse = max(x$rmse[good]), fraction_exact_1e8 = mean(x$max_error[good] < 1e-8))
}))
summary <- summary[order(-summary$n_complete, summary$median_rmse), ]
write.table(results, file.path(out, "processing_order_candidates.tsv"), sep = "\t",
  quote = FALSE, row.names = FALSE)
write.table(summary, file.path(out, "processing_order_summary.tsv"), sep = "\t",
  quote = FALSE, row.names = FALSE)
saveRDS(curves, file.path(out, "trimmed_derivative_activity_candidates.rds"))
print(head(summary, 12), row.names = FALSE)

# Fixed illustrative promoters; all eligible trajectories remain in the cache.
examples <- d$selected[d$selected$selection_rank == 1 & d$selected$cluster %in% c(2, 4, 7, 12), ]
plot_data <- do.call(rbind, lapply(examples$reporter, function(a) {
  eligible <- Filter(function(x) x$reporter == a && all(is.finite(x$scaled_activity)), curves)
  if (!length(eligible)) return(NULL)
  x <- eligible[[1]]
  data.frame(reporter = a, sample = seq_along(x$scaled_activity),
    deposited = x$deposited_activity, reconstructed = x$scaled_activity,
    unscaled = x$unscaled_activity, control = x$control, position = x$position)
}))
long <- rbind(transform(plot_data, value = deposited, series = "Deposited"),
  transform(plot_data, value = reconstructed, series = "Reconstructed"))
p <- ggplot2::ggplot(long, ggplot2::aes(sample, value, color = series, linetype = series)) +
  ggplot2::geom_line(linewidth = .65) + ggplot2::facet_wrap(~reporter, ncol = 2) +
  ggplot2::scale_color_manual(values = c(Deposited = "#222222", Reconstructed = "#DF7E32")) +
  ggplot2::scale_linetype_manual(values = c(Deposited = "solid", Reconstructed = "dashed")) +
  ggplot2::labs(x = "Deposited activity sample", y = "Min-max-scaled activity",
    color = NULL, linetype = NULL) + ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(legend.position = "bottom", strip.background = ggplot2::element_blank())
figure <- ggplot2::ggplotGrob(p)
ggplot2::ggsave(file.path(out, "activity_reconstruction_check.png"), figure,
  width = 9, height = 6, dpi = 200, bg = "white")
write.table(unique(plot_data[c("reporter", "control", "position")]),
  file.path(out, "activity_reconstruction_figure_sources.tsv"), sep = "\t",
  quote = FALSE, row.names = FALSE)
