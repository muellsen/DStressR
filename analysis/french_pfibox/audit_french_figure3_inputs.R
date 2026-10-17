# Trace a pre-specified high-membership subset without fitting or reclustering.
source(file.path("analysis", "_helpers.R"))
out <- analysis_output_dir("french_pfibox", "figure3_input_audit")
cache <- file.path(out, "selected_source_trajectories.rds")
write_tsv <- function(x, name) utils::write.table(x, file.path(out, name),
  sep = "\t", quote = FALSE, row.names = FALSE)
if (!file.exists(cache)) {
  tab <- suppressMessages(readxl::read_excel(analysis_path("data-raw", "french_pfibox", "source",
    "french_pfibox_table_s3_clusterdata_mmc2.xlsx"), sheet = "ClusterData", col_names = FALSE,
    .name_repair = "minimal"))
  major <- suppressWarnings(as.integer(tab[[2]]))
  keep <- !is.na(tab[[1]]) & !is.na(major) & major %in% 1:12
  annotation <- data.frame(reporter = as.character(tab[[1]][keep]), cluster = major[keep])
  membership <- sapply(tab[keep, 3:14], as.numeric)
  annotation$membership <- membership[cbind(seq_len(nrow(annotation)), annotation$cluster)]
  selected <- do.call(rbind, lapply(1:12, function(k) {
    x <- annotation[annotation$cluster == k, ]
    x <- x[order(-x$membership, x$reporter), ]
    x <- head(x, 2)
    x$selection_rank <- seq_len(nrow(x))
    x
  }))
  message("Loading source once for 24 pre-specified promoters.")
  e <- new.env(parent = emptyenv())
  load(analysis_path("data-raw", "french_pfibox", "source", "gfpData.Rdata"), envir = e)
  legend <- e$gfpData$legend
  positions <- which(as.character(legend$Promoter) %in% selected$reporter)
  ids <- cbind(position = positions, as.data.frame(legend[positions, ]))
  ids <- merge(ids, selected, by.x = "Promoter", by.y = "reporter", sort = FALSE)
  ids <- ids[order(ids$position), ]
  controls <- grep("-0$", names(e$gfpData), value = TRUE)
  traces <- lapply(controls, function(t) {
    x <- e$gfpData[[t]]
    lapply(x[c("raw", "raw.loess", "constit.normalized", "promoter.activity")],
      function(m) as.matrix(m)[ids$position, , drop = FALSE])
  })
  names(traces) <- controls
  saveRDS(list(selected = selected, ids = ids, traces = traces), cache)
  rm(e)
  invisible(gc())
}
d <- readRDS(cache)
write_tsv(d$selected, "preselected_promoters.tsv")
write_tsv(d$ids, "preselected_colony_positions.tsv")
scale01 <- function(x) {
  good <- is.finite(x)
  if (sum(good) < 2 || diff(range(x[good])) == 0) return(rep(NA_real_, length(x)))
  (x - min(x[good])) / diff(range(x[good]))
}
comparisons <- list()
ranges <- list()
for (t in names(d$traces)) {
  message("Checking derivative alignment: ", t)
  x <- d$traces[[t]]
  for (i in seq_len(nrow(d$ids))) {
    pa <- x$promoter.activity[i, ]
    if (!all(is.finite(pa))) next
    q <- x$constit.normalized[i, ]
    dq <- diff(q)
    xx <- seq_along(dq)
    candidates <- list(unsmoothed_difference = dq)
    for (span in c(0.1, 0.2, 0.5)) {
      good <- is.finite(dq)
      if (sum(good) < 30) next
      candidates[[paste0("loess_span_", span)]] <- suppressWarnings(as.numeric(stats::predict(
        stats::loess(y ~ time, data = data.frame(y = dq[good], time = xx[good]), span = span),
        newdata = data.frame(time = xx))))
    }
    ranges[[length(ranges) + 1L]] <- data.frame(control = t, reporter = d$ids$Promoter[i],
      position = d$ids$position[i], pa_min = min(pa), pa_max = max(pa),
      normalized_fluorescence_finite = sum(is.finite(q)), raw_samples = length(q), activity_samples = length(pa))
    for (method in names(candidates)) {
      v <- candidates[[method]]
      for (offset in 0:(length(v) - length(pa))) {
        candidate <- v[seq_along(pa) + offset]
        good <- is.finite(candidate) & is.finite(pa)
        if (sum(good) < 0.9 * length(pa)) next
        scaled <- scale01(candidate)
        comparisons[[length(comparisons) + 1L]] <- data.frame(
          control = t, reporter = d$ids$Promoter[i], position = d$ids$position[i], method = method,
          derivative_prefix_skipped = offset, n_compared = sum(good),
          correlation = suppressWarnings(stats::cor(candidate[good], pa[good])),
          scaled_rmse = sqrt(mean((scaled[good] - pa[good])^2)))
      }
    }
  }
}
comparisons <- do.call(rbind, comparisons)
comparisons <- comparisons[is.finite(comparisons$correlation), ]
write_tsv(comparisons, "derivative_alignment_candidates.tsv")
key <- paste(comparisons$control, comparisons$position, sep = "::")
comparisons <- comparisons[order(key, -comparisons$correlation, comparisons$scaled_rmse), ]
best <- comparisons[!duplicated(paste(comparisons$control, comparisons$position, sep = "::")), ]
write_tsv(best, "best_numerical_alignments.tsv")
write_tsv(do.call(rbind, ranges), "selected_activity_ranges.tsv")
summary <- data.frame(complete_colony_trajectories = nrow(best),
  median_best_correlation = median(best$correlation), minimum_best_correlation = min(best$correlation),
  fraction_correlation_above_099 = mean(best$correlation > 0.99),
  median_best_scaled_rmse = median(best$scaled_rmse))
print(summary)
print(as.data.frame(table(best$method, best$derivative_prefix_skipped)))
write_tsv(summary, "numerical_alignment_summary.tsv")

# Show raw stages for fixed examples: the highest-membership promoter in F2/F4/F7/F12.
examples <- d$selected[d$selected$selection_rank == 1 & d$selected$cluster %in% c(2, 4, 7, 12), ]
control <- "Ampicillin-0"
x <- d$traces[[control]]
long <- list()
for (a in examples$reporter) {
  ii <- which(d$ids$Promoter == a)
  for (stage in names(x)) for (i in ii) {
    v <- x[[stage]][i, ]
    if (!any(is.finite(v))) next
    long[[length(long) + 1L]] <- data.frame(reporter = a, stage = stage, sample = seq_along(v),
      value = v, colony = paste0("Position ", d$ids$position[i]))
  }
}
long <- do.call(rbind, long)
long$stage <- factor(long$stage, levels = c("raw", "raw.loess", "constit.normalized", "promoter.activity"),
  labels = c("Fluorescence", "Smoothed fluorescence", "Biomass-normalized fluorescence", "Deposited activity (0-1)"))
row_labels <- paste0(examples$reporter, " | French F", examples$cluster,
  " | membership ", sprintf("%.3f", examples$membership))
long$reporter <- factor(long$reporter, levels = examples$reporter, labels = row_labels)
p <- ggplot2::ggplot(long, ggplot2::aes(sample, value, group = colony)) +
  ggplot2::geom_line(color = "#287E9E", alpha = 0.65, linewidth = 0.45) +
  ggplot2::facet_wrap(reporter ~ stage, ncol = 4, scales = "free") +
  ggplot2::labs(x = "Sample index within each deposited stage (not aligned absolute time)", y = NULL) +
  ggplot2::theme_classic(base_size = 11) + ggplot2::theme(
    strip.background = ggplot2::element_blank(), strip.text = ggplot2::element_text(size = 9),
    panel.spacing = grid::unit(10, "pt"))
figure <- ggplot2::ggplotGrob(p)
ggplot2::ggsave(file.path(out, "selected_promoter_processing_stages.png"), figure,
  width = 18, height = 11, dpi = 200, bg = "white")

# Preserve individual colonies; do not pool controls into a claimed reconstruction.
strongest <- d$selected[d$selected$selection_rank == 1, ]
long <- list()
for (t in names(d$traces)) for (i in which(d$ids$Promoter %in% strongest$reporter)) {
  v <- d$traces[[t]]$promoter.activity[i, ]
  if (!all(is.finite(v)) || stats::sd(v) == 0) next
  a <- d$ids$Promoter[i]
  long[[length(long) + 1L]] <- data.frame(reporter = a, control = t,
    colony = paste(t, d$ids$position[i], sep = "::"), sample = seq_along(v),
    value = (v - mean(v)) / stats::sd(v))
}
long <- do.call(rbind, long)
long$panel <- factor(long$reporter, levels = strongest$reporter,
  labels = paste0("F", strongest$cluster, ": ", strongest$reporter,
    " (membership ", sprintf("%.3f", strongest$membership), ")"))
p <- ggplot2::ggplot(long, ggplot2::aes(sample, value, group = colony, color = control)) +
  ggplot2::geom_line(alpha = 0.35, linewidth = 0.35) +
  ggplot2::facet_wrap(~panel, ncol = 3, scales = "free_y") +
  ggplot2::scale_color_manual(values = grDevices::hcl.colors(length(d$traces), "Dark 3"), name = "Untreated source matrix") +
  ggplot2::labs(x = "Deposited activity sample", y = "Standardized activity (individual colonies)") +
  ggplot2::theme_classic(base_size = 11) + ggplot2::theme(strip.background = ggplot2::element_blank(),
    legend.position = "bottom", legend.text = ggplot2::element_text(size = 8)) +
  ggplot2::guides(color = ggplot2::guide_legend(nrow = 2))
figure <- ggplot2::ggplotGrob(p)
ggplot2::ggsave(file.path(out, "strongest_members_individual_untreated_colonies.png"), figure,
  width = 15, height = 12, dpi = 200, bg = "white")
