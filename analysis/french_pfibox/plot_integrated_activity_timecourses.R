# Prepare the retained 10-positive / 10-negative display selection and curves.
source(file.path("analysis", "_helpers.R"))
source(file.path("analysis", "french_pfibox", "activity_plotting.R"))
out <- analysis_output_dir("french_pfibox", "temporal_unscaled_integral_analysis", "timecourses")
tag <- "reporters2000_perturbations100_Tall_available_treatment_replicate_scalable_temporal_unscaled_integral"
res <- read.delim(file.path(dirname(out), paste0("pfibox_pa_destress_pair_results_", tag, ".tsv")))
eligible <- res[is.finite(res$total_pvalue) & res$total_padj_global < .05, ]
top <- do.call(rbind, lapply(c("positive", "negative"), function(direction) {
  x <- eligible[if (direction == "positive") eligible$total_effect > 0 else eligible$total_effect < 0, ]
  x$effect_percentile <- rank(abs(x$total_effect), ties.method = "average") / nrow(x)
  x$evidence_percentile <- rank(-log10(pmax(x$total_pvalue, .Machine$double.xmin)),
    ties.method = "average") / nrow(x)
  x$display_score <- sqrt(x$effect_percentile * x$evidence_percentile)
  x <- x[order(-x$display_score, x$total_pvalue, -abs(x$total_effect), x$reporter), ]
  selected <- integer()
  for (i in seq_len(nrow(x))) {
    drug <- sub("-[0-9]+$", "", x$perturbation[i])
    previous <- sub("-[0-9]+$", "", x$perturbation[selected])
    if (x$reporter[i] %in% x$reporter[selected] || sum(previous == drug) >= 2L) next
    selected <- c(selected, i)
    if (length(selected) == 10L) break
  }
  stopifnot(length(selected) == 10L)
  x <- x[selected, ]
  x$direction <- direction
  x$direction_rank <- seq_len(nrow(x))
  x
}))
top$rank <- seq_len(nrow(top))
write.table(top, file.path(out, "top10_each_direction_by_volcano_diverse_pairs.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE)

# Cache only the source trajectories needed by these pairs and all untreated
# references. The published source file itself remains untouched.
cache <- file.path(out, "selected_pair_source_trajectories.rds")
if (!file.exists(cache)) {
  message("Loading source once for the selected trajectory cache.")
  e <- new.env(parent = emptyenv())
  load(analysis_path("data-raw", "french_pfibox", "source", "gfpData.Rdata"), envir = e)
  legend <- as.data.frame(e$gfpData$legend)
  positions <- which(legend$Promoter %in% top$reporter)
  treatments <- unique(c(grep("-0$", names(e$gfpData), value = TRUE), top$perturbation))
  q <- lapply(treatments, function(t) as.matrix(e$gfpData[[t]]$constit.normalized)[positions, , drop = FALSE])
  names(q) <- treatments
  saveRDS(list(ids = cbind(position = positions, legend[positions, ]), selected = q), cache)
  rm(e)
  invisible(gc())
}
d <- readRDS(cache)
stopifnot(all(top$reporter %in% d$ids$Promoter), all(top$perturbation %in% names(d$selected)))
time <- seq_len(197L)
S <- vapply(time, function(k) {
  y <- as.numeric(time == k)
  as.numeric(predict(loess(y ~ time, span = .5, degree = 2)))
}, numeric(197L))
pa <- lapply(d$selected, function(q) {
  dq <- q[, 21:217, drop = FALSE] - q[, 20:216, drop = FALSE]
  z <- dq %*% t(S) * 12
  z[rowSums(is.finite(dq)) != 197L, ] <- NA_real_
  z
})
long <- list()
for (k in seq_len(nrow(top))) {
  positions <- which(d$ids$Promoter == top$reporter[k])
  sources <- c(grep("-0$", names(pa), value = TRUE), top$perturbation[k])
  for (s in sources) for (i in positions) {
    v <- pa[[s]][i, ]
    if (!all(is.finite(v))) next
    long[[length(long) + 1L]] <- data.frame(rank = k, time_h = (20:216) / 12,
      value = v, condition = if (s == top$perturbation[k]) "Treated" else "Untreated references",
      source = s, colony = paste(s, d$ids$position[i], sep = "::"))
  }
}
long <- do.call(rbind, long)
labels <- paste0(top$direction_rank, ". ", top$reporter, " | ", top$perturbation,
  "\nEffect = ", formatC(top$total_effect, format = "fg", digits = 4),
  "; q = ", formatC(top$total_padj_global, format = "e", digits = 1))
long$panel <- factor(long$rank, levels = top$rank, labels = labels)
long$condition <- factor(long$condition, levels = c("Untreated references", "Treated"))
med <- aggregate(value ~ panel + condition + time_h, long, median)
saveRDS(list(top_pairs = top, trajectories = long, medians = med),
  file.path(out, "top10_each_direction_by_volcano_diverse_timecourses.rds"))
for (direction in c("positive", "negative")) {
  panels <- labels[top$direction == direction]
  p <- plot_activity_trajectories(long[long$panel %in% panels, ],
    base_size = 12, summary_linewidth = .85)
  grob <- ggplot2::ggplotGrob(p)
  ggplot2::ggsave(file.path(out, paste0("top10_", direction, "_by_volcano_diverse_timecourses.png")),
    grob, width = 12, height = 14, dpi = 220, bg = "white")
}
message("Saved 20 display pairs and their colony trajectories in ", out)
