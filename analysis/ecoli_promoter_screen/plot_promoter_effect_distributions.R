source(file.path("analysis", "_helpers.R"))

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package `ggplot2` is required for promoter effect distribution plots.", call. = FALSE)
}
ggplot2 <- asNamespace("ggplot2")

out_dir <- analysis_output_dir("binsfeld_promoter_diagnostics")
result_file <- file.path(
  analysis_project_root(),
  "analysis", "outputs", "binsfeld_evc_calibrated", "evc_huber_pair_results.tsv"
)
if (!file.exists(result_file)) {
  stop("Missing Binsfeld EVC result table: ", result_file, call. = FALSE)
}

res <- read_tsv_base(result_file)
required <- c(
  "promoter", "compound", "specific_effect", "specific_pvalue",
  "specific_padj_by_reporter", "evc_huber_hit_class"
)
missing_cols <- setdiff(required, names(res))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
}

res$specific_effect <- as.numeric(res$specific_effect)
res$specific_pvalue <- as.numeric(res$specific_pvalue)
res$specific_padj_by_reporter <- as.numeric(res$specific_padj_by_reporter)
res <- res[is.finite(res$specific_effect), , drop = FALSE]
res$hit <- is.finite(res$specific_padj_by_reporter) & res$specific_padj_by_reporter <= 0.05
res$effect_class <- ifelse(
  res$hit & res$specific_effect > 0,
  "Positive hit",
  ifelse(res$hit & res$specific_effect < 0, "Negative hit", "Not significant")
)

promoter_summary <- do.call(rbind, lapply(split(res, res$promoter), function(tab) {
  data.frame(
    promoter = tab$promoter[1],
    n_perturbations = length(unique(tab$compound)),
    median_effect = stats::median(tab$specific_effect, na.rm = TRUE),
    mean_effect = mean(tab$specific_effect, na.rm = TRUE),
    median_abs_effect = stats::median(abs(tab$specific_effect), na.rm = TRUE),
    mean_abs_effect = mean(abs(tab$specific_effect), na.rm = TRUE),
    effect_iqr = stats::IQR(tab$specific_effect, na.rm = TRUE),
    effect_sd = stats::sd(tab$specific_effect, na.rm = TRUE),
    n_hits = sum(tab$hit, na.rm = TRUE),
    n_positive_hits = sum(tab$hit & tab$specific_effect > 0, na.rm = TRUE),
    n_negative_hits = sum(tab$hit & tab$specific_effect < 0, na.rm = TRUE),
    min_pvalue = min(tab$specific_pvalue, na.rm = TRUE),
    min_padj_by_reporter = min(tab$specific_padj_by_reporter, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))
paper_reporter_order <- c("acrABp", "marRABp", "micFp", "ompFp", "robp", "soxSp", "tolCp")
promoter_summary$paper_order <- match(promoter_summary$promoter, paper_reporter_order)
promoter_summary <- promoter_summary[
  order(promoter_summary$paper_order, promoter_summary$promoter),
  ,
  drop = FALSE
]
promoter_summary$rank_median_abs_effect <- rank(
  promoter_summary$median_abs_effect,
  ties.method = "first"
)

utils::write.table(
  promoter_summary,
  file.path(out_dir, "binsfeld_promoter_effect_distribution_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

promoter_levels <- rev(paper_reporter_order[paper_reporter_order %in% promoter_summary$promoter])
res$promoter <- factor(res$promoter, levels = promoter_levels)
res$effect_class <- factor(
  res$effect_class,
  levels = c("Negative hit", "Not significant", "Positive hit")
)

effect_limit <- stats::quantile(abs(res$specific_effect), 0.995, na.rm = TRUE)
effect_limit <- max(effect_limit, max(abs(res$specific_effect[res$hit]), na.rm = TRUE))
effect_limit <- ceiling(effect_limit * 10) / 10

distribution_plot <- ggplot2$ggplot(res, ggplot2$aes(specific_effect, promoter)) +
  ggplot2$geom_violin(
    ggplot2$aes(group = promoter),
    fill = "#E5E7EB",
    color = "#9CA3AF",
    linewidth = 0.25,
    scale = "width",
    trim = TRUE
  ) +
  ggplot2$geom_boxplot(
    width = 0.12,
    outlier.shape = NA,
    fill = "white",
    color = "#374151",
    linewidth = 0.25
  ) +
  ggplot2$geom_point(
    data = res[res$effect_class == "Not significant", , drop = FALSE],
    color = "#9CA3AF",
    alpha = 0.28,
    size = 0.85,
    position = ggplot2$position_jitter(height = 0.075, width = 0, seed = 3)
  ) +
  ggplot2$geom_point(
    data = res[res$effect_class != "Not significant", , drop = FALSE],
    ggplot2$aes(color = effect_class),
    alpha = 0.9,
    size = 1.45,
    position = ggplot2$position_jitter(height = 0.075, width = 0, seed = 4)
  ) +
  ggplot2$geom_vline(xintercept = 0, color = "#111827", linewidth = 0.3) +
  ggplot2$scale_color_manual(
    values = c("Negative hit" = "#1D4ED8", "Positive hit" = "#B91C1C"),
    drop = TRUE,
    name = NULL
  ) +
  ggplot2$coord_cartesian(xlim = c(-effect_limit, effect_limit)) +
  ggplot2$theme_light(base_size = 10) +
  ggplot2$theme(
    panel.grid.minor = ggplot2$element_blank(),
    panel.grid.major.y = ggplot2$element_blank(),
    axis.title = ggplot2$element_text(color = "#111827"),
    axis.text = ggplot2$element_text(color = "#4B5563"),
    legend.position = c(0.985, 0.985),
    legend.justification = c(1, 1),
    legend.background = ggplot2$element_rect(
      fill = grDevices::adjustcolor("white", alpha.f = 0.9),
      color = "#D1D5DB",
      linewidth = 0.25
    ),
    legend.key = ggplot2$element_blank(),
    legend.margin = ggplot2$margin(3, 5, 3, 5)
  ) +
  ggplot2$labs(
    x = expression(hat(Delta * y)[a * "\u00b7"]^"spec"),
    y = NULL
  )

ggplot2$ggsave(
  file.path(out_dir, "binsfeld_promoter_specific_effect_distributions.png"),
  distribution_plot,
  width = 6.6,
  height = 3.7,
  dpi = 320,
  bg = "white"
)
ggplot2$ggsave(
  file.path(out_dir, "binsfeld_promoter_specific_effect_distributions.pdf"),
  distribution_plot,
  width = 6.6,
  height = 3.7
)

message("Wrote Binsfeld promoter effect distribution diagnostics to: ", out_dir)
