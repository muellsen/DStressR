source(file.path("analysis", "_helpers.R"))
load_destress_package()

required_packages <- c("ggplot2", "gridExtra")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Missing required package(s): ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

ggplot2 <- asNamespace("ggplot2")
gridExtra <- asNamespace("gridExtra")
out_dir <- analysis_output_dir("dryad_global_regulators")

default_path <- file.path(out_dir, "dryad_weak_stress_windows_calibrated_alpha_pair_results.tsv")
rank1_path <- file.path(out_dir, "dryad_weak_stress_windows_calibrated_alpha_rank1_pair_results.tsv")
if (!file.exists(default_path) || !file.exists(rank1_path)) {
  stop(
    "Missing calibrated-alpha result tables. Run the calibrated-alpha and rank-1 workflows first.",
    call. = FALSE
  )
}

default <- utils::read.delim(default_path, check.names = FALSE)
rank1 <- utils::read.delim(rank1_path, check.names = FALSE)

make_target <- function(d, target, effect, pvalue, padj) {
  out <- d[, c("promoter", "compound", "base_promoter", "window", effect, pvalue, padj), drop = FALSE]
  names(out)[5:7] <- c("effect", "pvalue", "padj")
  out$target <- target
  out$hit <- is.finite(out$padj) & out$padj < 0.05
  out
}

targets <- rbind(
  make_target(
    default,
    "Total effect",
    "total_effect",
    "total_pvalue",
    "total_padj_by_promoter"
  ),
  make_target(
    default,
    "Default specific effect",
    "specific_effect",
    "specific_pvalue",
    "specific_padj_by_promoter"
  ),
  make_target(
    rank1,
    "Rank-adjusted total effect",
    "rank_adjusted_total_effect",
    "rank_adjusted_total_pvalue",
    "rank_adjusted_total_padj_by_promoter"
  )
)
targets$target <- factor(
  targets$target,
  levels = c("Total effect", "Default specific effect", "Rank-adjusted total effect")
)
targets$compound <- factor(targets$compound, levels = c("Iron", "Tetracycline", "H2O2", "Kanamycin"))
targets$promoter_label <- paste(targets$base_promoter, targets$window, sep = " | ")
targets$neglog10_padj <- -log10(pmax(targets$padj, .Machine$double.xmin))
targets$neglog10_padj_plot <- pmin(targets$neglog10_padj, 35)
targets$expected_pair <- paste(targets$base_promoter, targets$compound, sep = "|") %in%
  c("Fur|Iron", "MarA|Tetracycline", "SoxS|H2O2", "LexA|Kanamycin")

target_summary <- do.call(rbind, lapply(split(targets, targets$target), function(d) {
  data.frame(
    target = as.character(unique(d$target)),
    tested_pairs = nrow(d),
    significant_pairs = sum(d$hit, na.rm = TRUE),
    median_pvalue = stats::median(d$pvalue, na.rm = TRUE),
    min_pvalue = min(d$pvalue, na.rm = TRUE),
    max_pvalue = max(d$pvalue, na.rm = TRUE),
    median_padj = stats::median(d$padj, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))
target_summary <- target_summary[match(levels(targets$target), target_summary$target), , drop = FALSE]
utils::write.table(
  target_summary,
  file.path(out_dir, "dryad_weak_stress_windows_calibrated_alpha_three_target_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
utils::write.table(
  targets,
  file.path(out_dir, "dryad_weak_stress_windows_calibrated_alpha_three_target_results.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

p_hist <- ggplot2::ggplot(targets, ggplot2::aes(x = pvalue, fill = target)) +
  ggplot2::geom_histogram(
    breaks = seq(0, 1, by = 0.05),
    color = "white",
    linewidth = 0.25
  ) +
  ggplot2::facet_wrap(ggplot2::vars(target), nrow = 1) +
  ggplot2::scale_fill_manual(
    values = c(
      "Total effect" = "#D55E00",
      "Default specific effect" = "#009E73",
      "Rank-adjusted total effect" = "#0072B2"
    ),
    guide = "none"
  ) +
  ggplot2::theme_light(base_size = 9) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    strip.background = ggplot2::element_rect(fill = "grey95", color = "grey75"),
    strip.text = ggplot2::element_text(face = "bold", size = 9, color = "grey15"),
    legend.position = "none"
  ) +
  ggplot2::labs(
    x = "Raw p-value",
    y = "Number of tests"
  )
ggplot2::ggsave(
  file.path(out_dir, "dryad_weak_stress_windows_calibrated_alpha_three_target_pvalue_histograms.png"),
  p_hist,
  width = 11,
  height = 3.2,
  dpi = 300
)
ggplot2::ggsave(
  file.path(out_dir, "dryad_weak_stress_windows_calibrated_alpha_three_target_pvalue_histograms.pdf"),
  p_hist,
  width = 11,
  height = 3.2
)

make_label_df <- function(d) {
  d$label_priority <- ifelse(d$expected_pair, 0, 1)
  candidates <- d[d$hit | d$expected_pair, , drop = FALSE]
  candidates <- candidates[
    order(candidates$target, candidates$label_priority, candidates$padj, -abs(candidates$effect)),
    ,
    drop = FALSE
  ]
  top_by_target <- c(
    "Total effect" = 8,
    "Default specific effect" = 13,
    "Rank-adjusted total effect" = 14
  )
  candidates <- do.call(rbind, lapply(split(candidates, candidates$target), function(x) {
    utils::head(x, top_by_target[as.character(x$target[1])])
  }))
  candidates$label <- paste(candidates$base_promoter, candidates$window, candidates$compound)
  candidates
}

label_df <- make_label_df(targets)

p_volcano <- ggplot2::ggplot(targets, ggplot2::aes(effect, neglog10_padj_plot)) +
  ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "longdash", color = "#4B5563", linewidth = 0.35) +
  ggplot2::geom_vline(xintercept = 0, color = "#94A3B8", linewidth = 0.3) +
  ggplot2::geom_point(
    ggplot2::aes(color = base_promoter, alpha = hit),
    size = 1.9,
    stroke = 0.25
  ) +
  ggplot2::facet_wrap(ggplot2::vars(target), nrow = 1, scales = "free") +
  ggplot2::scale_color_manual(
    values = c("Fur" = "#D55E00", "MarA" = "#0072B2", "SoxS" = "#009E73", "LexA" = "#CC79A7")
  ) +
  ggplot2::scale_alpha_manual(values = c("FALSE" = 0.38, "TRUE" = 0.95), guide = "none") +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.03, 0.12))) +
  ggplot2::theme_light(base_size = 9) +
  ggplot2::theme(
    legend.position = "bottom",
    panel.grid.minor = ggplot2::element_blank(),
    strip.background = ggplot2::element_rect(fill = "grey95", color = "grey75"),
    strip.text = ggplot2::element_text(face = "bold", size = 9, color = "grey15"),
    plot.margin = ggplot2::margin(8, 34, 8, 8)
  ) +
  ggplot2::labs(
    x = "Estimated effect",
    y = "-log10 adjusted p-value (capped at 35)",
    color = "Reporter"
  )

if (nrow(label_df) > 0) {
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p_volcano <- p_volcano +
      ggrepel::geom_text_repel(
        data = label_df,
        ggplot2::aes(label = label),
        size = 2.25,
        min.segment.length = 0,
        box.padding = 0.24,
        point.padding = 0.12,
        max.overlaps = Inf,
        seed = 11,
        force = 2.5,
        max.iter = 8000,
        segment.size = 0.25,
        show.legend = FALSE
      )
  } else {
    p_volcano <- p_volcano +
      ggplot2::geom_text(
        data = label_df,
        ggplot2::aes(label = label),
        size = 2.15,
        vjust = -0.55,
        check_overlap = TRUE,
        show.legend = FALSE
      )
  }
}
ggplot2::ggsave(
  file.path(out_dir, "dryad_weak_stress_windows_calibrated_alpha_three_target_volcano.png"),
  p_volcano,
  width = 11,
  height = 4.6,
  dpi = 300
)
ggplot2::ggsave(
  file.path(out_dir, "dryad_weak_stress_windows_calibrated_alpha_three_target_volcano.pdf"),
  p_volcano,
  width = 11,
  height = 4.6
)

p_combined <- gridExtra::arrangeGrob(
  p_hist,
  p_volcano,
  ncol = 1,
  heights = c(0.42, 0.58)
)
ggplot2::ggsave(
  file.path(out_dir, "dryad_weak_stress_windows_calibrated_alpha_three_target_inference_combined.png"),
  p_combined,
  width = 11,
  height = 7.6,
  dpi = 300
)
ggplot2::ggsave(
  file.path(out_dir, "dryad_weak_stress_windows_calibrated_alpha_three_target_inference_combined.pdf"),
  p_combined,
  width = 11,
  height = 7.6
)

message("Wrote three-target inference plots to: ", out_dir)
