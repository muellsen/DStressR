#!/usr/bin/env Rscript

# Supplementary diagnostic for the PFIbox analysis: compare raw ordinary
# t-test p-values with empirical-Bayes moderated p-values for the same fitted
# treatment-level DStressR model.

source(file.path("analysis", "_helpers.R"))

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package `ggplot2` is required for the PFIbox moderation diagnostic.", call. = FALSE)
}
if (!requireNamespace("gridExtra", quietly = TRUE)) {
  stop("Package `gridExtra` is required for the PFIbox moderation diagnostic.", call. = FALSE)
}

out_dir <- analysis_output_dir("french_pfibox", "promoter_activity_subset")
tag_moderated <- "reporters2000_perturbations100_T2-18h_treatment_replicate_scalable"
tag_rawse <- paste0(tag_moderated, "_rawse")

moderated_file <- file.path(out_dir, paste0("pfibox_pa_destress_pair_results_", tag_moderated, ".tsv"))
rawse_file <- file.path(out_dir, paste0("pfibox_pa_destress_pair_results_", tag_rawse, ".tsv"))

if (!file.exists(moderated_file)) {
  stop("Missing moderated result table: ", moderated_file, call. = FALSE)
}
if (!file.exists(rawse_file)) {
  stop("Missing unmoderated result table: ", rawse_file, call. = FALSE)
}

moderated <- read_tsv_base(moderated_file)
rawse <- read_tsv_base(rawse_file)

required <- c("reporter", "perturbation", "specific_pvalue", "specific_padj_by_reporter")
missing_moderated <- setdiff(required, names(moderated))
missing_rawse <- setdiff(required, names(rawse))
if (length(missing_moderated)) {
  stop("Moderated table is missing columns: ", paste(missing_moderated, collapse = ", "), call. = FALSE)
}
if (length(missing_rawse)) {
  stop("Unmoderated table is missing columns: ", paste(missing_rawse, collapse = ", "), call. = FALSE)
}

make_panel_data <- function(tab, label) {
  data.frame(
    method = label,
    pvalue = as.numeric(tab$specific_pvalue),
    padj = as.numeric(tab$specific_padj_by_reporter),
    stringsAsFactors = FALSE
  )
}

plot_data <- rbind(
  make_panel_data(rawse, "Ordinary t"),
  make_panel_data(moderated, "Moderated t")
)
plot_data <- plot_data[is.finite(plot_data$pvalue) & plot_data$pvalue >= 0 & plot_data$pvalue <= 1, , drop = FALSE]
plot_data$method <- factor(plot_data$method, levels = c("Ordinary t", "Moderated t"))
plot_data$hit <- is.finite(plot_data$padj) & plot_data$padj < 0.05

summary_tab <- do.call(
  rbind,
  lapply(split(plot_data, plot_data$method), function(d) {
    data.frame(
      method = as.character(d$method[1]),
      n_finite = nrow(d),
      n_p_lt_0_001 = sum(d$pvalue < 0.001, na.rm = TRUE),
      n_p_lt_0_01 = sum(d$pvalue < 0.01, na.rm = TRUE),
      n_hits_reporter_fdr_0_05 = sum(d$hit, na.rm = TRUE),
      median_pvalue = stats::median(d$pvalue, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
)
utils::write.table(
  summary_tab,
  file.path(out_dir, "pfibox_pa_moderation_diagnostic_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

paper_theme <- function(base_size = 10) {
  ggplot2::theme_light(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(linewidth = 0.24, color = "#E5E7EB"),
      axis.title = ggplot2::element_text(color = "#111827"),
      axis.text = ggplot2::element_text(color = "#4B5563"),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", color = "#111827"),
      legend.position = "none",
      plot.title = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_blank()
    )
}

p <- ggplot2::ggplot(plot_data, ggplot2::aes(pvalue, fill = hit)) +
  ggplot2::geom_histogram(
    bins = 50,
    color = "white",
    linewidth = 0.14,
    boundary = 0
  ) +
  ggplot2::facet_wrap(~method, nrow = 1) +
  ggplot2::scale_fill_manual(
    values = c(`FALSE` = "#4E79A7", `TRUE` = "#1D4ED8")
  ) +
  ggplot2::scale_x_continuous(limits = c(0, 1), expand = ggplot2::expansion(mult = c(0, 0.01))) +
  ggplot2::labs(
    x = "Specific-effect p-value",
    y = "Number of reporter-treatment pairs"
  ) +
  paper_theme(base_size = 10)

png_path <- file.path(out_dir, "pfibox_pa_moderation_diagnostic_pvalue_histograms.png")
pdf_path <- file.path(out_dir, "pfibox_pa_moderation_diagnostic_pvalue_histograms.pdf")

ggplot2::ggsave(png_path, p, width = 7.4, height = 3.2, dpi = 320, bg = "white")
ggplot2::ggsave(pdf_path, p, width = 7.4, height = 3.2, device = grDevices::pdf, bg = "white")

message("Wrote PFIbox moderation diagnostic:")
message("  ", png_path)
message("  ", pdf_path)
