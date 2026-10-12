#!/usr/bin/env Rscript

# Paper-style diagnostics for the PFIbox biomass-proxy construction. The plots
# focus only on the response-exponent stage, before perturbation inference.

source(file.path("analysis", "_helpers.R"))
load_destress_package()

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package `ggplot2` is required for PFIbox biomass diagnostics.", call. = FALSE)
}
if (!requireNamespace("ggrepel", quietly = TRUE)) {
  stop("Package `ggrepel` is required for PFIbox biomass diagnostics.", call. = FALSE)
}
if (!requireNamespace("gridExtra", quietly = TRUE)) {
  stop("Package `gridExtra` is required for PFIbox biomass diagnostics.", call. = FALSE)
}

ggplot2 <- asNamespace("ggplot2")
gridExtra <- asNamespace("gridExtra")

processed_dir <- analysis_path("data-raw", "french_pfibox", "processed")
out_dir <- analysis_output_dir("french_pfibox")
auc_file <- file.path(processed_dir, "pfibox_colony_auc_for_destress.tsv")
if (!file.exists(auc_file)) {
  stop("Missing PFIbox AUC table: ", auc_file, call. = FALSE)
}

reporter_limit <- suppressWarnings(as.integer(Sys.getenv("PFIBOX_REPORTER_LIMIT", "40")))
run_tag <- paste0("full_reporters", reporter_limit, "_perturbations6")
alpha_file <- file.path(out_dir, paste0("pfibox_", run_tag, "_alpha_estimated_growth_exponents.tsv"))
if (!file.exists(alpha_file)) {
  stop("Missing alpha estimates for the smoke-test run: ", alpha_file, call. = FALSE)
}

auc <- utils::read.delim(auc_file, check.names = FALSE)
reporters <- sort(unique(auc$reporter))
keep_reporters <- reporters[seq_len(min(reporter_limit, length(reporters)))]
auc <- auc[auc$reporter %in% keep_reporters, , drop = FALSE]
alpha <- utils::read.delim(alpha_file, check.names = FALSE)

control <- auc[auc$perturbation == "Control", , drop = FALSE]
control$log_reporter_auc <- log2(control$fluorescence_auc + 1e-8)
control$log_biomass_auc <- log2(control$biomass_proxy_auc + 1e-8)
control <- control[
  is.finite(control$log_reporter_auc) &
    is.finite(control$log_biomass_auc),
  ,
  drop = FALSE
]

select_reporters <- unique(c(
  utils::head(alpha$reporter[order(alpha$alpha_shrunk)], 3),
  utils::head(alpha$reporter[order(-alpha$alpha_shrunk)], 3),
  utils::head(alpha$reporter[order(abs(alpha$alpha_shrunk))], 3)
))
control_examples <- control[control$reporter %in% select_reporters, , drop = FALSE]
alpha_examples <- alpha[match(select_reporters, alpha$reporter), , drop = FALSE]
control_examples$alpha_label <- paste0(
  control_examples$reporter,
  " (", sprintf("%.2f", alpha_examples$alpha_shrunk[
    match(control_examples$reporter, alpha_examples$reporter)
  ]), ")"
)
control_examples$alpha_label <- factor(
  control_examples$alpha_label,
  levels = unique(control_examples$alpha_label)
)

p_slopes <- ggplot2$ggplot(
  control_examples,
  ggplot2$aes(x = log_biomass_auc, y = log_reporter_auc)
) +
  ggplot2$geom_point(
    ggplot2$aes(color = drug),
    size = 1.35,
    alpha = 0.78
  ) +
  ggplot2$geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = FALSE,
    linewidth = 0.45,
    color = "#111827"
  ) +
  ggplot2$facet_wrap(ggplot2$vars(alpha_label), scales = "free", ncol = 3) +
  ggplot2$scale_color_manual(
    values = c(
      Ampicillin = "#E76F51",
      Azithromycin = "#2A9D8F",
      Cefotaxime = "#4E79A7",
      Chloramphenicol = "#F4A261",
      Ciprofloxacin = "#7C3AED",
      Clarithromycin = "#8A8F3A",
      Clindamycin = "#C4499A",
      Doxycycline = "#6B7280",
      Erythromycin = "#B91C1C",
      Gentamicin = "#1D4ED8",
      Kanamycin = "#0F766E",
      Lincomycin = "#D97706",
      Minocycline = "#9333EA",
      NalidixicAcid = "#374151",
      Neomycin = "#059669",
      PenicillinG = "#DC2626",
      PolymyxinB = "#2563EB",
      Sulfamethoxazole = "#65A30D",
      Trimethoprim = "#EA580C",
      TrimSulfa4 = "#0891B2"
    ),
    guide = "none"
  ) +
  ggplot2$labs(
    x = expression(log[2]~"biomass-proxy AUC"),
    y = expression(log[2]~"reporter fluorescence AUC")
  ) +
  ggplot2$theme_light(base_size = 10) +
  ggplot2$theme(
    panel.grid.minor = ggplot2$element_blank(),
    panel.grid.major = ggplot2$element_line(linewidth = 0.25, color = "#E5E7EB"),
    strip.background = ggplot2$element_rect(fill = "#111827", color = NA),
    strip.text = ggplot2$element_text(color = "white", face = "bold", size = 8.5),
    axis.title = ggplot2$element_text(color = "#111827", size = 9),
    axis.text = ggplot2$element_text(color = "#4B5563", size = 7.5),
    plot.margin = ggplot2$margin(4, 6, 4, 6)
  )

alpha$abs_alpha <- abs(alpha$alpha_shrunk)
alpha$log_growth_sd_plot <- alpha$log_growth_sd
label_alpha <- alpha[order(-alpha$abs_alpha), , drop = FALSE]
label_alpha <- utils::head(label_alpha, 8)

p_alpha <- ggplot2$ggplot(
  alpha,
  ggplot2$aes(x = log_growth_sd_plot, y = alpha_shrunk)
) +
  ggplot2$geom_hline(yintercept = 1, color = "#6B7280", linewidth = 0.35, linetype = "dashed") +
  ggplot2$geom_hline(yintercept = 0, color = "#9CA3AF", linewidth = 0.25) +
  ggplot2$geom_point(
    ggplot2$aes(size = control_n, color = alpha_shrunk_se),
    alpha = 0.84
  ) +
  ggrepel::geom_text_repel(
    data = label_alpha,
    ggplot2$aes(label = reporter),
    size = 3.0,
    fontface = "bold",
    min.segment.length = 0,
    box.padding = 0.22,
    point.padding = 0.12,
    max.overlaps = Inf,
    seed = 1,
    show.legend = FALSE
  ) +
  ggplot2$scale_color_gradient(
    low = "#2A9D8F",
    high = "#D95F02",
    name = expression("SE(" * hat(alpha)[a] * ")")
  ) +
  ggplot2$scale_size_continuous(range = c(1.1, 3.0), name = "Control colonies") +
  ggplot2$labs(
    x = expression("SD of " * log[2]~"biomass-proxy AUC in controls"),
    y = expression(hat(alpha)[a])
  ) +
  ggplot2$theme_light(base_size = 10) +
  ggplot2$theme(
    panel.grid.minor = ggplot2$element_blank(),
    panel.grid.major = ggplot2$element_line(linewidth = 0.25, color = "#E5E7EB"),
    axis.title = ggplot2$element_text(color = "#111827", size = 9),
    axis.text = ggplot2$element_text(color = "#4B5563", size = 8),
    legend.position = "right",
    legend.title = ggplot2$element_text(size = 8),
    legend.text = ggplot2$element_text(size = 7),
    plot.margin = ggplot2$margin(4, 6, 4, 6)
  )

combined <- gridExtra$arrangeGrob(
  p_slopes,
  p_alpha,
  ncol = 1,
  heights = c(1.25, 1)
)

ggplot2$ggsave(
  file.path(out_dir, paste0("pfibox_", run_tag, "_biomass_proxy_diagnostics.png")),
  combined,
  width = 10.5,
  height = 9.2,
  dpi = 280,
  bg = "white"
)
ggplot2$ggsave(
  file.path(out_dir, paste0("pfibox_", run_tag, "_biomass_proxy_diagnostics.pdf")),
  combined,
  width = 10.5,
  height = 9.2,
  bg = "white"
)

message("Wrote PFIbox biomass-proxy diagnostics to: ", out_dir)
