source(file.path("analysis", "_helpers.R"))
load_destress_package()

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package `ggplot2` is required for E. coli EVC-calibrated plots.", call. = FALSE)
}
if (!requireNamespace("MASS", quietly = TRUE)) {
  stop("Package `MASS` is required for Huber EVC calibration.", call. = FALSE)
}

ggplot2 <- asNamespace("ggplot2")

load(analysis_path("data", "binsfeld_reporter_data.rda"))

out_dir <- analysis_output_dir("binsfeld_evc_calibrated")
three_method_dir <- analysis_output_dir("binsfeld_three_method")

safe_neglog10 <- function(x) -log10(pmax(as.numeric(x), .Machine$double.xmin))

wt_auc <- binsfeld_reporter_auc[
  binsfeld_reporter_auc$strain == "WT" &
    binsfeld_reporter_auc$removed == "No",
]

assay <- prepare_assay(
  wt_auc,
  promoter = "promoter",
  compound = "compound",
  control = "Water",
  lux = "lux_auc",
  growth = "od_auc",
  growth_exponent = "estimate",
  batch = "concentration_index",
  replicate = "replicate",
  background_promoter = "EVC",
  background_method = "huber",
  background_by = c("compound", "concentration_index", "replicate")
)

fit <- fit_destress(
  assay,
  technical = c("replicate", "concentration_index"),
  empirical_bayes = TRUE,
  adjustment = "by_promoter",
  interaction = FALSE
)

res <- results(fit)
res$evc_huber_hit_class <- call_hits(
  res,
  fdr = 0.05,
  effect = "specific_effect",
  padj = "specific_padj_by_promoter"
)$hit
res$evc_huber_hit <- res$evc_huber_hit_class != "Not DE"

params <- model_parameters(fit)

utils::write.table(
  params$background_calibration,
  file.path(out_dir, "evc_huber_background_calibration.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
if (!is.null(params$growth_exponents)) {
  utils::write.table(
    params$growth_exponents,
    file.path(out_dir, "evc_huber_growth_exponents.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

pair_results <- res[, c(
  "promoter", "compound", "total_effect", "global_effect", "low_rank_effect",
  "specific_effect", "specific_se", "specific_statistic", "specific_pvalue",
  "specific_padj_by_promoter", "evc_huber_hit", "evc_huber_hit_class"
)]
utils::write.table(
  pair_results,
  file.path(out_dir, "evc_huber_pair_results.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
utils::write.table(
  pair_results[pair_results$evc_huber_hit, , drop = FALSE],
  file.path(out_dir, "evc_huber_significant_pairs.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

comparison_file <- file.path(three_method_dir, "binsfeld_three_method_all_pair_comparison.tsv")
if (file.exists(comparison_file)) {
  comparison <- utils::read.delim(comparison_file, check.names = FALSE)
  evc <- pair_results[, c(
    "promoter", "compound", "specific_effect", "specific_pvalue",
    "specific_padj_by_promoter", "evc_huber_hit", "evc_huber_hit_class"
  )]
  names(evc) <- c(
    "promoter", "compound", "evc_huber_effect", "evc_huber_pvalue",
    "evc_huber_padj_by_promoter", "evc_huber_hit", "evc_huber_hit_class"
  )
  comparison <- merge(comparison, evc, by = c("promoter", "compound"), all.x = TRUE, sort = FALSE)
  comparison$evc_huber_hit[is.na(comparison$evc_huber_hit)] <- FALSE
  comparison$binsfeld_hit <- as.logical(comparison$binsfeld_hit)
  comparison$modeled_hit <- as.logical(comparison$modeled_hit)
  comparison$standard_hit <- as.logical(comparison$standard_hit)

  comparison$evc_huber_class <- apply(
    comparison[, c("binsfeld_hit", "modeled_hit", "evc_huber_hit")],
    1,
    function(x) {
      names <- c("E. coli reference", "DStressR modeled", "DStressR EVC-Huber")[as.logical(x)]
      if (length(names) == 0) {
        "None"
      } else {
        paste(names, collapse = " + ")
      }
    }
  )
  comparison$neglog10_binsfeld <- safe_neglog10(comparison$binsfeld_pvalue)
  comparison$neglog10_binsfeld_padj <- safe_neglog10(comparison$binsfeld_padj)
  comparison$neglog10_modeled <- safe_neglog10(comparison$modeled_pvalue)
  comparison$neglog10_evc_huber <- safe_neglog10(comparison$evc_huber_pvalue)

  utils::write.table(
    comparison,
    file.path(out_dir, "evc_huber_comparison_to_binsfeld_and_default.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  primary_dstressr_only <- comparison[
    !comparison$binsfeld_hit & (comparison$modeled_hit | comparison$evc_huber_hit),
    ,
    drop = FALSE
  ]
  literature_file <- file.path(three_method_dir, "dstressr_only_literature_support.tsv")
  if (file.exists(literature_file)) {
    literature <- utils::read.delim(literature_file, check.names = FALSE)
    literature <- literature[, c("promoter", "compound", "literature_support", "support_note", "source")]
    primary_dstressr_only <- merge(
      primary_dstressr_only,
      literature,
      by = c("promoter", "compound"),
      all.x = TRUE,
      sort = FALSE
    )
    primary_dstressr_only$literature_support[is.na(primary_dstressr_only$literature_support)] <-
      "No direct prior support identified in current literature table"
    primary_dstressr_only$support_note[is.na(primary_dstressr_only$support_note)] <- ""
    primary_dstressr_only$source[is.na(primary_dstressr_only$source)] <- ""
  }
  utils::write.table(
    primary_dstressr_only,
    file.path(out_dir, "primary_dstressr_only_literature_support.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  summary <- data.frame(
    metric = c(
      "Promoter-compound pairs tested",
      "Reference significant pairs",
      "DStressR modeled-response significant pairs",
      "DStressR EVC-Huber significant pairs",
      "Reference and EVC-Huber overlap",
      "Modeled and EVC-Huber overlap",
      "All three overlap",
      "EVC-Huber only vs reference/default"
    ),
    count = c(
      nrow(comparison),
      sum(comparison$binsfeld_hit, na.rm = TRUE),
      sum(comparison$modeled_hit, na.rm = TRUE),
      sum(comparison$evc_huber_hit, na.rm = TRUE),
      sum(comparison$binsfeld_hit & comparison$evc_huber_hit, na.rm = TRUE),
      sum(comparison$modeled_hit & comparison$evc_huber_hit, na.rm = TRUE),
      sum(comparison$binsfeld_hit & comparison$modeled_hit & comparison$evc_huber_hit, na.rm = TRUE),
      sum(!comparison$binsfeld_hit & !comparison$modeled_hit & comparison$evc_huber_hit, na.rm = TRUE)
    )
  )
  utils::write.table(
    summary,
    file.path(out_dir, "evc_huber_hit_summary.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  region_counts <- c(
    binsfeld_only = sum(comparison$binsfeld_hit & !comparison$modeled_hit & !comparison$evc_huber_hit, na.rm = TRUE),
    modeled_only = sum(!comparison$binsfeld_hit & comparison$modeled_hit & !comparison$evc_huber_hit, na.rm = TRUE),
    evc_huber_only = sum(!comparison$binsfeld_hit & !comparison$modeled_hit & comparison$evc_huber_hit, na.rm = TRUE),
    binsfeld_modeled = sum(comparison$binsfeld_hit & comparison$modeled_hit & !comparison$evc_huber_hit, na.rm = TRUE),
    binsfeld_evc_huber = sum(comparison$binsfeld_hit & !comparison$modeled_hit & comparison$evc_huber_hit, na.rm = TRUE),
    modeled_evc_huber = sum(!comparison$binsfeld_hit & comparison$modeled_hit & comparison$evc_huber_hit, na.rm = TRUE),
    all_three = sum(comparison$binsfeld_hit & comparison$modeled_hit & comparison$evc_huber_hit, na.rm = TRUE)
  )
  set_counts <- c(
    binsfeld = sum(comparison$binsfeld_hit, na.rm = TRUE),
    modeled = sum(comparison$modeled_hit, na.rm = TRUE),
    evc_huber = sum(comparison$evc_huber_hit, na.rm = TRUE)
  )
  utils::write.table(
    data.frame(region = names(region_counts), count = as.integer(region_counts)),
    file.path(out_dir, "evc_huber_region_counts.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  circle_points <- function(cx, cy, r, n = 240) {
    theta <- seq(0, 2 * pi, length.out = n)
    data.frame(x = cx + r * cos(theta), y = cy + r * sin(theta))
  }
  venn_df <- rbind(
    cbind(circle_points(-0.55, 0.25, 0.9), method = "E. coli reference"),
    cbind(circle_points(0.55, 0.25, 0.9), method = "DStressR alpha_g"),
    cbind(circle_points(0, -0.48, 0.9), method = "DStressR alpha_g + EVC-Huber")
  )
  venn <- ggplot2$ggplot(venn_df, ggplot2$aes(x, y, fill = method, color = method)) +
    ggplot2$geom_polygon(alpha = 0.22, linewidth = 0.75) +
    ggplot2$annotate("text", x = -0.86, y = 0.55, label = region_counts[["binsfeld_only"]], size = 6, fontface = "bold") +
    ggplot2$annotate("text", x = 0.86, y = 0.55, label = region_counts[["modeled_only"]], size = 6, fontface = "bold") +
    ggplot2$annotate("text", x = 0, y = -1.05, label = region_counts[["evc_huber_only"]], size = 6, fontface = "bold") +
    ggplot2$annotate("text", x = 0, y = 0.62, label = region_counts[["binsfeld_modeled"]], size = 6, fontface = "bold") +
    ggplot2$annotate("text", x = -0.42, y = -0.22, label = region_counts[["binsfeld_evc_huber"]], size = 6, fontface = "bold") +
    ggplot2$annotate("text", x = 0.42, y = -0.22, label = region_counts[["modeled_evc_huber"]], size = 6, fontface = "bold") +
    ggplot2$annotate("text", x = 0, y = 0.1, label = region_counts[["all_three"]], size = 6.3, fontface = "bold") +
    ggplot2$annotate("text", x = -0.75, y = 1.35, label = paste0("E. coli reference\n", set_counts[["binsfeld"]], " hits"), size = 3.7, fontface = "bold") +
    ggplot2$annotate("text", x = 0.75, y = 1.35, label = paste0("DStressR alpha_g\n", set_counts[["modeled"]], " hits"), size = 3.7, fontface = "bold") +
    ggplot2$annotate("text", x = 0, y = -1.65, label = paste0("DStressR alpha_g + EVC-Huber\n", set_counts[["evc_huber"]], " hits"), size = 3.7, fontface = "bold") +
    ggplot2$scale_fill_manual(values = c(
      "E. coli reference" = "#2563eb",
      "DStressR alpha_g" = "#dc2626",
      "DStressR alpha_g + EVC-Huber" = "#059669"
    )) +
    ggplot2$scale_color_manual(values = c(
      "E. coli reference" = "#1d4ed8",
      "DStressR alpha_g" = "#b91c1c",
      "DStressR alpha_g + EVC-Huber" = "#047857"
    )) +
    ggplot2$coord_equal(xlim = c(-1.75, 1.75), ylim = c(-1.85, 1.65), expand = FALSE) +
    ggplot2$theme_void(base_size = 10) +
    ggplot2$theme(legend.position = "none") +
    ggplot2$labs(
      title = "Hit overlap across E. coli reference and DStressR response-construction analyses"
    )

  volcano_data <- rbind(
    data.frame(
      method = "E. coli reference",
      promoter = comparison$promoter,
      compound = comparison$compound,
      effect = comparison$mean_z,
      neglog10_pvalue = comparison$neglog10_binsfeld,
      hit = comparison$binsfeld_hit,
      stringsAsFactors = FALSE
    ),
    data.frame(
      method = "DStressR modeled",
      promoter = comparison$promoter,
      compound = comparison$compound,
      effect = comparison$modeled_effect,
      neglog10_pvalue = comparison$neglog10_modeled,
      hit = comparison$modeled_hit,
      stringsAsFactors = FALSE
    ),
    data.frame(
      method = "DStressR EVC-Huber",
      promoter = comparison$promoter,
      compound = comparison$compound,
      effect = comparison$evc_huber_effect,
      neglog10_pvalue = comparison$neglog10_evc_huber,
      hit = comparison$evc_huber_hit,
      stringsAsFactors = FALSE
    )
  )
  volcano_data$method <- factor(
    volcano_data$method,
    levels = c("E. coli reference", "DStressR modeled", "DStressR EVC-Huber")
  )
  volcano_data$status <- ifelse(volcano_data$hit, "Hit", "Not called")
  volcano_data$label <- paste(volcano_data$promoter, volcano_data$compound, sep = "-")
  volcano_data$label_hjust <- ifelse(volcano_data$effect > 0, 1.05, -0.05)

  top_volcano_labels <- do.call(rbind, lapply(split(volcano_data, volcano_data$method), function(d) {
    d <- d[d$hit & is.finite(d$effect) & is.finite(d$neglog10_pvalue), , drop = FALSE]
    d <- d[order(-d$neglog10_pvalue, -abs(d$effect)), , drop = FALSE]
    utils::head(d, 6)
  }))

  reference_volcano <- comparison[, c(
    "promoter", "compound", "mean_z", "binsfeld_padj",
    "neglog10_binsfeld_padj", "binsfeld_hit"
  )]
  reference_volcano$status <- ifelse(reference_volcano$binsfeld_hit, reference_volcano$promoter, "Not called")
  reference_volcano$label <- paste(reference_volcano$promoter, reference_volcano$compound, sep = "-")
  reference_volcano$label_hjust <- ifelse(reference_volcano$mean_z > 0, 1.05, -0.05)
  reference_labels <- reference_volcano[
    reference_volcano$binsfeld_hit &
      is.finite(reference_volcano$mean_z) &
      is.finite(reference_volcano$neglog10_binsfeld_padj),
    ,
    drop = FALSE
  ]
  reference_labels <- reference_labels[
    order(-reference_labels$neglog10_binsfeld_padj, -abs(reference_labels$mean_z)),
    ,
    drop = FALSE
  ]
  reference_labels <- utils::head(reference_labels, 12)

  promoter_colors <- c(
    acrABp = "#2563eb",
    marRABp = "#dc2626",
    micFp = "#059669",
    ompFp = "#7c3aed",
    robp = "#d97706",
    soxSp = "#0891b2",
    tolCp = "#be123c",
    "Not called" = "#d1d5db"
  )

  reference_plot <- ggplot2$ggplot(
    reference_volcano,
    ggplot2$aes(mean_z, neglog10_binsfeld_padj, color = status)
  ) +
    ggplot2$geom_hline(yintercept = -log10(0.05), color = "#9ca3af", linewidth = 0.35, linetype = "dashed") +
    ggplot2$geom_vline(xintercept = c(-1, 1), color = "#9ca3af", linewidth = 0.35, linetype = "dotted") +
    ggplot2$geom_vline(xintercept = 0, color = "#9ca3af", linewidth = 0.3) +
    ggplot2$geom_point(alpha = 0.82, size = 1.65) +
    ggplot2$geom_text(
      data = reference_labels,
      ggplot2$aes(label = label, hjust = label_hjust),
      color = "#111827",
      size = 2.45,
      vjust = -0.55,
      check_overlap = TRUE,
      show.legend = FALSE
    ) +
    ggplot2$scale_color_manual(values = promoter_colors, breaks = names(promoter_colors)) +
    ggplot2$scale_x_continuous(expand = ggplot2$expansion(mult = c(0.08, 0.1))) +
    ggplot2$scale_y_continuous(expand = ggplot2$expansion(mult = c(0.04, 0.16))) +
    ggplot2$theme_light(base_size = 10) +
    ggplot2$theme(
      panel.grid.minor = ggplot2$element_blank(),
      plot.title = ggplot2$element_text(face = "bold"),
      legend.position = "bottom"
    ) +
    ggplot2$labs(
      title = "E. coli reference volcano plot",
      subtitle = "Reference rule: promoter-wise BH-adjusted Wilcoxon p < 0.05 and absolute mean Z-score > 1",
      x = "Mean Z-score",
      y = "-log10 promoter-wise BH adjusted p-value",
      color = "Promoter"
    )

  volcano_plot <- ggplot2$ggplot(
    volcano_data,
    ggplot2$aes(effect, neglog10_pvalue, color = status)
  ) +
    ggplot2$geom_hline(yintercept = -log10(0.05), color = "#9ca3af", linewidth = 0.3, linetype = "dashed") +
    ggplot2$geom_vline(xintercept = 0, color = "#9ca3af", linewidth = 0.3) +
    ggplot2$geom_point(alpha = 0.78, size = 1.45) +
    ggplot2$geom_text(
      data = top_volcano_labels,
      ggplot2$aes(label = label, hjust = label_hjust),
      color = "#111827",
      size = 2.3,
      vjust = -0.55,
      check_overlap = TRUE,
      show.legend = FALSE
    ) +
    ggplot2$facet_wrap(ggplot2$vars(method), scales = "free_x", ncol = 3) +
    ggplot2$scale_color_manual(values = c("Not called" = "#d1d5db", "Hit" = "#111827")) +
    ggplot2$scale_x_continuous(expand = ggplot2$expansion(mult = c(0.08, 0.12))) +
    ggplot2$scale_y_continuous(expand = ggplot2$expansion(mult = c(0.04, 0.16))) +
    ggplot2$theme_light(base_size = 10) +
    ggplot2$theme(
      panel.grid.minor = ggplot2$element_blank(),
      plot.title = ggplot2$element_text(face = "bold"),
      legend.position = "bottom",
      plot.margin = ggplot2$margin(8, 12, 8, 8)
    ) +
    ggplot2$labs(
      title = "E. coli reference, modeled-response, and EVC-Huber volcano plots",
      x = "Effect score",
      y = "-log10 raw p-value",
      color = "Call"
    )

  ggplot2$ggsave(file.path(out_dir, "evc_huber_hit_overlap_venn.png"), venn, width = 7.2, height = 6.8, dpi = 220)
  ggplot2$ggsave(file.path(out_dir, "evc_huber_hit_overlap_venn.pdf"), venn, width = 7.2, height = 6.8)
  ggplot2$ggsave(file.path(out_dir, "ecoli_reference_volcano_plot.png"), reference_plot, width = 8.2, height = 5.6, dpi = 220)
  ggplot2$ggsave(file.path(out_dir, "ecoli_reference_volcano_plot.pdf"), reference_plot, width = 8.2, height = 5.6)
  ggplot2$ggsave(file.path(out_dir, "evc_huber_volcano_plots.png"), volcano_plot, width = 12, height = 4.8, dpi = 220)
  ggplot2$ggsave(file.path(out_dir, "evc_huber_volcano_plots.pdf"), volcano_plot, width = 12, height = 4.8)

  print(summary)
}

message("Wrote E. coli EVC-Huber analysis to: ", out_dir)
