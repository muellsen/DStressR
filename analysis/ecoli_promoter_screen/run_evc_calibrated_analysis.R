source(file.path("analysis", "_helpers.R"))
load_destress_package()

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package `ggplot2` is required for E. coli EVC-calibrated plots.", call. = FALSE)
}
if (!requireNamespace("MASS", quietly = TRUE)) {
  stop("Package `MASS` is required for Huber EVC calibration.", call. = FALSE)
}
if (!requireNamespace("gridExtra", quietly = TRUE)) {
  stop("Package `gridExtra` is required for E. coli multi-panel figures.", call. = FALSE)
}
if (!requireNamespace("ggrepel", quietly = TRUE)) {
  stop("Package `ggrepel` is required for E. coli volcano label placement.", call. = FALSE)
}

ggplot2 <- asNamespace("ggplot2")
gridExtra <- asNamespace("gridExtra")
ggrepel <- asNamespace("ggrepel")

panel_label <- function(label, size = 16) {
  grid::textGrob(
    label,
    x = grid::unit(0, "npc"),
    y = grid::unit(0.08, "npc"),
    hjust = 0,
    vjust = 0,
    gp = grid::gpar(fontsize = size, fontface = "bold", col = "#111827")
  )
}

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
  batch = "dose_level",
  replicate = "replicate",
  growth_covariates = "replicate",
  numeric_covariates = "dose_level",
  background_promoter = "EVC",
  background_method = "huber",
  background_by = c("compound", "dose_level", "replicate")
)

fit <- fit_destress(
  assay,
  technical = c("replicate", "dose_level"),
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
      names <- c("Binsfeld reference", "DStressR without EV", "DStressR with EV")[as.logical(x)]
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
      "DStressR without EV significant pairs",
      "DStressR with EV significant pairs",
      "Reference and DStressR with EV overlap",
      "DStressR workflows overlap",
      "All three overlap",
      "DStressR with EV only vs reference/without EV"
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
    cbind(circle_points(-0.55, 0.25, 0.9), method = "Binsfeld reference"),
    cbind(circle_points(0.55, 0.25, 0.9), method = "DStressR without EV"),
    cbind(circle_points(0, -0.48, 0.9), method = "DStressR with EV")
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
    ggplot2$annotate("text", x = -0.75, y = 1.35, label = paste0("Binsfeld reference\n", set_counts[["binsfeld"]], " hits"), size = 3.7, fontface = "bold") +
    ggplot2$annotate("text", x = 0.75, y = 1.35, label = paste0("DStressR without EV\n", set_counts[["modeled"]], " hits"), size = 3.7, fontface = "bold") +
    ggplot2$annotate("text", x = 0, y = -1.65, label = paste0("DStressR with EV\n", set_counts[["evc_huber"]], " hits"), size = 3.7, fontface = "bold") +
    ggplot2$scale_fill_manual(values = c(
      "Binsfeld reference" = "#6f7f8f",
      "DStressR without EV" = "#d8b56d",
      "DStressR with EV" = "#9b6a55"
    )) +
    ggplot2$scale_color_manual(values = c(
      "Binsfeld reference" = "#3f5264",
      "DStressR without EV" = "#9f7625",
      "DStressR with EV" = "#633b2d"
    )) +
    ggplot2$coord_equal(xlim = c(-1.75, 1.75), ylim = c(-1.85, 1.65), expand = FALSE) +
    ggplot2$theme_void(base_size = 10) +
    ggplot2$theme(legend.position = "none")
  venn_compact <- ggplot2$ggplot(venn_df, ggplot2$aes(x, y, fill = method, color = method)) +
    ggplot2$geom_polygon(alpha = 0.22, linewidth = 0.72) +
    ggplot2$annotate("text", x = -0.86, y = 0.55, label = region_counts[["binsfeld_only"]], size = 5.2, fontface = "bold") +
    ggplot2$annotate("text", x = 0.86, y = 0.55, label = region_counts[["modeled_only"]], size = 5.2, fontface = "bold") +
    ggplot2$annotate("text", x = 0, y = -1.05, label = region_counts[["evc_huber_only"]], size = 5.2, fontface = "bold") +
    ggplot2$annotate("text", x = 0, y = 0.62, label = region_counts[["binsfeld_modeled"]], size = 5.2, fontface = "bold") +
    ggplot2$annotate("text", x = -0.42, y = -0.22, label = region_counts[["binsfeld_evc_huber"]], size = 5.2, fontface = "bold") +
    ggplot2$annotate("text", x = 0.42, y = -0.22, label = region_counts[["modeled_evc_huber"]], size = 5.2, fontface = "bold") +
    ggplot2$annotate("text", x = 0, y = 0.1, label = region_counts[["all_three"]], size = 5.5, fontface = "bold") +
    ggplot2$annotate("text", x = -1.05, y = 1.35, label = paste0("Binsfeld reference\n", set_counts[["binsfeld"]], " hits"), size = 2.8, fontface = "bold") +
    ggplot2$annotate("text", x = 1.05, y = 1.35, label = paste0("DStressR without\nEV control\n", set_counts[["modeled"]], " hits"), size = 2.8, fontface = "bold") +
    ggplot2$annotate("text", x = 0, y = -1.66, label = paste0("DStressR with EV\n", set_counts[["evc_huber"]], " hits"), size = 2.8, fontface = "bold") +
    ggplot2$scale_fill_manual(values = c(
      "Binsfeld reference" = "#6f7f8f",
      "DStressR without EV" = "#d8b56d",
      "DStressR with EV" = "#9b6a55"
    )) +
    ggplot2$scale_color_manual(values = c(
      "Binsfeld reference" = "#3f5264",
      "DStressR without EV" = "#9f7625",
      "DStressR with EV" = "#633b2d"
    )) +
    ggplot2$coord_equal(xlim = c(-1.88, 1.88), ylim = c(-1.85, 1.68), expand = FALSE) +
    ggplot2$theme_void(base_size = 10) +
    ggplot2$theme(
      legend.position = "none",
      plot.margin = ggplot2$margin(2, 2, 2, 2)
    )

  volcano_data <- rbind(
    data.frame(
      method = "Binsfeld reference",
      promoter = comparison$promoter,
      compound = comparison$compound,
      effect = comparison$mean_z,
      neglog10_pvalue = comparison$neglog10_binsfeld,
      hit = comparison$binsfeld_hit,
      stringsAsFactors = FALSE
    ),
    data.frame(
      method = "DStressR without EV",
      promoter = comparison$promoter,
      compound = comparison$compound,
      effect = comparison$modeled_effect,
      neglog10_pvalue = comparison$neglog10_modeled,
      hit = comparison$modeled_hit,
      stringsAsFactors = FALSE
    ),
    data.frame(
      method = "DStressR with EV",
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
    levels = c("Binsfeld reference", "DStressR without EV", "DStressR with EV")
  )
  volcano_data$in_all_three <- comparison$binsfeld_hit & comparison$modeled_hit & comparison$evc_huber_hit
  volcano_data$method_hit_count <- ave(
    volcano_data$hit,
    volcano_data$method,
    FUN = function(x) sum(x, na.rm = TRUE)
  )
  volcano_data$status <- ifelse(volcano_data$hit, "Hit", "Not called")
  volcano_data$label <- paste(volcano_data$promoter, volcano_data$compound, sep = "-")
  volcano_data$label_hjust <- ifelse(volcano_data$effect > 0, 1.05, -0.05)
  volcano_data$rank_score <- volcano_data$neglog10_pvalue + 0.15 * abs(volcano_data$effect)
  volcano_y_max <- ceiling(max(
    volcano_data$neglog10_pvalue[is.finite(volcano_data$neglog10_pvalue)],
    na.rm = TRUE
  ) / 5) * 5
  volcano_x_bounds <- do.call(rbind, lapply(split(volcano_data, volcano_data$method), function(d) {
    xmax <- max(abs(d$effect[is.finite(d$effect)]), na.rm = TRUE)
    xmax <- ceiling(xmax * 1.12 * 10) / 10
    data.frame(
      method = unique(d$method),
      effect = c(-xmax, xmax),
      neglog10_pvalue = c(0, volcano_y_max),
      stringsAsFactors = FALSE
    )
  }))
  volcano_x_bounds$method <- factor(volcano_x_bounds$method, levels = levels(volcano_data$method))
  top_volcano_labels <- do.call(rbind, lapply(split(volcano_data, volcano_data$method), function(d) {
    d <- d[d$hit & is.finite(d$effect) & is.finite(d$neglog10_pvalue), , drop = FALSE]
    d <- d[order(-d$neglog10_pvalue, -abs(d$effect)), , drop = FALSE]
    utils::head(d, 6)
  }))

  pvalue_long <- rbind(
    data.frame(
      method = "Binsfeld reference",
      promoter = comparison$promoter,
      pvalue = comparison$binsfeld_pvalue,
      hit = comparison$binsfeld_hit,
      stringsAsFactors = FALSE
    ),
    data.frame(
      method = "DStressR without EV",
      promoter = comparison$promoter,
      pvalue = comparison$modeled_pvalue,
      hit = comparison$modeled_hit,
      stringsAsFactors = FALSE
    ),
    data.frame(
      method = "DStressR with EV",
      promoter = comparison$promoter,
      pvalue = comparison$evc_huber_pvalue,
      hit = comparison$evc_huber_hit,
      stringsAsFactors = FALSE
    )
  )
  pvalue_long$method <- factor(
    pvalue_long$method,
    levels = c("Binsfeld reference", "DStressR without EV", "DStressR with EV")
  )
  pvalue_long$status <- ifelse(pvalue_long$hit, "Called hit", "Not called")
  pvalue_long <- pvalue_long[is.finite(pvalue_long$pvalue), , drop = FALSE]
  promoter_levels <- c("All promoters", "acrABp", "marRABp", "micFp", "ompFp", "robp", "soxSp", "tolCp")
  pvalue_panel <- rbind(
    transform(pvalue_long, panel = "All promoters"),
    transform(pvalue_long, panel = as.character(promoter))
  )
  pvalue_panel$panel <- factor(pvalue_panel$panel, levels = promoter_levels)
  pvalue_hist <- ggplot2$ggplot(pvalue_panel, ggplot2$aes(pvalue)) +
    ggplot2$geom_histogram(
      bins = 30,
      fill = "#d1d5db",
      color = "white",
      linewidth = 0.2
    ) +
    ggplot2$facet_grid(ggplot2$vars(panel), ggplot2$vars(method), scales = "free_y") +
    ggplot2$scale_x_continuous(breaks = seq(0, 1, by = 0.25)) +
    ggplot2$coord_cartesian(xlim = c(0, 1)) +
    ggplot2$theme_light(base_size = 10) +
	  ggplot2$theme(
	    panel.grid.minor = ggplot2$element_blank(),
	    strip.text = ggplot2$element_text(face = "bold", hjust = 0.5)
	  ) +
	  ggplot2$labs(
	    x = "Raw p-value",
	    y = "Promoter-compound pairs"
	  )
  pvalue_promoter_panel <- pvalue_panel[pvalue_panel$panel != "All promoters", , drop = FALSE]
  pvalue_hist_by_promoter <- ggplot2$ggplot(pvalue_promoter_panel, ggplot2$aes(pvalue)) +
    ggplot2$geom_histogram(
      bins = 30,
      fill = "#d1d5db",
      color = "white",
      linewidth = 0.2
    ) +
    ggplot2$facet_grid(ggplot2$vars(panel), ggplot2$vars(method), scales = "free_y") +
    ggplot2$scale_x_continuous(breaks = seq(0, 1, by = 0.25)) +
    ggplot2$coord_cartesian(xlim = c(0, 1)) +
    ggplot2$theme_light(base_size = 10) +
    ggplot2$theme(
      panel.grid.minor = ggplot2$element_blank(),
      strip.text = ggplot2$element_text(face = "bold", hjust = 0.5)
    ) +
    ggplot2$labs(
      x = "Raw p-value",
      y = "Promoter-compound pairs"
    )

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
    acrABp = "#0072B2",
    marRABp = "#D55E00",
    micFp = "#CC79A7",
    ompFp = "#E69F00",
    robp = "#009E73",
    soxSp = "#56B4E9",
    tolCp = "#6A3D9A",
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
	    legend.position = "bottom"
	  ) +
	  ggplot2$labs(
	    x = "Mean Z-score",
	    y = "-log10 promoter-wise BH adjusted p-value",
	    color = "Promoter"
    )

  volcano_plot <- ggplot2$ggplot(
    volcano_data,
    ggplot2$aes(effect, neglog10_pvalue, color = status)
  ) +
    ggplot2$geom_blank(
      data = volcano_x_bounds,
      ggplot2$aes(effect, neglog10_pvalue),
      inherit.aes = FALSE
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
    ggplot2$scale_x_continuous(expand = ggplot2$expansion(mult = c(0.02, 0.02))) +
    ggplot2$scale_y_continuous(limits = c(0, volcano_y_max), expand = ggplot2$expansion(mult = c(0.02, 0.08))) +
    ggplot2$theme_light(base_size = 10) +
	  ggplot2$theme(
	    panel.grid.minor = ggplot2$element_blank(),
	    strip.text = ggplot2$element_text(face = "bold", hjust = 0.5),
	    legend.position = "bottom",
	    plot.margin = ggplot2$margin(8, 12, 8, 8)
	  ) +
	  ggplot2$labs(
	    x = "Method-specific effect estimate",
	    y = "-log10 raw p-value",
	    color = "Call"
    )

  intersection_labels <- volcano_data[
    volcano_data$in_all_three & is.finite(volcano_data$effect) & is.finite(volcano_data$neglog10_pvalue),
    ,
    drop = FALSE
  ]
  intersection_labels <- intersection_labels[
    order(intersection_labels$method, -intersection_labels$rank_score),
    ,
    drop = FALSE
  ]
  intersection_labels <- do.call(rbind, lapply(split(intersection_labels, intersection_labels$method), function(d) {
    utils::head(d, 12)
  }))
  intersection_plot <- ggplot2$ggplot(
    volcano_data,
    ggplot2$aes(effect, neglog10_pvalue)
  ) +
    ggplot2$geom_blank(
      data = volcano_x_bounds,
      ggplot2$aes(effect, neglog10_pvalue),
      inherit.aes = FALSE
    ) +
    ggplot2$geom_hline(yintercept = -log10(0.05), color = "#9ca3af", linewidth = 0.3, linetype = "dashed") +
    ggplot2$geom_vline(xintercept = 0, color = "#9ca3af", linewidth = 0.3) +
    ggplot2$geom_point(color = "#d1d5db", alpha = 0.55, size = 1.2) +
    ggplot2$geom_point(
      data = volcano_data[volcano_data$in_all_three, , drop = FALSE],
      ggplot2$aes(color = promoter),
      alpha = 0.9,
      size = 1.9
    ) +
    ggplot2$geom_text(
      data = intersection_labels,
      ggplot2$aes(label = label, color = promoter, hjust = label_hjust),
      size = 2.15,
      vjust = -0.55,
      check_overlap = TRUE,
      show.legend = FALSE
    ) +
    ggplot2$facet_wrap(ggplot2$vars(method), scales = "free_x", ncol = 3) +
    ggplot2$scale_color_manual(values = promoter_colors[names(promoter_colors) != "Not called"]) +
    ggplot2$scale_x_continuous(expand = ggplot2$expansion(mult = c(0.02, 0.02))) +
    ggplot2$scale_y_continuous(limits = c(0, volcano_y_max), expand = ggplot2$expansion(mult = c(0.02, 0.08))) +
    ggplot2$theme_light(base_size = 10) +
	  ggplot2$theme(
	    panel.grid.minor = ggplot2$element_blank(),
	    strip.text = ggplot2$element_text(face = "bold", hjust = 0.5),
	    legend.position = "bottom",
	    plot.margin = ggplot2$margin(8, 12, 8, 8)
	  ) +
	  ggplot2$labs(
	    x = "Method-specific effect estimate",
	    y = "-log10 raw p-value",
	    color = "Promoter"
    )

  non_intersection <- volcano_data[
    volcano_data$hit & !volcano_data$in_all_three &
      is.finite(volcano_data$effect) & is.finite(volcano_data$neglog10_pvalue),
    ,
    drop = FALSE
  ]
  top_non_intersection <- do.call(rbind, lapply(split(non_intersection, non_intersection$method), function(d) {
    d <- d[order(-d$rank_score), , drop = FALSE]
    utils::head(d, 18)
  }))
  non_intersection$panel_label <- paste0(
    as.character(non_intersection$method),
    "\n",
    ave(non_intersection$hit, non_intersection$method, FUN = function(x) length(x)),
    " non-intersection hits"
  )
  panel_lookup <- stats::setNames(
    as.character(non_intersection$panel_label),
    as.character(non_intersection$method)
  )
  top_non_intersection$panel_label <- unname(panel_lookup[as.character(top_non_intersection$method)])
  top_non_intersection$panel_label <- factor(
    top_non_intersection$panel_label,
    levels = unique(non_intersection$panel_label)
  )
  non_intersection$panel_label <- factor(
    non_intersection$panel_label,
    levels = unique(non_intersection$panel_label)
  )
  non_intersection_x_bounds <- merge(
    volcano_x_bounds,
    unique(non_intersection[, c("method", "panel_label"), drop = FALSE]),
    by = "method",
    all.y = TRUE,
    sort = FALSE
  )
  volcano_panel_data <- merge(
    volcano_data,
    unique(non_intersection[, c("method", "panel_label"), drop = FALSE]),
    by = "method",
    all.y = TRUE,
    sort = FALSE
  )
  non_intersection_plot <- ggplot2$ggplot(
    volcano_panel_data,
    ggplot2$aes(effect, neglog10_pvalue)
  ) +
    ggplot2$geom_blank(
      data = non_intersection_x_bounds,
      ggplot2$aes(effect, neglog10_pvalue),
      inherit.aes = FALSE
    ) +
    ggplot2$geom_hline(yintercept = -log10(0.05), color = "#9ca3af", linewidth = 0.3, linetype = "dashed") +
    ggplot2$geom_vline(xintercept = 0, color = "#9ca3af", linewidth = 0.3) +
    ggplot2$geom_point(color = "#e5e7eb", alpha = 0.45, size = 1.1) +
    ggplot2$geom_point(
      data = non_intersection,
      ggplot2$aes(color = promoter),
      alpha = 0.88,
      size = 1.75
    ) +
    ggplot2$geom_text(
      data = top_non_intersection,
      ggplot2$aes(label = label, color = promoter, hjust = label_hjust),
      size = 2.05,
      vjust = -0.55,
      check_overlap = TRUE,
      show.legend = FALSE
    ) +
    ggplot2$facet_wrap(ggplot2$vars(panel_label), scales = "free_x", ncol = 3) +
    ggplot2$scale_color_manual(values = promoter_colors[names(promoter_colors) != "Not called"]) +
    ggplot2$scale_x_continuous(expand = ggplot2$expansion(mult = c(0.02, 0.02))) +
    ggplot2$scale_y_continuous(limits = c(0, volcano_y_max), expand = ggplot2$expansion(mult = c(0.02, 0.08))) +
    ggplot2$theme_light(base_size = 10) +
	  ggplot2$theme(
	    panel.grid.minor = ggplot2$element_blank(),
	    strip.text = ggplot2$element_text(face = "bold", hjust = 0.5),
	    legend.position = "bottom",
	    plot.margin = ggplot2$margin(8, 12, 8, 8)
	  ) +
	  ggplot2$labs(
	    x = "Method-specific effect estimate",
	    y = "-log10 raw p-value",
	    color = "Promoter"
    )

  pvalue_hist_all <- ggplot2$ggplot(
    pvalue_long,
    ggplot2$aes(pvalue)
  ) +
    ggplot2$geom_histogram(
      bins = 30,
      fill = "#d1d5db",
      color = "white",
      linewidth = 0.2
    ) +
    ggplot2$facet_wrap(ggplot2$vars(method), nrow = 1, scales = "free_y") +
    ggplot2$scale_x_continuous(breaks = seq(0, 1, by = 0.5)) +
    ggplot2$coord_cartesian(xlim = c(0, 1)) +
    ggplot2$theme_light(base_size = 11) +
    ggplot2$theme(
      panel.grid.minor = ggplot2$element_blank(),
      strip.text = ggplot2$element_text(face = "bold", hjust = 0.5, size = 10),
      plot.margin = ggplot2$margin(4, 4, 4, 4)
    ) +
    ggplot2$labs(
      x = "Raw p-value",
      y = "Pairs"
    )

  overlap_volcano_labels <- volcano_data[
    volcano_data$in_all_three & is.finite(volcano_data$effect) & is.finite(volcano_data$neglog10_pvalue),
    ,
    drop = FALSE
  ]
  overlap_volcano_labels$short_label <- paste(overlap_volcano_labels$promoter, overlap_volcano_labels$compound, sep = "-")
  overlap_volcano_labels$label_side <- ifelse(overlap_volcano_labels$effect >= 0, 1, -1)
  overlap_volcano_labels$nudge_x <- ifelse(
    overlap_volcano_labels$method == "Binsfeld reference",
    0.32 * overlap_volcano_labels$label_side,
    0.10 * overlap_volcano_labels$label_side
  )
  overlap_volcano_labels$nudge_y <- ifelse(
    overlap_volcano_labels$method == "Binsfeld reference",
    pmax(5.5, volcano_y_max * 0.20),
    pmax(0.45, volcano_y_max * 0.035)
  )
  binsfeld_overlap_labels <- overlap_volcano_labels[
    overlap_volcano_labels$method == "Binsfeld reference",
    ,
    drop = FALSE
  ]
  dstressr_overlap_labels <- overlap_volcano_labels[
    overlap_volcano_labels$method != "Binsfeld reference",
    ,
    drop = FALSE
  ]

  volcano_overlap_plot <- ggplot2$ggplot(
    volcano_data,
    ggplot2$aes(effect, neglog10_pvalue)
  ) +
    ggplot2$geom_blank(
      data = volcano_x_bounds,
      ggplot2$aes(effect, neglog10_pvalue),
      inherit.aes = FALSE
    ) +
    ggplot2$geom_hline(yintercept = -log10(0.05), color = "#9ca3af", linewidth = 0.28, linetype = "dashed") +
    ggplot2$geom_vline(xintercept = 0, color = "#9ca3af", linewidth = 0.28) +
    ggplot2$geom_point(color = "#d1d5db", alpha = 0.42, size = 0.9) +
    ggplot2$geom_point(
      data = volcano_data[volcano_data$in_all_three, , drop = FALSE],
      ggplot2$aes(color = promoter),
      alpha = 0.95,
      size = 1.55
    ) +
    ggrepel$geom_text_repel(
      data = binsfeld_overlap_labels,
      ggplot2$aes(label = short_label, color = promoter),
      size = 1.55,
      min.segment.length = 0,
      segment.color = "#6b7280",
      segment.size = 0.16,
      segment.alpha = 0.78,
      box.padding = 0.20,
      point.padding = 0.10,
      force = 3.2,
      force_pull = 0.025,
      max.overlaps = Inf,
      max.time = 4,
      max.iter = 12000,
      nudge_x = binsfeld_overlap_labels$nudge_x,
      nudge_y = binsfeld_overlap_labels$nudge_y,
      seed = 3260,
      show.legend = FALSE
    ) +
    ggrepel$geom_text_repel(
      data = dstressr_overlap_labels,
      ggplot2$aes(label = short_label, color = promoter),
      size = 1.55,
      min.segment.length = 0,
      segment.color = "#6b7280",
      segment.size = 0.16,
      segment.alpha = 0.72,
      box.padding = 0.16,
      point.padding = 0.10,
      force = 1.7,
      force_pull = 0.20,
      max.overlaps = Inf,
      max.time = 3,
      max.iter = 10000,
      nudge_x = dstressr_overlap_labels$nudge_x,
      nudge_y = dstressr_overlap_labels$nudge_y,
      seed = 3261,
      show.legend = FALSE
    ) +
    ggplot2$facet_wrap(ggplot2$vars(method), scales = "free_x", nrow = 1) +
    ggplot2$scale_color_manual(values = promoter_colors[names(promoter_colors) != "Not called"]) +
    ggplot2$scale_x_continuous(expand = ggplot2$expansion(mult = c(0.04, 0.1))) +
    ggplot2$scale_y_continuous(limits = c(0, volcano_y_max), expand = ggplot2$expansion(mult = c(0.02, 0.12))) +
    ggplot2$theme_light(base_size = 11) +
    ggplot2$theme(
      panel.grid.minor = ggplot2$element_blank(),
      strip.text = ggplot2$element_text(face = "bold", hjust = 0.5, size = 10),
      legend.position = "none",
      plot.margin = ggplot2$margin(4, 4, 4, 4)
    ) +
    ggplot2$labs(
      x = "Method-specific effect estimate",
      y = "-log10 raw p-value"
    )

  pvalue_scatter_long <- rbind(
    data.frame(
      comparison = "DStressR without EV",
      promoter = comparison$promoter,
      compound = comparison$compound,
      binsfeld_pvalue = comparison$binsfeld_pvalue,
      destress_pvalue = comparison$modeled_pvalue,
      binsfeld_neglog10_pvalue = comparison$neglog10_binsfeld,
      destress_neglog10_pvalue = comparison$neglog10_modeled,
      binsfeld_hit = comparison$binsfeld_hit,
      modeled_hit = comparison$modeled_hit,
      evc_huber_hit = comparison$evc_huber_hit,
      destress_hit = comparison$modeled_hit,
      stringsAsFactors = FALSE
    ),
    data.frame(
      comparison = "DStressR with EV",
      promoter = comparison$promoter,
      compound = comparison$compound,
      binsfeld_pvalue = comparison$binsfeld_pvalue,
      destress_pvalue = comparison$evc_huber_pvalue,
      binsfeld_neglog10_pvalue = comparison$neglog10_binsfeld,
      destress_neglog10_pvalue = comparison$neglog10_evc_huber,
      binsfeld_hit = comparison$binsfeld_hit,
      modeled_hit = comparison$modeled_hit,
      evc_huber_hit = comparison$evc_huber_hit,
      destress_hit = comparison$evc_huber_hit,
      stringsAsFactors = FALSE
    )
  )
  pvalue_scatter_long$comparison <- factor(
    pvalue_scatter_long$comparison,
    levels = c("DStressR without EV", "DStressR with EV")
  )
  pvalue_scatter_long$disagreement <- abs(
    pvalue_scatter_long$destress_neglog10_pvalue - pvalue_scatter_long$binsfeld_neglog10_pvalue
  )
  pvalue_scatter_long$binsfeld_only_all <- pvalue_scatter_long$binsfeld_hit &
    !pvalue_scatter_long$modeled_hit & !pvalue_scatter_long$evc_huber_hit
  pvalue_scatter_long$dstressr_unique <- ifelse(
    pvalue_scatter_long$comparison == "DStressR without EV",
    !pvalue_scatter_long$binsfeld_hit & pvalue_scatter_long$modeled_hit & !pvalue_scatter_long$evc_huber_hit,
    !pvalue_scatter_long$binsfeld_hit & !pvalue_scatter_long$modeled_hit & pvalue_scatter_long$evc_huber_hit
  )
  pvalue_scatter_long$discordant_call <- ifelse(
    pvalue_scatter_long$binsfeld_only_all,
    "Binsfeld only",
    ifelse(pvalue_scatter_long$dstressr_unique, "DStressR unique", "Other")
  )
  pvalue_scatter_long$discordant_call <- factor(
    pvalue_scatter_long$discordant_call,
    levels = c("Binsfeld only", "DStressR unique", "Other")
  )
  pvalue_scatter_long$highlight <- pvalue_scatter_long$discordant_call != "Other"
  pvalue_scatter_long$label <- paste(pvalue_scatter_long$promoter, pvalue_scatter_long$compound, sep = "-")
  top_pvalue_scatter_labels <- do.call(rbind, lapply(split(pvalue_scatter_long, pvalue_scatter_long$comparison), function(d) {
    binsfeld_only_labels <- d[d$binsfeld_only_all & is.finite(d$disagreement), , drop = FALSE]
    binsfeld_only_labels <- binsfeld_only_labels[
      order(-binsfeld_only_labels$disagreement),
      ,
      drop = FALSE
    ]
    binsfeld_only_labels <- utils::head(binsfeld_only_labels, 16)
    dstressr_unique_labels <- d[d$dstressr_unique & is.finite(d$disagreement), , drop = FALSE]
    dstressr_unique_labels <- dstressr_unique_labels[
      order(-dstressr_unique_labels$disagreement),
      ,
      drop = FALSE
    ]
    rbind(binsfeld_only_labels, dstressr_unique_labels)
  }))
  pvalue_scatter_long$label_hjust <- ifelse(
    pvalue_scatter_long$binsfeld_pvalue < 0.15,
    -0.02,
    ifelse(
      pvalue_scatter_long$destress_pvalue > pvalue_scatter_long$binsfeld_pvalue,
      1.02,
      -0.02
    )
  )
  top_pvalue_scatter_labels$label_hjust <- ifelse(
    top_pvalue_scatter_labels$binsfeld_pvalue < 0.15,
    -0.02,
    ifelse(
      top_pvalue_scatter_labels$destress_pvalue > top_pvalue_scatter_labels$binsfeld_pvalue,
      1.02,
      -0.02
    )
  )
  pvalue_scatter_plot <- ggplot2$ggplot() +
    ggplot2$geom_abline(slope = 1, intercept = 0, color = "#111827", linewidth = 0.3, linetype = "dashed") +
    ggplot2$geom_point(
      data = pvalue_scatter_long,
      ggplot2$aes(binsfeld_pvalue, destress_pvalue),
      color = "#d1d5db",
      alpha = 0.56,
      size = 0.72
    ) +
    ggplot2$geom_point(
      data = pvalue_scatter_long[pvalue_scatter_long$highlight, , drop = FALSE],
      ggplot2$aes(
        binsfeld_pvalue,
        destress_pvalue,
        color = promoter,
        shape = discordant_call
      ),
      alpha = 0.9,
      size = 1.45
    ) +
    ggrepel$geom_text_repel(
      data = top_pvalue_scatter_labels,
      ggplot2$aes(
        binsfeld_pvalue,
        destress_pvalue,
        label = label,
        hjust = label_hjust,
        color = promoter
      ),
      size = 1.7,
      min.segment.length = 0,
      segment.color = "#6b7280",
      segment.size = 0.13,
      segment.alpha = 0.65,
      box.padding = 0.10,
      point.padding = 0.08,
      force = 1.6,
      force_pull = 0.12,
      max.overlaps = Inf,
      max.time = 3,
      max.iter = 8000,
      seed = 3262,
      show.legend = FALSE
    ) +
    ggplot2$facet_wrap(ggplot2$vars(comparison), nrow = 1) +
    ggplot2$coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2$scale_x_continuous(breaks = seq(0, 1, by = 0.25)) +
    ggplot2$scale_y_continuous(breaks = seq(0, 1, by = 0.25)) +
    ggplot2$scale_color_manual(values = promoter_colors[names(promoter_colors) != "Not called"]) +
    ggplot2$scale_shape_manual(values = c("Binsfeld only" = 17, "DStressR unique" = 16)) +
    ggplot2$theme_light(base_size = 10) +
    ggplot2$theme(
      panel.grid.minor = ggplot2$element_blank(),
      strip.text = ggplot2$element_text(face = "bold", hjust = 0.5, size = 9),
      legend.position = "none",
      plot.margin = ggplot2$margin(4, 4, 4, 4)
    ) +
    ggplot2$labs(
      x = "Binsfeld reference raw p-value",
      y = "DStressR raw p-value"
    )

  discovery_figure <- gridExtra$arrangeGrob(
    gridExtra$arrangeGrob(panel_label("a"), panel_label("b"), panel_label("c"), ncol = 3),
    pvalue_hist_all,
    gridExtra$arrangeGrob(panel_label("d"), panel_label("e"), panel_label("f"), ncol = 3),
    volcano_overlap_plot,
    gridExtra$arrangeGrob(
      panel_label("g"),
      gridExtra$arrangeGrob(panel_label("h"), panel_label("i"), ncol = 2),
      ncol = 2,
      widths = c(0.82, 1.18)
    ),
    gridExtra$arrangeGrob(
      venn_compact,
      pvalue_scatter_plot,
      ncol = 2,
      widths = c(0.82, 1.18)
    ),
    ncol = 1,
    heights = c(0.040, 0.85, 0.040, 1.08, 0.040, 1.08)
  )

  ggplot2$ggsave(file.path(out_dir, "evc_huber_hit_overlap_venn.png"), venn, width = 7.2, height = 6.8, dpi = 220)
  ggplot2$ggsave(file.path(out_dir, "evc_huber_hit_overlap_venn.pdf"), venn, width = 7.2, height = 6.8)
  ggplot2$ggsave(file.path(out_dir, "evc_huber_pvalue_histograms.png"), pvalue_hist, width = 11, height = 10.5, dpi = 220)
  ggplot2$ggsave(file.path(out_dir, "evc_huber_pvalue_histograms.pdf"), pvalue_hist, width = 11, height = 10.5)
  ggplot2$ggsave(file.path(out_dir, "evc_huber_pvalue_histograms_by_promoter.png"), pvalue_hist_by_promoter, width = 11, height = 9.4, dpi = 220)
  ggplot2$ggsave(file.path(out_dir, "evc_huber_pvalue_histograms_by_promoter.pdf"), pvalue_hist_by_promoter, width = 11, height = 9.4)
  ggplot2$ggsave(file.path(out_dir, "ecoli_reference_volcano_plot.png"), reference_plot, width = 8.2, height = 5.6, dpi = 220)
  ggplot2$ggsave(file.path(out_dir, "ecoli_reference_volcano_plot.pdf"), reference_plot, width = 8.2, height = 5.6)
  ggplot2$ggsave(file.path(out_dir, "evc_huber_volcano_plots.png"), volcano_plot, width = 12, height = 4.8, dpi = 220)
  ggplot2$ggsave(file.path(out_dir, "evc_huber_volcano_plots.pdf"), volcano_plot, width = 12, height = 4.8)
  ggplot2$ggsave(file.path(out_dir, "evc_huber_intersection_volcano_plots.png"), intersection_plot, width = 12, height = 5.2, dpi = 220)
  ggplot2$ggsave(file.path(out_dir, "evc_huber_intersection_volcano_plots.pdf"), intersection_plot, width = 12, height = 5.2)
  ggplot2$ggsave(file.path(out_dir, "evc_huber_nonintersection_volcano_plots.png"), non_intersection_plot, width = 12, height = 5.4, dpi = 220)
  ggplot2$ggsave(file.path(out_dir, "evc_huber_nonintersection_volcano_plots.pdf"), non_intersection_plot, width = 12, height = 5.4)
  ggplot2$ggsave(file.path(out_dir, "evc_huber_pvalue_histograms_all_promoters.png"), pvalue_hist_all, width = 11, height = 2.8, dpi = 220)
  ggplot2$ggsave(file.path(out_dir, "evc_huber_pvalue_histograms_all_promoters.pdf"), pvalue_hist_all, width = 11, height = 2.8)
  ggplot2$ggsave(file.path(out_dir, "evc_huber_pvalue_scatter_vs_binsfeld.png"), pvalue_scatter_plot, width = 8.2, height = 4.1, dpi = 220)
  ggplot2$ggsave(file.path(out_dir, "evc_huber_pvalue_scatter_vs_binsfeld.pdf"), pvalue_scatter_plot, width = 8.2, height = 4.1)
  ggplot2$ggsave(file.path(out_dir, "evc_huber_discovery_summary_figure.png"), discovery_figure, width = 12.5, height = 10.8, dpi = 220)
  ggplot2$ggsave(file.path(out_dir, "evc_huber_discovery_summary_figure.pdf"), discovery_figure, width = 12.5, height = 10.8)

  print(summary)
}

message("Wrote E. coli DStressR-with-EV analysis to: ", out_dir)
