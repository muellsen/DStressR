source(file.path("analysis", "_helpers.R"))
load_destress_package()

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package `ggplot2` is required for Binsfeld comparison plots.", call. = FALSE)
}

ggplot2 <- asNamespace("ggplot2")

load(analysis_path("data", "binsfeld_reporter_data.rda"))
out_dir <- analysis_output_dir("binsfeld")

hit_key <- function(x) paste(x$promoter, x$compound, sep = "::")
safe_neglog10 <- function(x) -log10(pmax(as.numeric(x), .Machine$double.xmin))

wt_auc <- binsfeld_reporter_auc[
  binsfeld_reporter_auc$strain == "WT" &
    binsfeld_reporter_auc$removed == "No",
]
wt_auc_model <- wt_auc[wt_auc$promoter != "EVC", ]

assay <- prepare_assay(
  wt_auc_model,
  promoter = "promoter",
  compound = "compound",
  control = "Water",
  lux = "lux_auc",
  growth = "od_auc",
  growth_exponent = "estimate",
  batch = "dose_level",
  replicate = "replicate",
  growth_covariates = "replicate",
  numeric_covariates = "dose_level"
)

fit <- fit_destress(
  assay,
  preset = "model",
  technical = c("replicate", "dose_level"),
  empirical_bayes = TRUE,
  adjustment = "by_promoter",
  interaction = FALSE
)

growth_exponents <- model_parameters(fit)$growth_exponents
if (!is.null(growth_exponents) && nrow(growth_exponents) > 0) {
  utils::write.table(
    growth_exponents,
    file.path(out_dir, "destress_default_growth_exponents.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
}

destress <- results(fit)
destress <- destress[destress$compound != "Water", ]
destress$destress_hit_class <- call_hits(
  destress,
  fdr = 0.05,
  effect = "specific_effect",
  padj = "specific_padj_by_promoter"
)$hit
destress$destress_hit <- destress$destress_hit_class != "Not DE"
destress$pair_id <- hit_key(destress)

wt_z <- binsfeld_reporter_scores[
  binsfeld_reporter_scores$strain == "WT" &
    binsfeld_reporter_scores$statistic == "Z_scores" &
    binsfeld_reporter_scores$promoter != "EVC",
]

author <- do.call(rbind, lapply(sort(unique(wt_z$promoter)), function(promoter) {
  promoter_z <- wt_z[wt_z$promoter == promoter, ]
  water <- promoter_z$value[grepl("^Water", promoter_z$drug)]
  rows <- lapply(sort(unique(promoter_z$drug)), function(drug) {
    values <- promoter_z$value[promoter_z$drug == drug]
    pvalue <- tryCatch(
      stats::wilcox.test(values, water)$p.value,
      error = function(e) NA_real_
    )
    data.frame(
      promoter = promoter,
      compound = ifelse(grepl("^Water_", drug), "Water", drug),
      mean_z = mean(values, na.rm = TRUE),
      binsfeld_pvalue = pvalue,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out$binsfeld_padj <- stats::p.adjust(out$binsfeld_pvalue, method = "BH")
  out
}))
author <- author[author$compound != "Water", ]
author$binsfeld_hit <- is.finite(author$binsfeld_padj) &
  author$binsfeld_padj < 0.05 &
  abs(author$mean_z) > 1
author$binsfeld_direction <- ifelse(
  !author$binsfeld_hit,
  "Not significant",
  ifelse(author$mean_z > 0, "Positive", "Negative")
)
author$pair_id <- hit_key(author)

comparison <- merge(
  author[, c(
    "promoter", "compound", "pair_id", "mean_z", "binsfeld_pvalue",
    "binsfeld_padj", "binsfeld_hit", "binsfeld_direction"
  )],
  destress[, c(
    "promoter", "compound", "pair_id", "specific_effect", "specific_pvalue",
    "specific_padj_by_promoter", "destress_hit", "destress_hit_class"
  )],
  by = c("promoter", "compound", "pair_id"),
  all = TRUE,
  sort = FALSE
)
comparison$overlap_class <- ifelse(
  comparison$binsfeld_hit & comparison$destress_hit,
  "Both",
  ifelse(
    comparison$binsfeld_hit,
    "Binsfeld only",
    ifelse(comparison$destress_hit, "DStressR only", "Neither")
  )
)

utils::write.table(
  comparison[order(comparison$promoter, comparison$compound), ],
  file.path(out_dir, "binsfeld_destress_all_pair_comparison.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

significant_union <- comparison[comparison$overlap_class != "Neither", c(
  "promoter", "compound", "overlap_class", "mean_z", "binsfeld_pvalue",
  "binsfeld_padj", "binsfeld_direction", "specific_effect", "specific_pvalue",
  "specific_padj_by_promoter", "destress_hit_class"
)]
significant_union <- significant_union[order(
  factor(significant_union$overlap_class, levels = c("Both", "Binsfeld only", "DStressR only")),
  significant_union$promoter,
  significant_union$compound
), ]
utils::write.table(
  significant_union,
  file.path(out_dir, "binsfeld_destress_significant_union.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
utils::write.table(
  significant_union[significant_union$overlap_class == "Binsfeld only", ],
  file.path(out_dir, "binsfeld_only_significant_pairs.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
utils::write.table(
  significant_union[significant_union$overlap_class == "DStressR only", ],
  file.path(out_dir, "destress_only_significant_pairs.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
utils::write.table(
  significant_union[significant_union$overlap_class == "Both", ],
  file.path(out_dir, "overlapping_significant_pairs.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

counts <- c(
  binsfeld_only = sum(comparison$overlap_class == "Binsfeld only"),
  overlap = sum(comparison$overlap_class == "Both"),
  destress_only = sum(comparison$overlap_class == "DStressR only")
)
binsfeld_hit_n <- sum(comparison$binsfeld_hit)
destress_hit_n <- sum(comparison$destress_hit)

venn_df <- data.frame(
  x = c(0, 1.1),
  y = c(0, 0),
  r = c(1, 1),
  method = c("Binsfeld", "DStressR")
)
circle_points <- do.call(rbind, lapply(seq_len(nrow(venn_df)), function(i) {
  theta <- seq(0, 2 * pi, length.out = 361)
  data.frame(
    x = venn_df$x[i] + venn_df$r[i] * cos(theta),
    y = venn_df$y[i] + venn_df$r[i] * sin(theta),
    method = venn_df$method[i],
    stringsAsFactors = FALSE
  )
}))

venn <- ggplot2$ggplot(circle_points, ggplot2$aes(x, y, fill = method, color = method)) +
  ggplot2$geom_polygon(alpha = 0.28, linewidth = 0.7) +
  ggplot2$annotate("text", x = -0.45, y = 0, label = counts[["binsfeld_only"]], size = 8, fontface = "bold") +
  ggplot2$annotate("text", x = 0.55, y = 0, label = counts[["overlap"]], size = 8, fontface = "bold") +
  ggplot2$annotate("text", x = 1.55, y = 0, label = counts[["destress_only"]], size = 8, fontface = "bold") +
  ggplot2$annotate("text", x = -0.38, y = 1.15, label = paste0("Binsfeld\n", binsfeld_hit_n, " hits"), size = 4.2, fontface = "bold") +
  ggplot2$annotate("text", x = 1.48, y = 1.15, label = paste0("DStressR\n", destress_hit_n, " hits"), size = 4.2, fontface = "bold") +
  ggplot2$scale_fill_manual(values = c("Binsfeld" = "#2563eb", "DStressR" = "#dc2626")) +
  ggplot2$scale_color_manual(values = c("Binsfeld" = "#1d4ed8", "DStressR" = "#b91c1c")) +
  ggplot2$coord_equal(xlim = c(-1.15, 2.25), ylim = c(-1.05, 1.35), expand = FALSE) +
  ggplot2$theme_void(base_size = 11) +
  ggplot2$theme(legend.position = "none", plot.title = ggplot2$element_text(face = "bold")) +
  ggplot2$labs(
    title = "WT reporter-screen hit overlap",
    subtitle = "Promoter-compound pairs; original Wilcoxon/Z-score rule vs DStressR model"
  )

pvalue_df <- rbind(
  data.frame(
    method = "Binsfeld Wilcoxon/Z-score",
    promoter = author$promoter,
    raw_pvalue = author$binsfeld_pvalue,
    adjusted_pvalue = author$binsfeld_padj,
    hit = author$binsfeld_hit,
    stringsAsFactors = FALSE
  ),
  data.frame(
    method = "DStressR default modeled response",
    promoter = destress$promoter,
    raw_pvalue = destress$specific_pvalue,
    adjusted_pvalue = destress$specific_padj_by_promoter,
    hit = destress$destress_hit,
    stringsAsFactors = FALSE
  )
)
pvalue_long <- rbind(
  data.frame(
    method = pvalue_df$method,
    promoter = pvalue_df$promoter,
    pvalue_type = "Raw p-value",
    pvalue = pvalue_df$raw_pvalue,
    hit = pvalue_df$hit,
    stringsAsFactors = FALSE
  ),
  data.frame(
    method = pvalue_df$method,
    promoter = pvalue_df$promoter,
    pvalue_type = "Promoter-wise BH adjusted",
    pvalue = pvalue_df$adjusted_pvalue,
    hit = pvalue_df$hit,
    stringsAsFactors = FALSE
  )
)
pvalue_long$method <- factor(pvalue_long$method, levels = c(
  "Binsfeld Wilcoxon/Z-score",
  "DStressR default modeled response"
))
pvalue_long$pvalue_type <- factor(pvalue_long$pvalue_type, levels = c(
  "Raw p-value",
  "Promoter-wise BH adjusted"
))

pvalue_hist <- ggplot2$ggplot(pvalue_long, ggplot2$aes(pvalue)) +
  ggplot2$geom_histogram(
    data = pvalue_long[!pvalue_long$hit, ],
    bins = 30,
    boundary = 0,
    fill = "#d1d5db",
    color = "white",
    linewidth = 0.1
  ) +
  ggplot2$geom_histogram(
    data = pvalue_long[pvalue_long$hit, ],
    ggplot2$aes(fill = method),
    bins = 30,
    boundary = 0,
    color = "white",
    linewidth = 0.1
  ) +
  ggplot2$facet_grid(pvalue_type ~ method) +
  ggplot2$scale_x_continuous(limits = c(0, 1), expand = c(0.01, 0.01)) +
  ggplot2$scale_fill_manual(values = c(
    "Binsfeld Wilcoxon/Z-score" = "#2563eb",
    "DStressR default modeled response" = "#dc2626"
  )) +
  ggplot2$theme_light(base_size = 10) +
  ggplot2$theme(
    legend.position = "none",
    panel.grid.minor = ggplot2$element_blank(),
    plot.title = ggplot2$element_text(face = "bold")
  ) +
  ggplot2$labs(
    title = "P-value distributions for the WT reporter screen",
    subtitle = "Grey bars are non-hit pairs; colored bars are hit pairs under each method",
    x = "p-value",
    y = "Promoter-compound pairs"
  )

effect_df <- rbind(
  data.frame(
    method = "Binsfeld mean Z-score",
    promoter = author$promoter,
    value = author$mean_z,
    hit = author$binsfeld_hit,
    stringsAsFactors = FALSE
  ),
  data.frame(
    method = "DStressR specific effect",
    promoter = destress$promoter,
    value = destress$specific_effect,
    hit = destress$destress_hit,
    stringsAsFactors = FALSE
  )
)
effect_df$method <- factor(effect_df$method, levels = c(
  "Binsfeld mean Z-score",
  "DStressR specific effect"
))

effect_hist <- ggplot2$ggplot(effect_df, ggplot2$aes(value)) +
  ggplot2$geom_histogram(
    data = effect_df[!effect_df$hit, ],
    bins = 45,
    fill = "#d1d5db",
    color = "white",
    linewidth = 0.1
  ) +
  ggplot2$geom_histogram(
    data = effect_df[effect_df$hit, ],
    ggplot2$aes(fill = method),
    bins = 45,
    color = "white",
    linewidth = 0.1
  ) +
  ggplot2$geom_vline(xintercept = 0, color = "#111827", linewidth = 0.35) +
  ggplot2$facet_wrap(ggplot2$vars(method), scales = "free_x", ncol = 1) +
  ggplot2$scale_fill_manual(values = c(
    "Binsfeld mean Z-score" = "#2563eb",
    "DStressR specific effect" = "#dc2626"
  )) +
  ggplot2$theme_light(base_size = 10) +
  ggplot2$theme(
    legend.position = "none",
    panel.grid.minor = ggplot2$element_blank(),
    plot.title = ggplot2$element_text(face = "bold")
  ) +
  ggplot2$labs(
    title = "Effect-score distributions",
    subtitle = "Colored bars mark method-specific hit pairs",
    x = "Effect score",
    y = "Promoter-compound pairs"
  )

comparison$neglog10_binsfeld <- safe_neglog10(comparison$binsfeld_pvalue)
comparison$neglog10_destress <- safe_neglog10(comparison$specific_pvalue)
comparison$overlap_class <- factor(
  comparison$overlap_class,
  levels = c("Neither", "Binsfeld only", "DStressR only", "Both")
)

scatter_effect <- ggplot2$ggplot(
  comparison,
  ggplot2$aes(mean_z, specific_effect, color = overlap_class)
) +
  ggplot2$geom_hline(yintercept = 0, color = "#6b7280", linewidth = 0.25) +
  ggplot2$geom_vline(xintercept = 0, color = "#6b7280", linewidth = 0.25) +
  ggplot2$geom_point(alpha = 0.82, size = 1.9) +
  ggplot2$scale_color_manual(values = c(
    "Neither" = "#9ca3af",
    "Binsfeld only" = "#2563eb",
    "DStressR only" = "#dc2626",
    "Both" = "#16a34a"
  )) +
  ggplot2$theme_light(base_size = 10) +
  ggplot2$theme(
    panel.grid.minor = ggplot2$element_blank(),
    plot.title = ggplot2$element_text(face = "bold"),
    legend.title = ggplot2$element_blank()
  ) +
  ggplot2$labs(
    title = "Effect comparison across all WT promoter-compound pairs",
    x = "Binsfeld mean Z-score",
    y = "DStressR specific effect"
  )

scatter_pvalue <- ggplot2$ggplot(
  comparison,
  ggplot2$aes(neglog10_binsfeld, neglog10_destress, color = overlap_class)
) +
  ggplot2$geom_point(alpha = 0.82, size = 1.9) +
  ggplot2$scale_color_manual(values = c(
    "Neither" = "#9ca3af",
    "Binsfeld only" = "#2563eb",
    "DStressR only" = "#dc2626",
    "Both" = "#16a34a"
  )) +
  ggplot2$theme_light(base_size = 10) +
  ggplot2$theme(
    panel.grid.minor = ggplot2$element_blank(),
    plot.title = ggplot2$element_text(face = "bold"),
    legend.title = ggplot2$element_blank()
  ) +
  ggplot2$labs(
    title = "Raw p-value comparison",
    x = "Binsfeld -log10 raw p-value",
    y = "DStressR -log10 raw p-value"
  )

scatter_pvalue_zoom <- scatter_pvalue +
  ggplot2$coord_cartesian(xlim = c(0, 6), ylim = c(0, 10)) +
  ggplot2$labs(
    title = "Raw p-value comparison, zoomed",
    subtitle = "Axes clipped to emphasize the main cloud; full plot is also saved"
  )

save_plot <- function(name, plot, width, height) {
  ggplot2$ggsave(file.path(out_dir, paste0(name, ".png")), plot, width = width, height = height, dpi = 300, bg = "white")
  ggplot2$ggsave(file.path(out_dir, paste0(name, ".pdf")), plot, width = width, height = height, bg = "white")
}

save_plot("hit_overlap_venn", venn, 6.5, 4.6)
save_plot("pvalue_histograms", pvalue_hist, 9.5, 6.5)
save_plot("effect_histograms", effect_hist, 8.5, 6.2)
save_plot("effect_scatter", scatter_effect, 7.2, 5.4)
save_plot("pvalue_scatter", scatter_pvalue, 7.2, 5.4)
save_plot("pvalue_scatter_zoom", scatter_pvalue_zoom, 7.2, 5.4)

plot_summary <- data.frame(
  figure = c(
    "hit_overlap_venn",
    "pvalue_histograms",
    "effect_histograms",
    "effect_scatter",
    "pvalue_scatter",
    "pvalue_scatter_zoom"
  ),
  png = file.path(out_dir, paste0(c(
    "hit_overlap_venn",
    "pvalue_histograms",
    "effect_histograms",
    "effect_scatter",
    "pvalue_scatter",
    "pvalue_scatter_zoom"
  ), ".png")),
  stringsAsFactors = FALSE
)
utils::write.table(
  plot_summary,
  file.path(out_dir, "binsfeld_plot_manifest.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

print(data.frame(
  metric = c("Binsfeld hits", "DStressR hits", "Overlap"),
  value = c(
    sum(comparison$binsfeld_hit),
    sum(comparison$destress_hit),
    sum(comparison$overlap_class == "Both")
  )
))
message("Wrote E. coli comparison plots to: ", out_dir)
