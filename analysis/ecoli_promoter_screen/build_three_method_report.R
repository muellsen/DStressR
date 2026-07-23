source(file.path("analysis", "_helpers.R"))
load_destress_package()

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package `ggplot2` is required for Binsfeld comparison plots.", call. = FALSE)
}
if (!requireNamespace("base64enc", quietly = TRUE)) {
  stop("Package `base64enc` is required to build the standalone HTML report.", call. = FALSE)
}

ggplot2 <- asNamespace("ggplot2")

load(analysis_path("data", "binsfeld_reporter_data.rda"))
out_dir <- analysis_output_dir("binsfeld_three_method")

hit_key <- function(x) paste(x$promoter, x$compound, sep = "::")
safe_neglog10 <- function(x) -log10(pmax(as.numeric(x), .Machine$double.xmin))
fmt_num <- function(x, digits = 3) {
  ifelse(is.na(x), "", format(signif(as.numeric(x), digits), scientific = TRUE, trim = TRUE))
}
html_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub('"', "&quot;", x, fixed = TRUE)
  x
}
img_uri <- function(path) {
  paste0("data:image/png;base64,", base64enc::base64encode(path))
}
html_table <- function(df, columns = names(df), max_rows = Inf) {
  d <- df[, columns, drop = FALSE]
  if (is.finite(max_rows) && nrow(d) > max_rows) {
    d <- d[seq_len(max_rows), , drop = FALSE]
  }
  header <- paste0("<th>", html_escape(names(d)), "</th>", collapse = "")
  rows <- apply(d, 1, function(row) {
    paste0("<tr>", paste0("<td>", html_escape(row), "</td>", collapse = ""), "</tr>")
  })
  paste0("<table><thead><tr>", header, "</tr></thead><tbody>", paste(rows, collapse = "\n"), "</tbody></table>")
}

wt_auc <- binsfeld_reporter_auc[
  binsfeld_reporter_auc$strain == "WT" &
    binsfeld_reporter_auc$removed == "No",
]
wt_auc_model <- wt_auc[wt_auc$promoter != "EVC", ]

run_destress <- function(growth_exponent, label) {
  assay <- prepare_assay(
    wt_auc_model,
    promoter = "promoter",
    compound = "compound",
    control = "Water",
    lux = "lux_auc",
    growth = "od_auc",
    growth_exponent = growth_exponent,
    batch = "dose_level",
    replicate = "replicate"
  )

  fit <- fit_destress(
    assay,
    preset = "model",
    technical = c("replicate", "dose_level"),
    empirical_bayes = TRUE,
    adjustment = "by_promoter",
    interaction = FALSE
  )

  out <- results(fit)
  out <- out[out$compound != "Water", ]
  hit_class <- call_hits(
    out,
    fdr = 0.05,
    effect = "specific_effect",
    padj = "specific_padj_by_promoter"
  )$hit
  out[[paste0(label, "_hit_class")]] <- hit_class
  out[[paste0(label, "_hit")]] <- hit_class != "Not DE"
  out$pair_id <- hit_key(out)
  keep <- c(
    "promoter", "compound", "pair_id", "specific_effect",
    "specific_pvalue", "specific_padj_by_promoter",
    paste0(label, "_hit"), paste0(label, "_hit_class")
  )
  out <- out[, keep]
  names(out)[names(out) == "specific_effect"] <- paste0(label, "_effect")
  names(out)[names(out) == "specific_pvalue"] <- paste0(label, "_pvalue")
  names(out)[names(out) == "specific_padj_by_promoter"] <- paste0(label, "_padj_by_promoter")

  list(results = out, parameters = model_parameters(fit))
}

modeled_fit <- run_destress("estimate", "modeled")
standard_fit <- run_destress(1, "standard")

growth_exponents <- modeled_fit$parameters$growth_exponents
if (!is.null(growth_exponents) && nrow(growth_exponents) > 0) {
  utils::write.table(
    growth_exponents,
    file.path(out_dir, "destress_modeled_growth_exponents.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
}

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
  modeled_fit$results,
  by = c("promoter", "compound", "pair_id"),
  all = TRUE,
  sort = FALSE
)
comparison <- merge(
  comparison,
  standard_fit$results,
  by = c("promoter", "compound", "pair_id"),
  all = TRUE,
  sort = FALSE
)

comparison$three_method_class <- apply(
  comparison[, c("binsfeld_hit", "modeled_hit", "standard_hit")],
  1,
  function(x) {
    names <- c("Binsfeld", "DStressR modeled", "DStressR alpha=1")[as.logical(x)]
    if (length(names) == 0) {
      "None"
    } else {
      paste(names, collapse = " + ")
    }
  }
)

comparison <- comparison[order(comparison$promoter, comparison$compound), ]
utils::write.table(
  comparison,
  file.path(out_dir, "binsfeld_three_method_all_pair_comparison.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

significant_union <- comparison[comparison$three_method_class != "None", c(
  "promoter", "compound", "three_method_class",
  "mean_z", "binsfeld_pvalue", "binsfeld_padj", "binsfeld_direction",
  "modeled_effect", "modeled_pvalue", "modeled_padj_by_promoter", "modeled_hit_class",
  "standard_effect", "standard_pvalue", "standard_padj_by_promoter", "standard_hit_class"
)]
class_levels <- c(
  "Binsfeld + DStressR modeled + DStressR alpha=1",
  "Binsfeld + DStressR modeled",
  "Binsfeld + DStressR alpha=1",
  "DStressR modeled + DStressR alpha=1",
  "Binsfeld",
  "DStressR modeled",
  "DStressR alpha=1"
)
significant_union <- significant_union[order(
  factor(significant_union$three_method_class, levels = class_levels),
  significant_union$promoter,
  significant_union$compound
), ]
utils::write.table(
  significant_union,
  file.path(out_dir, "binsfeld_three_method_significant_union.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

for (class_name in class_levels) {
  slug <- tolower(gsub("[^A-Za-z0-9]+", "_", class_name))
  slug <- gsub("^_|_$", "", slug)
  utils::write.table(
    significant_union[significant_union$three_method_class == class_name, ],
    file.path(out_dir, paste0(slug, "_significant_pairs.tsv")),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
}

set_counts <- c(
  binsfeld = sum(comparison$binsfeld_hit),
  modeled = sum(comparison$modeled_hit),
  standard = sum(comparison$standard_hit)
)
region_counts <- c(
  binsfeld_only = sum(comparison$binsfeld_hit & !comparison$modeled_hit & !comparison$standard_hit),
  modeled_only = sum(!comparison$binsfeld_hit & comparison$modeled_hit & !comparison$standard_hit),
  standard_only = sum(!comparison$binsfeld_hit & !comparison$modeled_hit & comparison$standard_hit),
  binsfeld_modeled = sum(comparison$binsfeld_hit & comparison$modeled_hit & !comparison$standard_hit),
  binsfeld_standard = sum(comparison$binsfeld_hit & !comparison$modeled_hit & comparison$standard_hit),
  modeled_standard = sum(!comparison$binsfeld_hit & comparison$modeled_hit & comparison$standard_hit),
  all_three = sum(comparison$binsfeld_hit & comparison$modeled_hit & comparison$standard_hit)
)

count_summary <- data.frame(
  Metric = c(
    "Promoter-compound pairs tested",
    "Binsfeld-style significant pairs",
    "DStressR modeled-response significant pairs",
    "DStressR alpha=1 significant pairs",
    "All three",
    "Binsfeld + modeled only",
    "Binsfeld + alpha=1 only",
    "Modeled + alpha=1 only",
    "Binsfeld only",
    "Modeled only",
    "Alpha=1 only",
    "Union significant by at least one method",
    "None"
  ),
  Count = c(
    nrow(comparison),
    set_counts[["binsfeld"]],
    set_counts[["modeled"]],
    set_counts[["standard"]],
    region_counts[["all_three"]],
    region_counts[["binsfeld_modeled"]],
    region_counts[["binsfeld_standard"]],
    region_counts[["modeled_standard"]],
    region_counts[["binsfeld_only"]],
    region_counts[["modeled_only"]],
    region_counts[["standard_only"]],
    nrow(significant_union),
    sum(comparison$three_method_class == "None")
  ),
  stringsAsFactors = FALSE
)
utils::write.table(
  count_summary,
  file.path(out_dir, "binsfeld_three_method_hit_summary.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

circle_points <- do.call(rbind, lapply(seq_len(3), function(i) {
  centers <- data.frame(
    x = c(0, 1.25, 0.62),
    y = c(0.3, 0.3, -0.78),
    method = c("Binsfeld", "DStressR modeled", "DStressR alpha=1")
  )
  theta <- seq(0, 2 * pi, length.out = 361)
  data.frame(
    x = centers$x[i] + cos(theta),
    y = centers$y[i] + sin(theta),
    method = centers$method[i],
    stringsAsFactors = FALSE
  )
}))

venn <- ggplot2$ggplot(circle_points, ggplot2$aes(x, y, fill = method, color = method)) +
  ggplot2$geom_polygon(alpha = 0.22, linewidth = 0.75) +
  ggplot2$annotate("text", x = -0.45, y = 0.45, label = region_counts[["binsfeld_only"]], size = 6.2, fontface = "bold") +
  ggplot2$annotate("text", x = 1.72, y = 0.45, label = region_counts[["modeled_only"]], size = 6.2, fontface = "bold") +
  ggplot2$annotate("text", x = 0.62, y = -1.22, label = region_counts[["standard_only"]], size = 6.2, fontface = "bold") +
  ggplot2$annotate("text", x = 0.62, y = 0.5, label = region_counts[["binsfeld_modeled"]], size = 6.2, fontface = "bold") +
  ggplot2$annotate("text", x = 0.18, y = -0.35, label = region_counts[["binsfeld_standard"]], size = 6.2, fontface = "bold") +
  ggplot2$annotate("text", x = 1.07, y = -0.35, label = region_counts[["modeled_standard"]], size = 6.2, fontface = "bold") +
  ggplot2$annotate("text", x = 0.62, y = -0.02, label = region_counts[["all_three"]], size = 6.8, fontface = "bold") +
  ggplot2$annotate("text", x = -0.42, y = 1.45, label = paste0("Binsfeld\n", set_counts[["binsfeld"]], " hits"), size = 3.8, fontface = "bold") +
  ggplot2$annotate("text", x = 1.67, y = 1.45, label = paste0("Modeled\n", set_counts[["modeled"]], " hits"), size = 3.8, fontface = "bold") +
  ggplot2$annotate("text", x = 0.62, y = -2.03, label = paste0("Alpha=1\n", set_counts[["standard"]], " hits"), size = 3.8, fontface = "bold") +
  ggplot2$scale_fill_manual(values = c(
    "Binsfeld" = "#2563eb",
    "DStressR modeled" = "#dc2626",
    "DStressR alpha=1" = "#059669"
  )) +
  ggplot2$scale_color_manual(values = c(
    "Binsfeld" = "#1d4ed8",
    "DStressR modeled" = "#b91c1c",
    "DStressR alpha=1" = "#047857"
  )) +
  ggplot2$coord_equal(xlim = c(-1.2, 2.45), ylim = c(-2.25, 1.7), expand = FALSE) +
  ggplot2$theme_void(base_size = 11) +
  ggplot2$theme(legend.position = "none", plot.title = ggplot2$element_text(face = "bold")) +
  ggplot2$labs(
    title = "Three-method WT reporter-screen hit overlap",
    subtitle = "Promoter-compound pairs called by Binsfeld-style, DStressR modeled-response, or DStressR alpha=1 analyses"
  )

pvalue_long <- rbind(
  data.frame(
    method = "Binsfeld Wilcoxon/Z-score",
    pvalue_type = "Raw p-value",
    pvalue = comparison$binsfeld_pvalue,
    hit = comparison$binsfeld_hit,
    stringsAsFactors = FALSE
  ),
  data.frame(
    method = "Binsfeld Wilcoxon/Z-score",
    pvalue_type = "Promoter-wise BH adjusted",
    pvalue = comparison$binsfeld_padj,
    hit = comparison$binsfeld_hit,
    stringsAsFactors = FALSE
  ),
  data.frame(
    method = "DStressR modeled response",
    pvalue_type = "Raw p-value",
    pvalue = comparison$modeled_pvalue,
    hit = comparison$modeled_hit,
    stringsAsFactors = FALSE
  ),
  data.frame(
    method = "DStressR modeled response",
    pvalue_type = "Promoter-wise BH adjusted",
    pvalue = comparison$modeled_padj_by_promoter,
    hit = comparison$modeled_hit,
    stringsAsFactors = FALSE
  ),
  data.frame(
    method = "DStressR alpha=1 response",
    pvalue_type = "Raw p-value",
    pvalue = comparison$standard_pvalue,
    hit = comparison$standard_hit,
    stringsAsFactors = FALSE
  ),
  data.frame(
    method = "DStressR alpha=1 response",
    pvalue_type = "Promoter-wise BH adjusted",
    pvalue = comparison$standard_padj_by_promoter,
    hit = comparison$standard_hit,
    stringsAsFactors = FALSE
  )
)
pvalue_long$method <- factor(pvalue_long$method, levels = c(
  "Binsfeld Wilcoxon/Z-score",
  "DStressR modeled response",
  "DStressR alpha=1 response"
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
    "DStressR modeled response" = "#dc2626",
    "DStressR alpha=1 response" = "#059669"
  )) +
  ggplot2$theme_light(base_size = 10) +
  ggplot2$theme(
    legend.position = "none",
    panel.grid.minor = ggplot2$element_blank(),
    plot.title = ggplot2$element_text(face = "bold")
  ) +
  ggplot2$labs(
    title = "P-value distributions for the WT reporter screen",
    subtitle = "Grey bars are non-hit pairs; colored bars are method-specific hits",
    x = "p-value",
    y = "Promoter-compound pairs"
  )

effect_long <- rbind(
  data.frame(
    method = "Binsfeld mean Z-score",
    value = comparison$mean_z,
    hit = comparison$binsfeld_hit,
    stringsAsFactors = FALSE
  ),
  data.frame(
    method = "DStressR modeled specific effect",
    value = comparison$modeled_effect,
    hit = comparison$modeled_hit,
    stringsAsFactors = FALSE
  ),
  data.frame(
    method = "DStressR alpha=1 specific effect",
    value = comparison$standard_effect,
    hit = comparison$standard_hit,
    stringsAsFactors = FALSE
  )
)
effect_long$method <- factor(effect_long$method, levels = c(
  "Binsfeld mean Z-score",
  "DStressR modeled specific effect",
  "DStressR alpha=1 specific effect"
))
effect_hist <- ggplot2$ggplot(effect_long, ggplot2$aes(value)) +
  ggplot2$geom_histogram(
    data = effect_long[!effect_long$hit, ],
    bins = 45,
    fill = "#d1d5db",
    color = "white",
    linewidth = 0.1
  ) +
  ggplot2$geom_histogram(
    data = effect_long[effect_long$hit, ],
    ggplot2$aes(fill = method),
    bins = 45,
    color = "white",
    linewidth = 0.1
  ) +
  ggplot2$geom_vline(xintercept = 0, color = "#111827", linewidth = 0.35) +
  ggplot2$facet_wrap(ggplot2$vars(method), scales = "free_x", ncol = 1) +
  ggplot2$scale_fill_manual(values = c(
    "Binsfeld mean Z-score" = "#2563eb",
    "DStressR modeled specific effect" = "#dc2626",
    "DStressR alpha=1 specific effect" = "#059669"
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
comparison$neglog10_modeled <- safe_neglog10(comparison$modeled_pvalue)
comparison$neglog10_standard <- safe_neglog10(comparison$standard_pvalue)

effect_reference <- rbind(
  data.frame(
    comparison = "DStressR modeled response",
    binsfeld_mean_z = comparison$mean_z,
    destress_effect = comparison$modeled_effect,
    hit_class = comparison$three_method_class,
    stringsAsFactors = FALSE
  ),
  data.frame(
    comparison = "DStressR alpha=1 response",
    binsfeld_mean_z = comparison$mean_z,
    destress_effect = comparison$standard_effect,
    hit_class = comparison$three_method_class,
    stringsAsFactors = FALSE
  )
)
effect_reference$comparison <- factor(effect_reference$comparison, levels = c(
  "DStressR modeled response",
  "DStressR alpha=1 response"
))

pvalue_reference <- rbind(
  data.frame(
    comparison = "DStressR modeled response",
    binsfeld_neglog10_pvalue = comparison$neglog10_binsfeld,
    destress_neglog10_pvalue = comparison$neglog10_modeled,
    hit_class = comparison$three_method_class,
    stringsAsFactors = FALSE
  ),
  data.frame(
    comparison = "DStressR alpha=1 response",
    binsfeld_neglog10_pvalue = comparison$neglog10_binsfeld,
    destress_neglog10_pvalue = comparison$neglog10_standard,
    hit_class = comparison$three_method_class,
    stringsAsFactors = FALSE
  )
)
pvalue_reference$comparison <- factor(pvalue_reference$comparison, levels = c(
  "DStressR modeled response",
  "DStressR alpha=1 response"
))

volcano_data <- rbind(
  data.frame(
    method = "Binsfeld-style",
    effect = comparison$mean_z,
    neglog10_pvalue = comparison$neglog10_binsfeld,
    hit = comparison$binsfeld_hit,
    stringsAsFactors = FALSE
  ),
  data.frame(
    method = "DStressR modeled response",
    effect = comparison$modeled_effect,
    neglog10_pvalue = comparison$neglog10_modeled,
    hit = comparison$modeled_hit,
    stringsAsFactors = FALSE
  )
)
volcano_data$method <- factor(volcano_data$method, levels = c(
  "Binsfeld-style",
  "DStressR modeled response"
))
volcano_data$status <- ifelse(volcano_data$hit, "Hit", "Not called")

volcano_plot <- ggplot2$ggplot(
  volcano_data,
  ggplot2$aes(effect, neglog10_pvalue, color = status)
) +
  ggplot2$geom_hline(yintercept = -log10(0.05), color = "#9ca3af", linewidth = 0.3, linetype = "dashed") +
  ggplot2$geom_vline(xintercept = 0, color = "#9ca3af", linewidth = 0.3) +
  ggplot2$geom_point(alpha = 0.78, size = 1.6) +
  ggplot2$facet_wrap(ggplot2$vars(method), scales = "free_x", ncol = 2) +
  ggplot2$scale_color_manual(values = c("Not called" = "#d1d5db", "Hit" = "#111827")) +
  ggplot2$theme_light(base_size = 10) +
  ggplot2$theme(
    panel.grid.minor = ggplot2$element_blank(),
    plot.title = ggplot2$element_text(face = "bold"),
    legend.position = "bottom"
  ) +
  ggplot2$labs(
    title = "Volcano plots for Binsfeld-style and modeled-response DStressR analyses",
    x = "Effect score",
    y = "-log10 raw p-value",
    color = "Call"
  )

reference_colors <- c(
  "None" = "#d1d5db",
  "Binsfeld" = "#2563eb",
  "DStressR modeled" = "#dc2626",
  "DStressR alpha=1" = "#059669",
  "Binsfeld + DStressR modeled" = "#7c3aed",
  "Binsfeld + DStressR alpha=1" = "#0891b2",
  "DStressR modeled + DStressR alpha=1" = "#ea580c",
  "Binsfeld + DStressR modeled + DStressR alpha=1" = "#111827"
)

response_scatter <- ggplot2$ggplot(
  effect_reference,
  ggplot2$aes(binsfeld_mean_z, destress_effect, color = hit_class)
) +
  ggplot2$geom_hline(yintercept = 0, color = "#9ca3af", linewidth = 0.3) +
  ggplot2$geom_vline(xintercept = 0, color = "#9ca3af", linewidth = 0.3) +
  ggplot2$geom_point(alpha = 0.8, size = 1.8) +
  ggplot2$facet_wrap(ggplot2$vars(comparison), ncol = 2) +
  ggplot2$scale_color_manual(values = reference_colors, drop = FALSE) +
  ggplot2$theme_light(base_size = 10) +
  ggplot2$theme(
    panel.grid.minor = ggplot2$element_blank(),
    plot.title = ggplot2$element_text(face = "bold"),
    legend.position = "bottom"
  ) +
  ggplot2$labs(
    title = "Binsfeld-referenced effect comparison",
    subtitle = "Each panel compares the Binsfeld-style mean Z-score with one DStressR response model",
    x = "Binsfeld-style mean Z-score",
    y = "DStressR specific effect",
    color = "Hit class"
  )

pvalue_scatter <- ggplot2$ggplot(
  pvalue_reference,
  ggplot2$aes(binsfeld_neglog10_pvalue, destress_neglog10_pvalue, color = hit_class)
) +
  ggplot2$geom_point(alpha = 0.8, size = 1.8) +
  ggplot2$facet_wrap(ggplot2$vars(comparison), ncol = 2) +
  ggplot2$scale_color_manual(values = reference_colors, drop = FALSE) +
  ggplot2$theme_light(base_size = 10) +
  ggplot2$theme(
    panel.grid.minor = ggplot2$element_blank(),
    plot.title = ggplot2$element_text(face = "bold"),
    legend.position = "bottom"
  ) +
  ggplot2$labs(
    title = "Binsfeld-referenced raw p-value comparison",
    subtitle = "Each panel compares the Binsfeld-style Wilcoxon p-value with one DStressR response model",
    x = "-log10 Binsfeld-style raw p-value",
    y = "-log10 DStressR raw p-value",
    color = "Hit class"
  )

ggplot2$ggsave(file.path(out_dir, "three_method_hit_overlap_venn.png"), venn, width = 7.2, height = 6.8, dpi = 220)
ggplot2$ggsave(file.path(out_dir, "three_method_pvalue_histograms.png"), pvalue_hist, width = 11, height = 6.2, dpi = 180)
ggplot2$ggsave(file.path(out_dir, "three_method_effect_histograms.png"), effect_hist, width = 9, height = 7, dpi = 180)
ggplot2$ggsave(file.path(out_dir, "binsfeld_modeled_volcano_plots.png"), volcano_plot, width = 10, height = 5.2, dpi = 220)
ggplot2$ggsave(file.path(out_dir, "binsfeld_modeled_volcano_plots.pdf"), volcano_plot, width = 10, height = 5.2)
ggplot2$ggsave(file.path(out_dir, "binsfeld_reference_effect_scatter.png"), response_scatter, width = 10, height = 5.8, dpi = 180)
ggplot2$ggsave(file.path(out_dir, "binsfeld_reference_pvalue_scatter.png"), pvalue_scatter, width = 10, height = 5.8, dpi = 180)

appendix <- significant_union
for (col in c(
  "mean_z", "binsfeld_pvalue", "binsfeld_padj", "modeled_effect",
  "modeled_pvalue", "modeled_padj_by_promoter", "standard_effect",
  "standard_pvalue", "standard_padj_by_promoter"
)) {
  appendix[[col]] <- fmt_num(appendix[[col]])
}

growth_table <- growth_exponents
for (col in intersect(c("a_raw", "a_raw_se", "alpha_raw", "alpha_raw_se", "alpha_shrunk"), names(growth_table))) {
  growth_table[[col]] <- fmt_num(growth_table[[col]])
}

pairwise_summary <- data.frame(
  Comparison = c(
    "Binsfeld vs DStressR modeled",
    "Binsfeld vs DStressR alpha=1",
    "DStressR modeled vs DStressR alpha=1"
  ),
  First = c(
    sum(comparison$binsfeld_hit & !comparison$modeled_hit),
    sum(comparison$binsfeld_hit & !comparison$standard_hit),
    sum(comparison$modeled_hit & !comparison$standard_hit)
  ),
  Overlap = c(
    sum(comparison$binsfeld_hit & comparison$modeled_hit),
    sum(comparison$binsfeld_hit & comparison$standard_hit),
    sum(comparison$modeled_hit & comparison$standard_hit)
  ),
  Second = c(
    sum(!comparison$binsfeld_hit & comparison$modeled_hit),
    sum(!comparison$binsfeld_hit & comparison$standard_hit),
    sum(!comparison$modeled_hit & comparison$standard_hit)
  ),
  stringsAsFactors = FALSE
)

dstressr_only_lit <- comparison[
  !comparison$binsfeld_hit & (comparison$modeled_hit | comparison$standard_hit),
  c(
    "promoter", "compound", "three_method_class",
    "modeled_effect", "modeled_padj_by_promoter",
    "standard_effect", "standard_padj_by_promoter"
  )
]
dstressr_only_lit$literature_support <- "Hypothesis-generating; no direct prior support identified in this quick review"
dstressr_only_lit$support_note <- "Not a Binsfeld-style hit; keep as candidate pending biological follow-up."
dstressr_only_lit$source <- "DStressR analysis"

set_support <- function(rows, category, note, source) {
  if (any(rows)) {
    dstressr_only_lit$literature_support[rows] <<- category
    dstressr_only_lit$support_note[rows] <<- note
    dstressr_only_lit$source[rows] <<- source
  }
}

set_support(
  dstressr_only_lit$promoter == "marRABp" &
    dstressr_only_lit$compound %in% c("Chloramphenicol", "Vanillin"),
  "Direct same-pair support",
  "Binsfeld et al. explicitly note these as known marRABp interactions missed by their stringent screen; DStressR recovers them.",
  "Binsfeld et al. 2025; https://doi.org/10.1371/journal.pbio.3003260"
)
set_support(
  dstressr_only_lit$promoter == "marRABp" &
    dstressr_only_lit$compound == "Acetylsalicylic acid",
  "Strong analog support",
  "Acetylsalicylic acid is salicylate-related; salicylate is the canonical marRABp cue through MarR de-repression.",
  "Alekshun and Levy review; https://pmc.ncbi.nlm.nih.gov/articles/PMC4031692/"
)
set_support(
  dstressr_only_lit$promoter == "robp" &
    dstressr_only_lit$compound == "Paraquat",
  "Direct regulatory support",
  "Paraquat/SoxRS-dependent down-regulation of rob transcription has been reported; DStressR detects robp repression.",
  "Michan et al. 2002 / Duval and Lister 2013; https://pmc.ncbi.nlm.nih.gov/articles/PMC3430332/"
)
set_support(
  dstressr_only_lit$promoter == "soxSp" &
    dstressr_only_lit$compound == "Quercetin",
  "Mechanistic oxidative-stress support",
  "Quercetin can generate oxidative/redox stress in E. coli; SoxRS is the canonical redox-stress system.",
  "Kwun and Lee 2024; https://pmc.ncbi.nlm.nih.gov/articles/PMC11294654/"
)
set_support(
  dstressr_only_lit$promoter == "acrABp" &
    dstressr_only_lit$compound %in% c(
      "Chloramphenicol", "Fusidic acid", "Levofloxacin", "Moxifloxacin",
      "Spiramycin", "Nitrofurantoin", "Procaine"
    ),
  "Same pathway/substrate support",
  "AcrAB-TolC is a broad multidrug efflux system regulated by MarA/SoxS/Rob; several of these compounds or their drug classes are established AcrAB-TolC substrates.",
  "Efflux reviews; https://pmc.ncbi.nlm.nih.gov/articles/PMC1471989/ and https://pmc.ncbi.nlm.nih.gov/articles/PMC10924465/"
)
set_support(
  dstressr_only_lit$compound %in% c(
    "Chloramphenicol", "Fusidic acid", "Levofloxacin", "Moxifloxacin",
    "Novobiocin", "Rifampicin", "Trimethoprim", "Ethidiumbromide",
    "Tetracycline", "Doxycycline"
  ) &
    dstressr_only_lit$promoter %in% c("marRABp", "micFp", "ompFp", "robp", "soxSp"),
  "Regulon/substrate-class support",
  "The compound is an AcrAB-TolC substrate or belongs to a class coupled to the mar-sox-rob permeability network; promoter-specific direction still needs validation.",
  "Mar-Sox-Rob and efflux reviews; https://pmc.ncbi.nlm.nih.gov/articles/PMC4031692/ and https://pmc.ncbi.nlm.nih.gov/articles/PMC1471989/"
)
set_support(
  dstressr_only_lit$promoter == "marRABp" &
    dstressr_only_lit$compound %in% c("Gentamicin", "Meropenem", "Procaine", "Puromycin", "Quercetin"),
  "Growth/stress-response support",
  "Binsfeld et al. report that marRABp activity is strongly growth/stress-correlated, supporting some additional marRABp calls as biologically plausible but not exact prior interactions.",
  "Binsfeld et al. 2025; https://doi.org/10.1371/journal.pbio.3003260"
)
set_support(
  dstressr_only_lit$promoter == "soxSp" &
    dstressr_only_lit$compound %in% c("2-2-Bipyridyl", "Phleomycin", "Tobramycin", "Tunicamycin"),
  "Stress-regulon support",
  "SoxRS responds to redox-active and stress-generating cues; this is plausible network support rather than exact pair-level evidence.",
  "Mar-Sox-Rob review; https://pmc.ncbi.nlm.nih.gov/articles/PMC4031692/"
)
set_support(
  dstressr_only_lit$promoter == "marRABp" &
    dstressr_only_lit$compound %in% c("Chloramphenicol", "Vanillin"),
  "Direct same-pair support",
  "Binsfeld et al. explicitly note these as known marRABp interactions missed by their stringent screen; DStressR recovers them.",
  "Binsfeld et al. 2025; https://doi.org/10.1371/journal.pbio.3003260"
)
set_support(
  dstressr_only_lit$promoter == "robp" &
    dstressr_only_lit$compound == "Paraquat",
  "Direct regulatory support",
  "Paraquat/SoxRS-dependent down-regulation of rob transcription has been reported; DStressR detects robp repression.",
  "Michan et al. 2002 / Duval and Lister 2013; https://pmc.ncbi.nlm.nih.gov/articles/PMC3430332/"
)

dstressr_only_lit <- dstressr_only_lit[order(
  factor(
    dstressr_only_lit$literature_support,
    levels = c(
      "Direct same-pair support",
      "Direct regulatory support",
      "Strong analog support",
      "Same pathway/substrate support",
      "Regulon/substrate-class support",
      "Mechanistic oxidative-stress support",
      "Growth/stress-response support",
      "Stress-regulon support",
      "Hypothesis-generating; no direct prior support identified in this quick review"
    )
  ),
  dstressr_only_lit$promoter,
  dstressr_only_lit$compound
), ]
utils::write.table(
  dstressr_only_lit,
  file.path(out_dir, "dstressr_only_literature_support.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

sections <- lapply(class_levels, function(class_name) {
  label <- paste0("<h3>", html_escape(class_name), "</h3>")
  rows <- appendix[appendix$three_method_class == class_name, ]
  paste0(
    label,
    html_table(
      rows,
      columns = c(
        "promoter", "compound", "mean_z", "binsfeld_padj",
        "modeled_effect", "modeled_padj_by_promoter",
        "standard_effect", "standard_padj_by_promoter",
        "modeled_hit_class", "standard_hit_class"
      )
    )
  )
})

report_file <- file.path(out_dir, "binsfeld_three_method_shareable_report.html")
created <- format(Sys.time(), "%Y-%m-%d %H:%M %Z")

html <- paste0(
'<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Three-method DStressR comparison on the Binsfeld et al. reporter screen</title>
<style>
body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; line-height: 1.48; color: #111827; margin: 0; background: #f8fafc; }
main { max-width: 1120px; margin: 0 auto; padding: 34px 28px 56px; background: white; }
h1 { font-size: 32px; margin: 0 0 8px; }
h2 { margin-top: 34px; border-top: 1px solid #e5e7eb; padding-top: 24px; }
h3 { margin-top: 24px; }
p, li { font-size: 15px; }
.meta { color: #4b5563; margin-bottom: 26px; }
.callout { background: #eff6ff; border-left: 4px solid #2563eb; padding: 12px 16px; margin: 18px 0; }
.grid { display: grid; grid-template-columns: repeat(2, minmax(0, 1fr)); gap: 18px; }
.figure { margin: 22px 0; }
.figure img { width: 100%; border: 1px solid #e5e7eb; background: white; }
.caption { color: #4b5563; font-size: 13px; margin-top: 6px; }
table { border-collapse: collapse; width: 100%; margin: 12px 0 22px; font-size: 13px; }
th, td { border: 1px solid #e5e7eb; padding: 6px 8px; vertical-align: top; }
th { background: #f3f4f6; text-align: left; }
code { background: #f3f4f6; padding: 1px 4px; border-radius: 3px; }
.small { font-size: 13px; color: #4b5563; }
@media print { body { background: white; } main { max-width: none; padding: 20px; } .grid { grid-template-columns: 1fr; } }
</style>
</head>
<body><main>
<h1>Three-method comparison on the Binsfeld et al. reporter screen</h1>
<div class="meta">Self-contained local report generated ', html_escape(created), ' from the DStressR repository.</div>

<div class="callout">
<strong>Headline result.</strong> The reproduced Binsfeld-style WT analysis calls ', set_counts[["binsfeld"]], ' hits, DStressR with the default modeled response calls ', set_counts[["modeled"]], ' hits, and DStressR with fixed alpha=1 response calls ', set_counts[["standard"]], ' hits. The three-way intersection contains ', region_counts[["all_three"]], ' promoter-compound pairs.
</div>

<h2>Data and Analysis Scope</h2>
<p>This report uses the public <em>E. coli</em> reporter-screen data from Binsfeld et al. The comparison is limited to WT rows. The author-style analysis uses WT Z-scores, while both DStressR analyses use AUC rows with <code>removed == "No"</code>, water controls, technical terms for replicate and dose level, empirical-Bayes moderation, promoter-wise FDR adjustment, and the EVC reporter excluded from the tested promoter set.</p>

<h2>Methods Compared</h2>
<ul>
<li><strong>Binsfeld-style rule:</strong> Wilcoxon tests comparing WT Z-score replicates/concentrations against water controls, promoter-wise BH adjustment, and hits at adjusted p-value &lt; 0.05 with absolute mean Z-score &gt; 1.</li>
<li><strong>DStressR modeled response:</strong> <code>log2(lux_auc) - alpha_g * log2(od_auc)</code>, with promoter-specific <code>alpha_g</code> estimated from water controls.</li>
<li><strong>DStressR alpha=1 response:</strong> fixed <code>log2(lux_auc) - log2(od_auc)</code>, otherwise using the same DStressR model settings as the modeled-response run.</li>
</ul>

<h3>Estimated Growth-Response Exponents</h3>',
html_table(
  growth_table,
  columns = intersect(c("promoter", "alpha_raw", "alpha_raw_se", "alpha_shrunk", "alpha_covariates"), names(growth_table))
),

'<h2>Numerical Summary</h2>',
html_table(count_summary),
'<h3>Pairwise Overlaps Within the Three-method Run</h3>',
html_table(pairwise_summary),

'<h2>Literature Review of DStressR-only Hits</h2>
<p>The table below summarizes DStressR-only calls against prior biological knowledge. The review is intentionally conservative: exact same-pair support is separated from broader support through shared drug classes, AcrAB-TolC substrates, or the mar-sox-rob stress-regulatory network.</p>',
html_table(
  dstressr_only_lit,
  columns = c(
    "promoter", "compound", "three_method_class",
    "literature_support", "support_note", "source"
  ),
  max_rows = 40
),
'<p class="small">The complete table is written as <code>analysis/outputs/binsfeld_three_method/dstressr_only_literature_support.tsv</code>.</p>',

'<h2>Figures</h2>
<div class="figure"><img alt="Three-method hit overlap Venn diagram" src="', img_uri(file.path(out_dir, "three_method_hit_overlap_venn.png")), '"><div class="caption">Figure 1. Exact promoter-compound hit overlap across the three hit sets.</div></div>
<div class="figure"><img alt="Three-method p-value histograms" src="', img_uri(file.path(out_dir, "three_method_pvalue_histograms.png")), '"><div class="caption">Figure 2. Raw and promoter-wise BH-adjusted p-value distributions. Colored bars are method-specific hits.</div></div>
<div class="figure"><img alt="Three-method effect histograms" src="', img_uri(file.path(out_dir, "three_method_effect_histograms.png")), '"><div class="caption">Figure 3. Effect-score distributions for the Binsfeld-style mean Z-score and the two DStressR specific-effect estimates.</div></div>
<div class="figure"><img alt="Binsfeld and modeled-response DStressR volcano plots" src="', img_uri(file.path(out_dir, "binsfeld_modeled_volcano_plots.png")), '"><div class="caption">Figure 4. Volcano plots for the Binsfeld-style analysis and the default modeled-response DStressR analysis.</div></div>
<div class="grid">
<div class="figure"><img alt="Binsfeld-referenced DStressR effect scatter" src="', img_uri(file.path(out_dir, "binsfeld_reference_effect_scatter.png")), '"><div class="caption">Figure 5. Binsfeld-style mean Z-score as reference, with DStressR modeled-response and alpha=1 specific effects shown in separate panels.</div></div>
<div class="figure"><img alt="Binsfeld-referenced DStressR p-value scatter" src="', img_uri(file.path(out_dir, "binsfeld_reference_pvalue_scatter.png")), '"><div class="caption">Figure 6. Binsfeld-style raw Wilcoxon p-value as reference, with DStressR modeled-response and alpha=1 raw p-values shown in separate panels.</div></div>
</div>

<h2>Appendix: Union of Significant Pairings</h2>
<p>The full union table is written as <code>analysis/outputs/binsfeld_three_method/binsfeld_three_method_significant_union.tsv</code>. The pairwise report files under <code>analysis/outputs/binsfeld/</code> are not regenerated or modified by this three-method report script.</p>',
paste(sections, collapse = "\n"),

'<h3>Full Three-method Union</h3>',
html_table(
  appendix,
  columns = c(
    "promoter", "compound", "three_method_class",
    "mean_z", "binsfeld_padj", "binsfeld_direction",
    "modeled_effect", "modeled_padj_by_promoter",
    "standard_effect", "standard_padj_by_promoter",
    "modeled_hit_class", "standard_hit_class"
  )
),

'<h2>Reproducibility</h2>
<p>Regenerate this report from the repository root with:</p>
<pre><code>Rscript analysis/ecoli_promoter_screen/build_three_method_report.R</code></pre>
<p class="small">Primary public sources: PLOS Biology article DOI 10.1371/journal.pbio.3003260 and Zenodo DOI 10.5281/zenodo.15600688.</p>
</main></body></html>'
)

writeLines(html, report_file, useBytes = TRUE)

message("Wrote three-method report: ", report_file)
print(count_summary)
