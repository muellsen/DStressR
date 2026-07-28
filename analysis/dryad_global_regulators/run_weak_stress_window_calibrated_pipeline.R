source(file.path("analysis", "_helpers.R"))
load_destress_package()

required_packages <- c("readxl", "ggplot2")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Missing required package(s): ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

readxl <- asNamespace("readxl")
ggplot2 <- asNamespace("ggplot2")

raw_data_dir <- analysis_path("analysis", "dryad_global_regulators", "raw", "Data")
out_dir <- analysis_output_dir("dryad_global_regulators")

weak_input_path <- file.path(out_dir, "dryad_weak_stress_windows_alpha1_input.tsv")
if (!file.exists(weak_input_path)) {
  stop(
    "Missing weak-stress window input. Run ",
    "analysis/dryad_global_regulators/run_weak_stress_window_pipeline.R first.",
    call. = FALSE
  )
}

gfp_path <- file.path(
  raw_data_dir,
  "Spectrophotometry", "Fluorescence", "Growth_Transition", "Raw_GFP_data.xlsx"
)
od_path <- file.path(
  raw_data_dir,
  "Spectrophotometry", "OD", "M9_Glucose_Growth_curves.xlsx"
)

if (!file.exists(gfp_path)) {
  stop("Missing growth-transition GFP workbook: ", gfp_path, call. = FALSE)
}
if (!file.exists(od_path)) {
  stop("Missing M9 glucose OD workbook: ", od_path, call. = FALSE)
}

canonical_reporter <- function(x) {
  key <- toupper(gsub("[^A-Za-z0-9]+", "", x))
  aliases <- c(
    "CRP" = "CRP", "FNR" = "Fnr", "HNS" = "HNS", "FIS" = "Fis",
    "NARL" = "NarL", "FUR" = "Fur", "LRP" = "Lrp", "NSRR" = "NsrR",
    "CRA" = "Cra", "FLHD" = "FlhD", "FLHDC" = "FlhD", "NARP" = "NarP",
    "PHOB" = "PhoB", "LEXA" = "LexA", "PHOP" = "PhoP", "MARA" = "MarA",
    "CPXR" = "CpxR", "PDHR" = "PdhR", "SOXS" = "SoxS", "MG1655" = "MG1655"
  )
  unname(ifelse(key %in% names(aliases), aliases[key], x))
}

read_raw_sheet <- function(path, sheet) {
  as.data.frame(
    readxl::read_excel(path, sheet = sheet, col_names = FALSE, .name_repair = "minimal"),
    stringsAsFactors = FALSE
  )
}

read_single_sheet_timeseries <- function(path, sheet, value_name, subtract_background = NULL) {
  raw <- read_raw_sheet(path, sheet)
  time <- suppressWarnings(as.numeric(unlist(raw[-1, 2], use.names = FALSE)))
  out <- do.call(rbind, lapply(1:3, function(k) {
    value <- suppressWarnings(as.numeric(unlist(raw[-1, 2 + k], use.names = FALSE)))
    if (!is.null(subtract_background)) {
      bg <- subtract_background[[paste0("sample_", k)]]
      value <- value - stats::approx(bg$time, bg$value, xout = time, rule = 2, ties = mean)$y
    }
    data.frame(
      promoter = canonical_reporter(sheet),
      replicate = paste0("sample_", k),
      time = time,
      value = value,
      stringsAsFactors = FALSE
    )
  }))
  names(out)[names(out) == "value"] <- value_name
  out[is.finite(out$time) & is.finite(out[[value_name]]), , drop = FALSE]
}

read_growth_transition_gfp <- function(path) {
  sheets <- readxl::excel_sheets(path)
  control <- read_single_sheet_timeseries(path, "Control", "background_gfp")
  background <- split(control[, c("time", "background_gfp")], control$replicate)
  background <- lapply(background, function(x) {
    names(x) <- c("time", "value")
    x
  })
  reporter_sheets <- setdiff(sheets, "Control")
  do.call(rbind, lapply(reporter_sheets, function(sheet) {
    read_single_sheet_timeseries(path, sheet, "gfp", subtract_background = background)
  }))
}

read_growth_transition_od <- function(path) {
  raw <- read_raw_sheet(path, "OD600nm")
  labels <- as.character(unlist(raw[1, ], use.names = FALSE))
  label_cols <- seq(1, ncol(raw), by = 7)
  label_cols <- label_cols[label_cols <= ncol(raw)]
  label_cols <- label_cols[!is.na(labels[label_cols])]

  out <- do.call(rbind, lapply(label_cols, function(label_col) {
    reporter <- canonical_reporter(labels[label_col])
    time <- suppressWarnings(as.numeric(unlist(raw[-1, label_col + 1], use.names = FALSE)))
    do.call(rbind, lapply(1:3, function(k) {
      value <- suppressWarnings(as.numeric(unlist(raw[-1, label_col + 1 + k], use.names = FALSE)))
      data.frame(
        promoter = reporter,
        replicate = paste0("sample_", k),
        time = time,
        od = value,
        stringsAsFactors = FALSE
      )
    }))
  }))
  out[is.finite(out$time) & is.finite(out$od), , drop = FALSE]
}

auc_trapezoid <- function(time, value) {
  ok <- is.finite(time) & is.finite(value)
  time <- time[ok]
  value <- value[ok]
  if (length(time) < 2) return(NA_real_)
  ord <- order(time)
  time <- time[ord]
  value <- value[ord]
  sum(diff(time) * (head(value, -1) + tail(value, -1)) / 2)
}

window_definitions <- data.frame(
  window = c("Early", "Middle", "Late"),
  start_min = c(60, 120, 200),
  end_min = c(100, 180, 240),
  stringsAsFactors = FALSE
)
promoter_levels <- c("Fur", "MarA", "SoxS", "LexA")
pseudo_reporter_levels <- unlist(
  lapply(promoter_levels, function(promoter) paste(promoter, window_definitions$window, sep = " | ")),
  use.names = FALSE
)

summarize_windows <- function(d, group_cols, value_col, out_col) {
  rows <- list()
  idx <- 1
  groups <- unique(d[group_cols])
  for (i in seq_len(nrow(groups))) {
    g <- groups[i, , drop = FALSE]
    keep <- rep(TRUE, nrow(d))
    for (nm in group_cols) {
      keep <- keep & d[[nm]] == g[[nm]]
    }
    gd <- d[keep, , drop = FALSE]
    for (w in seq_len(nrow(window_definitions))) {
      wd <- window_definitions[w, ]
      x <- gd[gd$time >= wd$start_min & gd$time <= wd$end_min, , drop = FALSE]
      ans <- as.list(g)
      ans$window <- wd$window
      ans$start_min <- wd$start_min
      ans$end_min <- wd$end_min
      ans[[out_col]] <- auc_trapezoid(x$time, x[[value_col]])
      ans$n_time <- length(unique(x$time[is.finite(x$time)]))
      rows[[idx]] <- as.data.frame(ans, stringsAsFactors = FALSE)
      idx <- idx + 1
    }
  }
  out <- do.call(rbind, rows)
  row.names(out) <- NULL
  out
}

gfp <- read_growth_transition_gfp(gfp_path)
od <- read_growth_transition_od(od_path)

calib_gfp_auc <- summarize_windows(gfp, c("promoter", "replicate"), "gfp", "gfp_auc")
calib_od_auc <- summarize_windows(od, c("promoter", "replicate"), "od", "od_auc")
calib <- merge(
  calib_gfp_auc,
  calib_od_auc[, c("promoter", "replicate", "window", "od_auc")],
  by = c("promoter", "replicate", "window"),
  all = FALSE,
  sort = FALSE
)
calib <- calib[is.finite(calib$gfp_auc) & is.finite(calib$od_auc), , drop = FALSE]
calib <- calib[calib$gfp_auc > 0 & calib$od_auc > 0, , drop = FALSE]
calib$compound <- "Calibration"
calib$window <- factor(calib$window, levels = window_definitions$window)

alpha_fit <- estimate_growth_exponents(
  calib,
  promoter = "promoter",
  compound = "compound",
  lux = "gfp_auc",
  growth = "od_auc",
  covariates = "window",
  controls = "Calibration",
  min_control_n = 6,
  shrink = TRUE
)

calib$pseudo_reporter <- paste(calib$promoter, calib$window, sep = " | ")
alpha_window_fit <- estimate_growth_exponents(
  calib,
  promoter = "pseudo_reporter",
  compound = "compound",
  lux = "gfp_auc",
  growth = "od_auc",
  controls = "Calibration",
  min_control_n = 3,
  shrink = TRUE
)
alpha_window_parts <- do.call(rbind, strsplit(alpha_window_fit$promoter, " \\| "))
alpha_window_fit$base_promoter <- alpha_window_parts[, 1]
alpha_window_fit$window <- alpha_window_parts[, 2]
alpha_window_fit$window <- factor(alpha_window_fit$window, levels = window_definitions$window)

utils::write.table(
  calib,
  file.path(out_dir, "dryad_growth_transition_matching_windows_calibration_input.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
utils::write.table(
  alpha_fit,
  file.path(out_dir, "dryad_growth_transition_matching_windows_growth_exponents.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
utils::write.table(
  alpha_window_fit,
  file.path(out_dir, "dryad_growth_transition_matching_windows_growth_exponents_by_window.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

alpha_promoter_long <- rbind(
  data.frame(
    promoter = alpha_fit$promoter,
    panel = "Raw growth exponent",
    estimate = alpha_fit$alpha_raw,
    se = alpha_fit$alpha_raw_se,
    control_n = alpha_fit$control_n,
    stringsAsFactors = FALSE
  ),
  data.frame(
    promoter = alpha_fit$promoter,
    panel = "Shrunken growth exponent",
    estimate = alpha_fit$alpha_shrunk,
    se = alpha_fit$alpha_shrunk_se,
    control_n = alpha_fit$control_n,
    stringsAsFactors = FALSE
  )
)
alpha_promoter_long$panel <- factor(
  alpha_promoter_long$panel,
  levels = c("Raw growth exponent", "Shrunken growth exponent")
)
promoter_alpha_order <- alpha_fit[order(alpha_fit$alpha_shrunk), "promoter"]
alpha_promoter_long$promoter <- factor(alpha_promoter_long$promoter, levels = promoter_alpha_order)
alpha_global <- unique(alpha_fit$alpha_global)
alpha_global <- alpha_global[is.finite(alpha_global)][1]

p_promoter_alpha <- ggplot2::ggplot(
  alpha_promoter_long,
  ggplot2::aes(x = estimate, y = promoter)
) +
  ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "#94A3B8", linewidth = 0.35) +
  ggplot2::geom_vline(xintercept = alpha_global, linetype = "dotted", color = "#475569", linewidth = 0.4) +
  ggplot2::geom_errorbar(
    ggplot2::aes(xmin = estimate - 1.96 * se, xmax = estimate + 1.96 * se),
    width = 0,
    color = "#64748B",
    linewidth = 0.45
  ) +
  ggplot2::geom_point(color = "#1D4ED8", fill = "#93C5FD", size = 2.2) +
  ggplot2::facet_wrap(
    ggplot2::vars(panel),
    nrow = 1
  ) +
  ggplot2::theme_light(base_size = 10) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    strip.background = ggplot2::element_rect(fill = "#F1F5F9", color = "#CBD5E1"),
    strip.text = ggplot2::element_text(color = "#111827"),
    legend.position = "none"
  ) +
  ggplot2::labs(
    x = expression("Growth-response exponent estimate " * hat(alpha)[a]),
    y = "Reporter"
  )
ggplot2::ggsave(file.path(out_dir, "dryad_growth_transition_matching_windows_alpha_shrinkage_estimates.png"), p_promoter_alpha, width = 8, height = 5.8, dpi = 300)
ggplot2::ggsave(file.path(out_dir, "dryad_growth_transition_matching_windows_alpha_shrinkage_estimates.pdf"), p_promoter_alpha, width = 8, height = 5.8)

alpha_window_long <- rbind(
  data.frame(
    base_promoter = alpha_window_fit$base_promoter,
    window = alpha_window_fit$window,
    estimate_type = "Raw",
    alpha = alpha_window_fit$alpha_raw,
    se = alpha_window_fit$alpha_raw_se,
    stringsAsFactors = FALSE
  ),
  data.frame(
    base_promoter = alpha_window_fit$base_promoter,
    window = alpha_window_fit$window,
    estimate_type = "Shrunken",
    alpha = alpha_window_fit$alpha_shrunk,
    se = alpha_window_fit$alpha_shrunk_se,
    stringsAsFactors = FALSE
  )
)
alpha_window_long$estimate_type <- factor(alpha_window_long$estimate_type, levels = c("Raw", "Shrunken"))
promoter_order <- stats::aggregate(alpha_shrunk ~ base_promoter, alpha_window_fit, median)
promoter_order <- promoter_order[order(promoter_order$alpha_shrunk), "base_promoter"]
alpha_window_long$base_promoter <- factor(alpha_window_long$base_promoter, levels = promoter_order)
alpha_limit <- max(abs(alpha_window_long$alpha), na.rm = TRUE)

p_alpha_heatmap <- ggplot2::ggplot(
  alpha_window_long,
  ggplot2::aes(window, base_promoter, fill = pmax(pmin(alpha, alpha_limit), -alpha_limit))
) +
  ggplot2::geom_tile(color = "white", linewidth = 0.3) +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", alpha)), size = 2.6) +
  ggplot2::facet_wrap(ggplot2::vars(estimate_type), nrow = 1) +
  ggplot2::scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-alpha_limit, alpha_limit)
  ) +
  ggplot2::theme_light(base_size = 9) +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    strip.background = ggplot2::element_rect(fill = "grey95", color = "grey75"),
    strip.text = ggplot2::element_text(face = "bold", size = 9, color = "grey15")
  ) +
  ggplot2::labs(
    x = "Time window",
    y = "Reporter",
    fill = expression(hat(alpha)[a])
  )
ggplot2::ggsave(file.path(out_dir, "dryad_growth_transition_matching_windows_alpha_heatmap.png"), p_alpha_heatmap, width = 7.5, height = 7, dpi = 300)
ggplot2::ggsave(file.path(out_dir, "dryad_growth_transition_matching_windows_alpha_heatmap.pdf"), p_alpha_heatmap, width = 7.5, height = 7)

alpha_curve <- alpha_window_fit
alpha_curve$base_promoter <- factor(alpha_curve$base_promoter, levels = promoter_order)
p_alpha_curves <- ggplot2::ggplot(alpha_curve, ggplot2::aes(window, alpha_raw, group = base_promoter)) +
  ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "grey55", linewidth = 0.25) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = alpha_raw - alpha_raw_se, ymax = alpha_raw + alpha_raw_se),
    width = 0.12,
    color = "grey65",
    linewidth = 0.25
  ) +
  ggplot2::geom_point(color = "grey30", size = 1.4) +
  ggplot2::geom_line(ggplot2::aes(y = alpha_shrunk), color = "#009E73", linewidth = 0.55) +
  ggplot2::geom_point(ggplot2::aes(y = alpha_shrunk), color = "#009E73", size = 1.4) +
  ggplot2::facet_wrap(ggplot2::vars(base_promoter), ncol = 3) +
  ggplot2::coord_cartesian(ylim = stats::quantile(c(alpha_curve$alpha_raw, alpha_curve$alpha_shrunk), c(0.02, 0.98), na.rm = TRUE)) +
  ggplot2::theme_light(base_size = 8) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(angle = 30, hjust = 1),
    strip.background = ggplot2::element_rect(fill = "grey95", color = "grey75"),
    strip.text = ggplot2::element_text(face = "bold", size = 8, color = "grey15")
  ) +
  ggplot2::labs(
    x = "Time window",
    y = expression("Growth-response exponent estimate " * hat(alpha)[a])
  )
ggplot2::ggsave(file.path(out_dir, "dryad_growth_transition_matching_windows_alpha_curves.png"), p_alpha_curves, width = 8, height = 8, dpi = 300)
ggplot2::ggsave(file.path(out_dir, "dryad_growth_transition_matching_windows_alpha_curves.pdf"), p_alpha_curves, width = 8, height = 8)

base_alpha <- stats::setNames(alpha_fit$alpha_shrunk, alpha_fit$promoter)
stress_alpha <- unlist(
  lapply(promoter_levels, function(promoter) {
    stats::setNames(rep(base_alpha[[promoter]], nrow(window_definitions)), paste(promoter, window_definitions$window, sep = " | "))
  }),
  use.names = TRUE
)
if (any(!is.finite(stress_alpha))) {
  stop("Could not assign calibrated alpha values to all weak-stress pseudo-reporters.", call. = FALSE)
}

screen <- utils::read.delim(weak_input_path, check.names = FALSE)
screen$promoter <- factor(screen$promoter, levels = promoter_levels)
screen$window <- factor(screen$window, levels = window_definitions$window)
screen$pseudo_reporter <- factor(screen$pseudo_reporter, levels = pseudo_reporter_levels)
screen$compound <- factor(screen$compound, levels = c("Standard", "Iron", "Tetracycline", "H2O2", "Kanamycin"))

assay <- prepare_assay(
  screen,
  promoter = "pseudo_reporter",
  compound = "compound",
  control = "Standard",
  lux = "gfp_auc",
  growth = "od_auc",
  growth_exponent = stress_alpha,
  replicate = "replicate"
)
fit <- fit_destress(
  assay,
  technical = "replicate",
  empirical_bayes = TRUE,
  adjustment = "by_promoter",
  interaction = FALSE
)

res <- results(fit)
parts <- do.call(rbind, strsplit(as.character(res$promoter), " \\| "))
res$base_promoter <- parts[, 1]
res$window <- parts[, 2]
res$base_promoter <- factor(res$base_promoter, levels = promoter_levels)
res$window <- factor(res$window, levels = window_definitions$window)
hit_table <- call_hits(
  res,
  fdr = 0.05,
  effect = "specific_effect",
  padj = "specific_padj_by_promoter"
)
res$hit_class <- hit_table$hit
res$hit <- res$hit_class != "Not DE"

utils::write.table(
  screen,
  file.path(out_dir, "dryad_weak_stress_windows_calibrated_alpha_input.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
utils::write.table(
  data.frame(pseudo_reporter = names(stress_alpha), alpha_calibrated = as.numeric(stress_alpha), stringsAsFactors = FALSE),
  file.path(out_dir, "dryad_weak_stress_windows_calibrated_alpha_supplied.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
utils::write.table(
  res,
  file.path(out_dir, "dryad_weak_stress_windows_calibrated_alpha_pair_results.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
utils::write.table(
  res[res$hit, , drop = FALSE],
  file.path(out_dir, "dryad_weak_stress_windows_calibrated_alpha_significant_pairs.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

summary <- data.frame(
  metric = c(
    "calibration_reporters",
    "calibration_rows",
    "base_stress_reporters",
    "time_windows",
    "pseudo_reporters",
    "conditions_including_reference",
    "reference_condition",
    "growth_exponent",
    "weak_stress_rows",
    "tested_pseudo_reporter_condition_pairs",
    "significant_pairs"
  ),
  value = c(
    length(unique(calib$promoter)),
    nrow(calib),
    length(unique(screen$promoter)),
    length(unique(screen$window)),
    length(unique(screen$pseudo_reporter)),
    length(unique(screen$compound)),
    "Standard",
    "calibrated from no-stress growth-transition reporter panel",
    nrow(screen),
    nrow(res),
    sum(res$hit, na.rm = TRUE)
  ),
  stringsAsFactors = FALSE
)
utils::write.table(
  summary,
  file.path(out_dir, "dryad_weak_stress_windows_calibrated_alpha_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

plot_window_heatmap <- function(tab, effect_col, fill_label) {
  d <- tab
  d$compound <- factor(d$compound, levels = c("Iron", "Tetracycline", "H2O2", "Kanamycin"))
  d$pseudo_reporter <- factor(d$promoter, levels = rev(pseudo_reporter_levels))
  limit <- max(abs(d[[effect_col]]), na.rm = TRUE)
  ggplot2::ggplot(d, ggplot2::aes(compound, pseudo_reporter, fill = .data[[effect_col]])) +
    ggplot2::geom_tile(color = "white", linewidth = 0.35) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", .data[[effect_col]])), size = 2.4) +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-limit, limit)
    ) +
    ggplot2::theme_light(base_size = 9) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 30, hjust = 1)
    ) +
    ggplot2::labs(x = "Weak-stress condition", y = "Reporter-window unit", fill = fill_label)
}

p_total <- plot_window_heatmap(
  res,
  "total_effect",
  "Estimated total effect"
)
ggplot2::ggsave(file.path(out_dir, "dryad_weak_stress_windows_calibrated_alpha_total_effect_heatmap.png"), p_total, width = 7.5, height = 6, dpi = 300)
ggplot2::ggsave(file.path(out_dir, "dryad_weak_stress_windows_calibrated_alpha_total_effect_heatmap.pdf"), p_total, width = 7.5, height = 6)

p_specific <- plot_window_heatmap(
  res,
  "specific_effect",
  "Estimated specific effect"
)
ggplot2::ggsave(file.path(out_dir, "dryad_weak_stress_windows_calibrated_alpha_specific_effect_heatmap.png"), p_specific, width = 7.5, height = 6, dpi = 300)
ggplot2::ggsave(file.path(out_dir, "dryad_weak_stress_windows_calibrated_alpha_specific_effect_heatmap.pdf"), p_specific, width = 7.5, height = 6)

p_hist <- ggplot2::ggplot(res, ggplot2::aes(x = specific_pvalue)) +
  ggplot2::geom_histogram(breaks = seq(0, 1, by = 0.05), fill = "#009E73", color = "white", linewidth = 0.25) +
  ggplot2::theme_light(base_size = 10) +
  ggplot2::labs(
    x = "Raw p-value",
    y = "Number of tests"
  )
ggplot2::ggsave(file.path(out_dir, "dryad_weak_stress_windows_calibrated_alpha_pvalue_histogram_combined.png"), p_hist, width = 7, height = 4.5, dpi = 300)
ggplot2::ggsave(file.path(out_dir, "dryad_weak_stress_windows_calibrated_alpha_pvalue_histogram_combined.pdf"), p_hist, width = 7, height = 4.5)

p_volcano <- plot_volcano(
  res,
  effect = "specific_effect",
  padj = "specific_padj_by_promoter",
  title = NULL,
  label_by = "pair",
  top_n = 15,
  top_promoters = 8,
  xlab = "Estimated specific effect",
  ylab = "-log10 adjusted p-value"
)
ggplot2::ggsave(file.path(out_dir, "dryad_weak_stress_windows_calibrated_alpha_volcano.png"), p_volcano, width = 8, height = 5.5, dpi = 300)
ggplot2::ggsave(file.path(out_dir, "dryad_weak_stress_windows_calibrated_alpha_volcano.pdf"), p_volcano, width = 8, height = 5.5)

comparison <- merge(
  utils::read.delim(file.path(out_dir, "dryad_weak_stress_windows_alpha1_pair_results.tsv"))[
    , c("promoter", "compound", "total_effect", "specific_effect", "specific_pvalue", "specific_padj_by_promoter", "hit")
  ],
  res[, c("promoter", "compound", "total_effect", "specific_effect", "specific_pvalue", "specific_padj_by_promoter", "hit")],
  by = c("promoter", "compound"),
  suffixes = c("_alpha1", "_calibrated_alpha"),
  all = TRUE,
  sort = FALSE
)
utils::write.table(
  comparison,
  file.path(out_dir, "dryad_weak_stress_windows_alpha1_vs_calibrated_alpha_comparison.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Wrote Dryad calibrated-alpha weak-stress outputs to: ", out_dir)
