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

fluorescence_path <- file.path(
  raw_data_dir,
  "Spectrophotometry", "Fluorescence", "Weak_Stresses", "Weak stresses.xlsx"
)
od_path <- file.path(
  raw_data_dir,
  "Spectrophotometry", "OD", "Stress_growth_curves.xlsx"
)

if (!file.exists(fluorescence_path)) {
  stop("Missing fluorescence workbook: ", fluorescence_path, call. = FALSE)
}
if (!file.exists(od_path)) {
  stop("Missing OD workbook: ", od_path, call. = FALSE)
}

clean_condition <- function(x) {
  z <- tolower(gsub("[^A-Za-z0-9]+", "", x))
  if (z %in% c("standard", "mg1655", "nostress", "control")) return("Standard")
  if (grepl("tetra|^tet$", z)) return("Tetracycline")
  if (grepl("kana|kan", z)) return("Kanamycin")
  if (grepl("h2o2|peroxide|oxd|oxid", z)) return("H2O2")
  if (grepl("iron|fe", z)) return("Iron")
  x
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

read_raw_sheet <- function(path, sheet) {
  as.data.frame(
    readxl::read_excel(path, sheet = sheet, col_names = FALSE, .name_repair = "minimal"),
    stringsAsFactors = FALSE
  )
}

extract_block <- function(raw, label_col, time_col, sample_cols, promoter, compound, sheet, background = NULL) {
  time <- suppressWarnings(as.numeric(unlist(raw[-c(1, 2), time_col], use.names = FALSE)))
  do.call(rbind, lapply(seq_along(sample_cols), function(k) {
    value <- suppressWarnings(as.numeric(unlist(raw[-c(1, 2), sample_cols[k]], use.names = FALSE)))
    if (!is.null(background)) {
      value <- value - stats::approx(
        x = background$time,
        y = background$value,
        xout = time,
        rule = 2,
        ties = mean
      )$y
    }
    data.frame(
      promoter = promoter,
      compound = compound,
      replicate = paste0("sample_", k),
      sheet = sheet,
      time = time,
      value = value,
      source_label = as.character(raw[1, label_col]),
      stringsAsFactors = FALSE
    )
  }))
}

read_fluorescence <- function(path) {
  sheets <- readxl::excel_sheets(path)
  out <- list()
  diagnostics <- data.frame()

  for (sheet in sheets) {
    raw <- read_raw_sheet(path, sheet)
    labels <- as.character(unlist(raw[1, ], use.names = FALSE))
    label_cols <- which(!is.na(labels) & grepl("^(Fur|MarA|SoxS|LexA)(\\+.*)?$", labels))
    background_time <- suppressWarnings(as.numeric(unlist(raw[-c(1, 2), 2], use.names = FALSE)))
    background_value <- suppressWarnings(as.numeric(unlist(raw[-c(1, 2), 3], use.names = FALSE)))
    background <- data.frame(time = background_time, value = background_value)
    background <- background[is.finite(background$time) & is.finite(background$value), , drop = FALSE]
    diagnostics <- rbind(
      diagnostics,
      data.frame(
        workbook = basename(path),
        sheet = sheet,
        labels = paste(labels[label_cols], collapse = "; "),
        stringsAsFactors = FALSE
      )
    )

    for (label_col in label_cols) {
      label <- labels[label_col]
      promoter <- sub("\\+.*$", "", label)
      compound <- if (grepl("\\+", label)) clean_condition(sub("^.*\\+", "", label)) else "Standard"
      out[[length(out) + 1]] <- extract_block(
        raw = raw,
        label_col = label_col,
        time_col = label_col - 1,
        sample_cols = label_col + 0:2,
        promoter = promoter,
        compound = compound,
        sheet = sheet,
        background = background
      )
    }
  }

  list(data = do.call(rbind, out), diagnostics = diagnostics)
}

read_od <- function(path) {
  raw <- read_raw_sheet(path, "Weak Stress")
  labels <- as.character(unlist(raw[1, ], use.names = FALSE))
  label_cols <- which(!is.na(labels) & grepl("^MG1655", labels))
  out <- list()

  for (label_col in label_cols) {
    label <- labels[label_col]
    compound <- if (grepl("\\+", label)) clean_condition(sub("^.*\\+", "", label)) else "Standard"
    time_col <- label_col + 1
    sample_cols <- label_col + 2:4
    time <- suppressWarnings(as.numeric(unlist(raw[-1, time_col], use.names = FALSE)))
    out[[length(out) + 1]] <- do.call(rbind, lapply(seq_along(sample_cols), function(k) {
      value <- suppressWarnings(as.numeric(unlist(raw[-1, sample_cols[k]], use.names = FALSE)))
      data.frame(
        compound = compound,
        replicate = paste0("sample_", k),
        time = time,
        od = value,
        source_label = label,
        stringsAsFactors = FALSE
      )
    }))
  }

  data <- do.call(rbind, out)
  diagnostics <- data.frame(
    workbook = basename(path),
    sheet = "Weak Stress",
    labels = paste(labels[label_cols], collapse = "; "),
    stringsAsFactors = FALSE
  )
  list(data = data, diagnostics = diagnostics)
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

summarize_windows <- function(d, by, value_col, out_col) {
  rows <- list()
  idx <- 1
  groups <- unique(d[by])
  for (i in seq_len(nrow(groups))) {
    g <- groups[i, , drop = FALSE]
    keep <- rep(TRUE, nrow(d))
    for (nm in by) {
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

average_standard_controls <- function(d) {
  standard <- d[d$compound == "Standard", , drop = FALSE]
  stressed <- d[d$compound != "Standard", , drop = FALSE]
  if (!nrow(standard)) return(d)

  by <- c("promoter", "compound", "replicate", "window", "start_min", "end_min")
  keys <- interaction(standard[by], drop = TRUE)
  averaged <- do.call(rbind, lapply(split(standard, keys), function(x) {
    ans <- as.list(x[1, by, drop = FALSE])
    ans$sheet <- "averaged_standard"
    ans$gfp_auc <- mean(x$gfp_auc, na.rm = TRUE)
    ans$n_time <- min(x$n_time, na.rm = TRUE)
    ans$n_standard_sheets <- length(unique(x$sheet))
    as.data.frame(ans, stringsAsFactors = FALSE)
  }))
  stressed$n_standard_sheets <- NA_integer_
  row.names(averaged) <- NULL
  row.names(stressed) <- NULL
  rbind(averaged, stressed[, names(averaged), drop = FALSE])
}

fluorescence <- read_fluorescence(fluorescence_path)
od <- read_od(od_path)

utils::write.table(
  fluorescence$diagnostics,
  file.path(out_dir, "fluorescence_workbook_diagnostics.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
utils::write.table(
  od$diagnostics,
  file.path(out_dir, "od_workbook_diagnostics.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

fluorescence_auc <- summarize_windows(
  fluorescence$data,
  by = c("promoter", "compound", "replicate", "sheet"),
  value_col = "value",
  out_col = "gfp_auc"
)
fluorescence_auc <- average_standard_controls(fluorescence_auc)
od_auc <- summarize_windows(
  od$data,
  by = c("compound", "replicate"),
  value_col = "od",
  out_col = "od_auc"
)

screen <- merge(
  fluorescence_auc,
  od_auc[, c("compound", "replicate", "window", "od_auc")],
  by = c("compound", "replicate", "window"),
  all.x = TRUE,
  sort = FALSE
)
screen <- screen[is.finite(screen$gfp_auc) & is.finite(screen$od_auc), , drop = FALSE]
screen <- screen[screen$gfp_auc > 0 & screen$od_auc > 0, , drop = FALSE]
screen$promoter <- factor(screen$promoter, levels = promoter_levels)
screen$window <- factor(screen$window, levels = window_definitions$window)
screen$pseudo_reporter <- paste(screen$promoter, screen$window, sep = " | ")
screen$pseudo_reporter <- factor(
  screen$pseudo_reporter,
  levels = pseudo_reporter_levels
)
screen$compound <- factor(screen$compound, levels = c("Standard", "Iron", "Kanamycin", "Tetracycline", "H2O2"))

assay <- prepare_assay(
  screen,
  promoter = "pseudo_reporter",
  compound = "compound",
  control = "Standard",
  lux = "gfp_auc",
  growth = "od_auc",
  growth_exponent = 1,
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
  file.path(out_dir, "dryad_weak_stress_windows_alpha1_input.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
utils::write.table(
  res,
  file.path(out_dir, "dryad_weak_stress_windows_alpha1_pair_results.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
utils::write.table(
  res[res$hit, , drop = FALSE],
  file.path(out_dir, "dryad_weak_stress_windows_alpha1_significant_pairs.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

summary <- data.frame(
  metric = c(
    "base_reporters",
    "time_windows",
    "pseudo_reporters",
    "conditions_including_reference",
    "reference_condition",
    "growth_exponent",
    "well_window_summaries",
    "tested_pseudo_reporter_condition_pairs",
    "significant_pairs"
  ),
  value = c(
    length(unique(screen$promoter)),
    length(unique(screen$window)),
    length(unique(screen$pseudo_reporter)),
    length(unique(screen$compound)),
    "Standard",
    "fixed alpha = 1",
    nrow(screen),
    nrow(res),
    sum(res$hit, na.rm = TRUE)
  ),
  stringsAsFactors = FALSE
)
utils::write.table(
  summary,
  file.path(out_dir, "dryad_weak_stress_windows_alpha1_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

p_hist <- ggplot2::ggplot(res, ggplot2::aes(x = specific_pvalue)) +
  ggplot2::geom_histogram(bins = 20, fill = "#009E73", color = "white", linewidth = 0.2) +
  ggplot2::facet_grid(base_promoter ~ window) +
  ggplot2::theme_light(base_size = 9) +
  ggplot2::labs(
    title = "Dryad weak-stress reporter windows: DStressR p-values with alpha = 1",
    x = "Raw p-value",
    y = "Number of reporter-condition pairs"
  )
ggplot2::ggsave(file.path(out_dir, "dryad_weak_stress_windows_alpha1_pvalue_histogram.png"), p_hist, width = 8, height = 6, dpi = 300)
ggplot2::ggsave(file.path(out_dir, "dryad_weak_stress_windows_alpha1_pvalue_histogram.pdf"), p_hist, width = 8, height = 6)

p_hist_combined <- ggplot2::ggplot(res, ggplot2::aes(x = specific_pvalue)) +
  ggplot2::geom_histogram(
    breaks = seq(0, 1, by = 0.05),
    fill = "#009E73",
    color = "white",
    linewidth = 0.25,
    boundary = 0
  ) +
  ggplot2::theme_light(base_size = 10) +
  ggplot2::labs(
    title = "Dryad weak-stress reporter windows: combined p-values",
    subtitle = paste0(nrow(res), " reporter-window/stress tests; ", sum(res$hit), " significant at within-reporter FDR 0.05"),
    x = "Raw p-value",
    y = "Number of tests"
  )
ggplot2::ggsave(file.path(out_dir, "dryad_weak_stress_windows_alpha1_pvalue_histogram_combined.png"), p_hist_combined, width = 7, height = 4.5, dpi = 300)
ggplot2::ggsave(file.path(out_dir, "dryad_weak_stress_windows_alpha1_pvalue_histogram_combined.pdf"), p_hist_combined, width = 7, height = 4.5)

p_volcano <- plot_volcano(
  res,
  effect = "specific_effect",
  padj = "specific_padj_by_promoter",
  title = "Dryad weak-stress reporter windows: DStressR volcano with alpha = 1",
  label_by = "pair",
  top_n = 15,
  top_promoters = 8
)
ggplot2::ggsave(file.path(out_dir, "dryad_weak_stress_windows_alpha1_volcano.png"), p_volcano, width = 8, height = 5.5, dpi = 300)
ggplot2::ggsave(file.path(out_dir, "dryad_weak_stress_windows_alpha1_volcano.pdf"), p_volcano, width = 8, height = 5.5)

plot_window_heatmap <- function(tab, effect_col = "specific_effect", title = NULL, fill_label = NULL) {
  d <- tab
  d$compound <- factor(d$compound, levels = c("Iron", "Tetracycline", "H2O2", "Kanamycin"))
  d$pseudo_reporter <- factor(
    d$promoter,
    levels = rev(pseudo_reporter_levels)
  )
  limit <- max(abs(d[[effect_col]]), na.rm = TRUE)
  if (is.null(title)) {
    title <- paste("Dryad weak-stress reporter windows: DStressR", effect_col, "with alpha = 1")
  }
  if (is.null(fill_label)) {
    fill_label <- effect_col
  }
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
      axis.text.x = ggplot2::element_text(angle = 30, hjust = 1),
      plot.title.position = "plot"
    ) +
    ggplot2::labs(
      title = title,
      x = "Weak-stress condition",
      y = "Reporter-window",
      fill = fill_label
    )
}
p_heatmap <- plot_window_heatmap(
  res,
  effect_col = "specific_effect",
  title = "Dryad weak-stress reporter windows: DStressR specific effects with alpha = 1",
  fill_label = "Specific effect"
)
ggplot2::ggsave(file.path(out_dir, "dryad_weak_stress_windows_alpha1_effect_heatmap.png"), p_heatmap, width = 7.5, height = 6, dpi = 300)
ggplot2::ggsave(file.path(out_dir, "dryad_weak_stress_windows_alpha1_effect_heatmap.pdf"), p_heatmap, width = 7.5, height = 6)

p_total_heatmap <- plot_window_heatmap(
  res,
  effect_col = "total_effect",
  title = "Dryad weak-stress reporter windows: DStressR total effects with alpha = 1",
  fill_label = "Total effect"
)
ggplot2::ggsave(file.path(out_dir, "dryad_weak_stress_windows_alpha1_total_effect_heatmap.png"), p_total_heatmap, width = 7.5, height = 6, dpi = 300)
ggplot2::ggsave(file.path(out_dir, "dryad_weak_stress_windows_alpha1_total_effect_heatmap.pdf"), p_total_heatmap, width = 7.5, height = 6)

message("Wrote Dryad weak-stress windowed DStressR outputs to: ", out_dir)
