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
  time <- suppressWarnings(as.numeric(raw[-1, 2]))
  out <- do.call(rbind, lapply(1:3, function(k) {
    value <- suppressWarnings(as.numeric(raw[-1, 2 + k]))
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
  background <- split(
    control[, c("time", "background_gfp")],
    control$replicate
  )
  background <- lapply(background, function(x) {
    names(x) <- c("time", "value")
    x
  })

  reporter_sheets <- setdiff(sheets, "Control")
  out <- do.call(rbind, lapply(reporter_sheets, function(sheet) {
    read_single_sheet_timeseries(path, sheet, "gfp", subtract_background = background)
  }))
  out
}

read_growth_transition_od <- function(path) {
  raw <- read_raw_sheet(path, "OD600nm")
  labels <- as.character(unlist(raw[1, ], use.names = FALSE))
  label_cols <- seq(1, ncol(raw), by = 7)
  label_cols <- label_cols[label_cols <= ncol(raw)]
  label_cols <- label_cols[!is.na(labels[label_cols])]

  out <- do.call(rbind, lapply(label_cols, function(label_col) {
    reporter <- canonical_reporter(labels[label_col])
    time <- suppressWarnings(as.numeric(raw[-1, label_col + 1]))
    do.call(rbind, lapply(1:3, function(k) {
      value <- suppressWarnings(as.numeric(raw[-1, label_col + 1 + k]))
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
  compound = c("Baseline", "Early transition", "Late transition", "Stationary"),
  start = c(0, 320, 620, 920),
  end = c(300, 600, 900, 1260),
  stringsAsFactors = FALSE
)

summarize_windows <- function(d, value_col, out_col) {
  rows <- list()
  idx <- 1
  groups <- unique(d[, c("promoter", "replicate")])
  for (i in seq_len(nrow(groups))) {
    g <- groups[i, , drop = FALSE]
    gd <- d[d$promoter == g$promoter & d$replicate == g$replicate, , drop = FALSE]
    for (w in seq_len(nrow(window_definitions))) {
      wd <- window_definitions[w, ]
      x <- gd[gd$time >= wd$start & gd$time <= wd$end, , drop = FALSE]
      rows[[idx]] <- data.frame(
        promoter = g$promoter,
        replicate = g$replicate,
        compound = wd$compound,
        start_min = wd$start,
        end_min = wd$end,
        value = auc_trapezoid(x$time, x[[value_col]]),
        n_time = length(unique(x$time[is.finite(x$time)])),
        stringsAsFactors = FALSE
      )
      idx <- idx + 1
    }
  }
  out <- do.call(rbind, rows)
  names(out)[names(out) == "value"] <- out_col
  out
}

gfp <- read_growth_transition_gfp(gfp_path)
od <- read_growth_transition_od(od_path)

gfp_auc <- summarize_windows(gfp, "gfp", "gfp_auc")
od_auc <- summarize_windows(od, "od", "od_auc")

screen <- merge(
  gfp_auc,
  od_auc[, c("promoter", "replicate", "compound", "od_auc")],
  by = c("promoter", "replicate", "compound"),
  all = FALSE,
  sort = FALSE
)
screen <- screen[is.finite(screen$gfp_auc) & is.finite(screen$od_auc), , drop = FALSE]
screen <- screen[screen$gfp_auc > 0 & screen$od_auc > 0, , drop = FALSE]
screen$compound <- factor(screen$compound, levels = window_definitions$compound)

assay <- prepare_assay(
  screen,
  promoter = "promoter",
  compound = "compound",
  control = "Baseline",
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
  file.path(out_dir, "dryad_growth_transition_auc_input.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
utils::write.table(
  res,
  file.path(out_dir, "dryad_growth_transition_alpha1_pair_results.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
utils::write.table(
  res[res$hit, , drop = FALSE],
  file.path(out_dir, "dryad_growth_transition_alpha1_significant_pairs.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

summary <- data.frame(
  metric = c(
    "reporters",
    "growth_phase_windows_including_reference",
    "reference_window",
    "growth_exponent",
    "well_window_summaries",
    "tested_reporter_window_pairs",
    "significant_pairs"
  ),
  value = c(
    length(unique(screen$promoter)),
    length(unique(screen$compound)),
    "Baseline",
    "fixed alpha = 1",
    nrow(screen),
    nrow(res),
    sum(res$hit, na.rm = TRUE)
  ),
  stringsAsFactors = FALSE
)
utils::write.table(
  summary,
  file.path(out_dir, "dryad_growth_transition_alpha1_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

p_hist <- ggplot2::ggplot(res, ggplot2::aes(x = specific_pvalue)) +
  ggplot2::geom_histogram(bins = 30, fill = "#009E73", color = "white", linewidth = 0.2) +
  ggplot2::theme_light(base_size = 10) +
  ggplot2::labs(
    title = "Dryad growth-transition reporters: DStressR p-values with alpha = 1",
    x = "Raw p-value",
    y = "Number of reporter-window pairs"
  )
ggplot2::ggsave(file.path(out_dir, "dryad_growth_transition_alpha1_pvalue_histogram.png"), p_hist, width = 7, height = 4.5, dpi = 300)
ggplot2::ggsave(file.path(out_dir, "dryad_growth_transition_alpha1_pvalue_histogram.pdf"), p_hist, width = 7, height = 4.5)

p_volcano <- plot_volcano(
  res,
  effect = "specific_effect",
  padj = "specific_padj_by_promoter",
  title = "Dryad growth-transition reporters: DStressR volcano with alpha = 1",
  label_by = "pair",
  top_n = 15,
  top_promoters = 8
)
ggplot2::ggsave(file.path(out_dir, "dryad_growth_transition_alpha1_volcano.png"), p_volcano, width = 8, height = 5.5, dpi = 300)
ggplot2::ggsave(file.path(out_dir, "dryad_growth_transition_alpha1_volcano.pdf"), p_volcano, width = 8, height = 5.5)

plot_transition_heatmap <- function(tab) {
  d <- tab
  d$compound <- factor(d$compound, levels = window_definitions$compound[-1])
  reporter_order <- stats::aggregate(abs(d$specific_effect), by = list(promoter = d$promoter), FUN = max)
  reporter_order <- reporter_order[order(reporter_order$x), "promoter"]
  d$promoter <- factor(d$promoter, levels = reporter_order)
  limit <- stats::quantile(abs(d$specific_effect), 0.98, na.rm = TRUE)
  ggplot2::ggplot(d, ggplot2::aes(compound, promoter, fill = pmax(pmin(specific_effect, limit), -limit))) +
    ggplot2::geom_tile(color = "white", linewidth = 0.25) +
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
      axis.text.x = ggplot2::element_text(angle = 25, hjust = 1),
      plot.title.position = "plot"
    ) +
    ggplot2::labs(
      title = "Dryad growth-transition reporters: DStressR effects with alpha = 1",
      x = "Growth-phase window",
      y = "Reporter",
      fill = "Specific effect"
    )
}
p_heatmap <- plot_transition_heatmap(res)
ggplot2::ggsave(file.path(out_dir, "dryad_growth_transition_alpha1_effect_heatmap.png"), p_heatmap, width = 7.5, height = 6, dpi = 300)
ggplot2::ggsave(file.path(out_dir, "dryad_growth_transition_alpha1_effect_heatmap.pdf"), p_heatmap, width = 7.5, height = 6)

message("Wrote Dryad growth-transition DStressR outputs to: ", out_dir)
