source(file.path("analysis", "_helpers.R"))
source(file.path("analysis", "_paper_artifacts.R"))
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
  reporter = "promoter",
  perturbation = "compound",
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
  reporter = "pseudo_reporter",
  perturbation = "compound",
  lux = "gfp_auc",
  growth = "od_auc",
  controls = "Calibration",
  min_control_n = 3,
  shrink = TRUE
)
alpha_window_parts <- do.call(rbind, strsplit(alpha_window_fit$reporter, " \\| "))
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

base_alpha <- stats::setNames(alpha_fit$alpha_shrunk, alpha_fit$reporter)
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
  reporter = "pseudo_reporter",
  perturbation = "compound",
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
  adjustment = "by_reporter",
  interaction = FALSE
)

saveRDS(fit, file.path(out_dir, "dryad_weak_stress_windows_calibrated_alpha_fit.rds"))
res <- results(fit)
res$promoter <- res$reporter
res$compound <- res$perturbation
parts <- do.call(rbind, strsplit(as.character(res$promoter), " \\| "))
res$base_promoter <- parts[, 1]
res$window <- parts[, 2]
res$base_promoter <- factor(res$base_promoter, levels = promoter_levels)
res$window <- factor(res$window, levels = window_definitions$window)
hit_table <- call_hits(
  res,
  fdr = 0.05,
  effect = "specific_effect",
  padj = "specific_padj_by_reporter"
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

comparison <- merge(
  utils::read.delim(file.path(out_dir, "dryad_weak_stress_windows_alpha1_pair_results.tsv"))[
    , c("promoter", "compound", "total_effect", "specific_effect", "specific_pvalue", "specific_padj_by_reporter", "hit")
  ],
  res[, c("promoter", "compound", "total_effect", "specific_effect", "specific_pvalue", "specific_padj_by_reporter", "hit")],
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
write_paper_artifact(file.path(out_dir, "dash_calibration_paper.rds"),
  results = list(default = res), parameters = list(alpha_fit = alpha_fit, alpha_window_fit = alpha_window_fit),
  annotations = list(promoter_levels = promoter_levels, window_definitions = window_definitions,
                     pseudo_reporter_levels = pseudo_reporter_levels))
