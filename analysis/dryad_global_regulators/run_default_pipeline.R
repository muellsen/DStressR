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
  time <- suppressWarnings(as.numeric(raw[-c(1, 2), time_col]))
  do.call(rbind, lapply(seq_along(sample_cols), function(k) {
    value <- suppressWarnings(as.numeric(raw[-c(1, 2), sample_cols[k]]))
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
  promoters <- c("Fur", "MarA", "SoxS", "LexA")
  out <- list()
  diagnostics <- data.frame()

  for (sheet in sheets) {
    raw <- read_raw_sheet(path, sheet)
    labels <- as.character(unlist(raw[1, ], use.names = FALSE))
    label_cols <- which(!is.na(labels) & labels %in% c(promoters, paste0(rep(promoters, each = 1), "+", sheet)))
    label_cols <- which(!is.na(labels) & grepl("^(Fur|MarA|SoxS|LexA)(\\+.*)?$", labels))
    background_time <- suppressWarnings(as.numeric(raw[-c(1, 2), 2]))
    background_value <- suppressWarnings(as.numeric(raw[-c(1, 2), 3]))
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
    time <- suppressWarnings(as.numeric(raw[-1, time_col]))
    out[[length(out) + 1]] <- do.call(rbind, lapply(seq_along(sample_cols), function(k) {
      value <- suppressWarnings(as.numeric(raw[-1, sample_cols[k]]))
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

summarize_auc <- function(d, by, value_col, out_col) {
  key <- interaction(d[by], drop = TRUE)
  pieces <- split(d, key)
  out <- do.call(rbind, lapply(pieces, function(x) {
    ans <- as.list(x[1, by, drop = FALSE])
    ans[[out_col]] <- auc_trapezoid(x$time, x[[value_col]])
    ans$n_time <- length(unique(x$time[is.finite(x$time)]))
    as.data.frame(ans, stringsAsFactors = FALSE)
  }))
  row.names(out) <- NULL
  out
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

fluorescence_auc <- summarize_auc(
  fluorescence$data,
  by = c("promoter", "compound", "replicate", "sheet"),
  value_col = "value",
  out_col = "gfp_auc"
)
od_auc <- summarize_auc(
  od$data,
  by = c("compound", "replicate"),
  value_col = "od",
  out_col = "od_auc"
)

screen <- merge(
  fluorescence_auc,
  od_auc[, c("compound", "replicate", "od_auc")],
  by = c("compound", "replicate"),
  all.x = TRUE,
  sort = FALSE
)
screen <- screen[is.finite(screen$gfp_auc) & is.finite(screen$od_auc), , drop = FALSE]
screen$promoter <- factor(screen$promoter, levels = c("Fur", "MarA", "SoxS", "LexA"))
screen$compound <- factor(screen$compound, levels = c("Standard", "Iron", "Kanamycin", "Tetracycline", "H2O2"))

background_labels <- c("EVC", "empty_vector", "emptyvector", "no_promoter")
background_promoter <- intersect(background_labels, unique(as.character(screen$promoter)))[1]
if (is.na(background_promoter)) {
  background_promoter <- NULL
}

assay_args <- list(
  data = screen,
  promoter = "promoter",
  compound = "compound",
  control = "Standard",
  lux = "gfp_auc",
  growth = "od_auc",
  growth_exponent = "estimate",
  replicate = "replicate",
  batch = "sheet",
  growth_covariates = "replicate"
)
if (!is.null(background_promoter)) {
  assay_args$background_promoter <- background_promoter
  assay_args$background_method <- "huber"
  assay_args$background_by <- c("compound", "replicate", "sheet")
}

assay <- do.call(prepare_assay, assay_args)
fit <- fit_destress(
  assay,
  technical = c("replicate", "sheet"),
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
params <- model_parameters(fit)

utils::write.table(
  screen,
  file.path(out_dir, "dryad_weak_stress_auc_input.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
utils::write.table(
  res,
  file.path(out_dir, "dryad_default_pair_results.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
utils::write.table(
  res[res$hit, , drop = FALSE],
  file.path(out_dir, "dryad_default_significant_pairs.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
if (!is.null(params$growth_exponents)) {
  utils::write.table(
    params$growth_exponents,
    file.path(out_dir, "dryad_default_growth_exponents.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}
if (!is.null(params$background_calibration)) {
  utils::write.table(
    params$background_calibration,
    file.path(out_dir, "dryad_default_background_calibration.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

summary <- data.frame(
  metric = c(
    "promoters",
    "conditions_including_reference",
    "reference_condition",
    "background_promoter",
    "well_summaries",
    "tested_pairs",
    "significant_pairs"
  ),
  value = c(
    length(unique(screen$promoter)),
    length(unique(screen$compound)),
    "Standard",
    if (is.null(background_promoter)) "none detected" else background_promoter,
    nrow(screen),
    nrow(res),
    sum(res$hit, na.rm = TRUE)
  )
)
utils::write.table(
  summary,
  file.path(out_dir, "dryad_default_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

p_hist <- ggplot2::ggplot(res, ggplot2::aes(x = specific_pvalue)) +
  ggplot2::geom_histogram(bins = 20, fill = "#3B82F6", color = "white", linewidth = 0.2) +
  ggplot2::facet_wrap(ggplot2::vars(promoter), scales = "free_y") +
  ggplot2::theme_light(base_size = 10) +
  ggplot2::labs(
    title = "Dryad global-regulator reporters: default DStressR p-values",
    x = "Raw p-value",
    y = "Number of promoter-condition pairs"
  )
ggplot2::ggsave(file.path(out_dir, "dryad_default_pvalue_histogram.png"), p_hist, width = 8, height = 5, dpi = 300)
ggplot2::ggsave(file.path(out_dir, "dryad_default_pvalue_histogram.pdf"), p_hist, width = 8, height = 5)

p_volcano <- plot_volcano(
  res,
  effect = "specific_effect",
  padj = "specific_padj_by_promoter",
  title = "Dryad global-regulator reporters: default DStressR volcano",
  label_by = "pair",
  top_n = 12
)
ggplot2::ggsave(file.path(out_dir, "dryad_default_volcano.png"), p_volcano, width = 7.5, height = 5.5, dpi = 300)
ggplot2::ggsave(file.path(out_dir, "dryad_default_volcano.pdf"), p_volcano, width = 7.5, height = 5.5)

plot_dryad_heatmap <- function(tab, title) {
  d <- tab
  d$compound <- factor(d$compound, levels = c("Iron", "Kanamycin", "Tetracycline", "H2O2"))
  d$promoter <- factor(d$promoter, levels = rev(c("Fur", "MarA", "SoxS", "LexA")))
  limit <- max(abs(d$specific_effect), na.rm = TRUE)
  ggplot2::ggplot(d, ggplot2::aes(compound, promoter, fill = specific_effect)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", specific_effect)), size = 3) +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-limit, limit)
    ) +
    ggplot2::theme_light(base_size = 10) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 30, hjust = 1),
      plot.title.position = "plot"
    ) +
    ggplot2::labs(
      title = title,
      x = "Weak-stress condition",
      y = "Reporter",
      fill = "Specific effect"
    )
}

p_heatmap <- plot_dryad_heatmap(res, "Dryad global-regulator reporters: default DStressR effects")
ggplot2::ggsave(file.path(out_dir, "dryad_default_effect_heatmap.png"), p_heatmap, width = 7, height = 4.5, dpi = 300)
ggplot2::ggsave(file.path(out_dir, "dryad_default_effect_heatmap.pdf"), p_heatmap, width = 7, height = 4.5)

assay_args_alpha1 <- assay_args
assay_args_alpha1$growth_exponent <- 1
assay_alpha1 <- do.call(prepare_assay, assay_args_alpha1)
fit_alpha1 <- fit_destress(
  assay_alpha1,
  technical = c("replicate", "sheet"),
  empirical_bayes = TRUE,
  adjustment = "by_promoter",
  interaction = FALSE
)
res_alpha1 <- results(fit_alpha1)
hit_table_alpha1 <- call_hits(
  res_alpha1,
  fdr = 0.05,
  effect = "specific_effect",
  padj = "specific_padj_by_promoter"
)
res_alpha1$hit_class <- hit_table_alpha1$hit
res_alpha1$hit <- res_alpha1$hit_class != "Not DE"
params_alpha1 <- model_parameters(fit_alpha1)

utils::write.table(
  res_alpha1,
  file.path(out_dir, "dryad_alpha1_pair_results.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
utils::write.table(
  res_alpha1[res_alpha1$hit, , drop = FALSE],
  file.path(out_dir, "dryad_alpha1_significant_pairs.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
if (!is.null(params_alpha1$background_calibration)) {
  utils::write.table(
    params_alpha1$background_calibration,
    file.path(out_dir, "dryad_alpha1_background_calibration.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

summary_alpha1 <- data.frame(
  metric = c(
    "promoters",
    "conditions_including_reference",
    "reference_condition",
    "background_promoter",
    "growth_exponent",
    "well_summaries",
    "tested_pairs",
    "significant_pairs"
  ),
  value = c(
    length(unique(screen$promoter)),
    length(unique(screen$compound)),
    "Standard",
    if (is.null(background_promoter)) "none detected" else background_promoter,
    "fixed alpha = 1",
    nrow(screen),
    nrow(res_alpha1),
    sum(res_alpha1$hit, na.rm = TRUE)
  )
)
utils::write.table(
  summary_alpha1,
  file.path(out_dir, "dryad_alpha1_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

p_hist_alpha1 <- ggplot2::ggplot(res_alpha1, ggplot2::aes(x = specific_pvalue)) +
  ggplot2::geom_histogram(bins = 20, fill = "#009E73", color = "white", linewidth = 0.2) +
  ggplot2::facet_wrap(ggplot2::vars(promoter), scales = "free_y") +
  ggplot2::theme_light(base_size = 10) +
  ggplot2::labs(
    title = "Dryad global-regulator reporters: DStressR p-values with alpha = 1",
    x = "Raw p-value",
    y = "Number of promoter-condition pairs"
  )
ggplot2::ggsave(file.path(out_dir, "dryad_alpha1_pvalue_histogram.png"), p_hist_alpha1, width = 8, height = 5, dpi = 300)
ggplot2::ggsave(file.path(out_dir, "dryad_alpha1_pvalue_histogram.pdf"), p_hist_alpha1, width = 8, height = 5)

p_volcano_alpha1 <- plot_volcano(
  res_alpha1,
  effect = "specific_effect",
  padj = "specific_padj_by_promoter",
  title = "Dryad global-regulator reporters: DStressR volcano with alpha = 1",
  label_by = "pair",
  top_n = 12
)
ggplot2::ggsave(file.path(out_dir, "dryad_alpha1_volcano.png"), p_volcano_alpha1, width = 7.5, height = 5.5, dpi = 300)
ggplot2::ggsave(file.path(out_dir, "dryad_alpha1_volcano.pdf"), p_volcano_alpha1, width = 7.5, height = 5.5)

p_heatmap_alpha1 <- plot_dryad_heatmap(res_alpha1, "Dryad global-regulator reporters: DStressR effects with alpha = 1")
ggplot2::ggsave(file.path(out_dir, "dryad_alpha1_effect_heatmap.png"), p_heatmap_alpha1, width = 7, height = 4.5, dpi = 300)
ggplot2::ggsave(file.path(out_dir, "dryad_alpha1_effect_heatmap.pdf"), p_heatmap_alpha1, width = 7, height = 4.5)

comparison <- merge(
  res[, c("promoter", "compound", "specific_effect", "specific_pvalue", "specific_padj_by_promoter", "hit")],
  res_alpha1[, c("promoter", "compound", "specific_effect", "specific_pvalue", "specific_padj_by_promoter", "hit")],
  by = c("promoter", "compound"),
  suffixes = c("_estimated_alpha", "_alpha1"),
  all = TRUE,
  sort = FALSE
)
utils::write.table(
  comparison,
  file.path(out_dir, "dryad_estimated_alpha_vs_alpha1_comparison.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Wrote Dryad default and alpha=1 DStressR outputs to: ", out_dir)
