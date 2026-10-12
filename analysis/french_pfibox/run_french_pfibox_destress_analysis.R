#!/usr/bin/env Rscript

# Replicate-level DStressR analysis scaffold for the French et al. PFIbox
# promoter-reporter data. The public object stores colony-level GFP traces and
# constitutive-normalized traces. We use those summaries to construct a
# Binsfeld-like response analysis with and without reporter-specific response
# exponent estimation.

source(file.path("analysis", "_helpers.R"))
load_destress_package()

stopifnot(requireNamespace("ggplot2", quietly = TRUE))

source_dir <- analysis_path("data-raw", "french_pfibox", "source")
processed_dir <- analysis_path("data-raw", "french_pfibox", "processed")
out_dir <- analysis_output_dir("french_pfibox")
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

auc_file <- file.path(processed_dir, "pfibox_colony_auc_for_destress.tsv")
proxy_map_file <- file.path(processed_dir, "pfibox_ampicillin0_proxy_maps.tsv")
rebuild_auc <- identical(tolower(Sys.getenv("PFIBOX_REBUILD_AUC", "false")), "true")
reporter_limit <- suppressWarnings(as.integer(Sys.getenv("PFIBOX_REPORTER_LIMIT", NA_character_)))
perturbation_limit <- suppressWarnings(as.integer(Sys.getenv("PFIBOX_PERTURBATION_LIMIT", NA_character_)))
run_tag <- "full"

parse_alpha_bounds <- function(x) {
  if (!nzchar(x) || tolower(x) %in% c("none", "null", "unbounded")) {
    return(NULL)
  }
  parts <- strsplit(x, ",", fixed = TRUE)[[1]]
  if (length(parts) != 2) {
    stop("PFIBOX_ALPHA_BOUNDS must be `none` or two comma-separated numbers.", call. = FALSE)
  }
  bounds <- as.numeric(parts)
  if (any(!is.finite(bounds)) || bounds[1] >= bounds[2]) {
    stop("PFIBOX_ALPHA_BOUNDS must define finite increasing bounds.", call. = FALSE)
  }
  bounds
}

alpha_bounds <- parse_alpha_bounds(Sys.getenv("PFIBOX_ALPHA_BOUNDS", "none"))

gfp_file <- file.path(source_dir, "gfpData.Rdata")

parse_treatment <- function(x) {
  m <- regexec("^(.*)-([0-9]+)$", x)
  r <- regmatches(x, m)[[1]]
  if (length(r) == 3) {
    data.frame(
      treatment = x,
      drug = r[2],
      concentration_code = r[3],
      perturbation = if (identical(r[3], "0")) "Control" else x,
      is_control = identical(r[3], "0"),
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      treatment = x,
      drug = x,
      concentration_code = NA_character_,
      perturbation = x,
      is_control = FALSE,
      stringsAsFactors = FALSE
    )
  }
}

auc_one_treatment <- function(treatment_name) {
  treatment <- gfpData[[treatment_name]]
  if (is.null(treatment$raw.loess)) {
    return(NULL)
  }
  raw <- as.matrix(treatment$raw.loess)
  legend <- as.data.frame(gfpData$legend, stringsAsFactors = FALSE)
  if (nrow(raw) != nrow(legend)) {
    stop("Treatment `", treatment_name, "` does not match the common PFIbox legend.", call. = FALSE)
  }
  proxy_map <- utils::read.delim(proxy_map_file, check.names = FALSE)
  keep <- !is.na(legend$Promoter) & nzchar(legend$Promoter) &
    !legend$Promoter %in% c("Empty", "U66", "U139")
  position <- seq_len(nrow(legend))
  keep_map <- proxy_map[match(position[keep], proxy_map$position), , drop = FALSE]
  complete_map <- !is.na(keep_map$constit_position) & !is.na(keep_map$background_position)
  keep_idx <- which(keep)[complete_map]
  keep_map <- keep_map[complete_map, , drop = FALSE]
  legend <- legend[keep_idx, , drop = FALSE]

  reporter_trace <- raw[keep_idx, , drop = FALSE] -
    raw[keep_map$background_position, , drop = FALSE]
  biomass_proxy_trace <- raw[keep_map$constit_position, , drop = FALSE] -
    raw[keep_map$background_position, , drop = FALSE]
  reporter_trace[!is.finite(reporter_trace) | reporter_trace <= 0] <- NA_real_
  biomass_proxy_trace[!is.finite(biomass_proxy_trace) | biomass_proxy_trace <= 0] <- NA_real_

  fluorescence_auc <- rowSums(reporter_trace, na.rm = TRUE)
  biomass_proxy_auc <- rowSums(biomass_proxy_trace, na.rm = TRUE)
  biomass_proxy_auc[!is.finite(biomass_proxy_auc) | biomass_proxy_auc <= 0] <- NA_real_
  fluorescence_auc[!is.finite(fluorescence_auc) | fluorescence_auc <= 0] <- NA_real_

  meta <- treatment_metadata[treatment_metadata$treatment == treatment_name, , drop = FALSE]
  d <- data.frame(
    reporter = as.character(legend$Promoter),
    perturbation = meta$perturbation,
    treatment = treatment_name,
    drug = meta$drug,
    concentration_code = as.integer(meta$concentration_code),
    is_control = meta$is_control,
    colony_position = keep_idx,
    constitutive_position = keep_map$constit_position,
    background_position = keep_map$background_position,
    plate1536 = legend$Plate1536,
    row1536 = legend$Row1536,
    column1536 = legend$Column1536,
    plate6144 = legend$Plate6144,
    row6144 = legend$Row6144,
    column6144 = legend$Column6144,
    replicate = ave(
      seq_along(legend$Promoter),
      legend$Promoter,
      FUN = seq_along
    ),
    fluorescence_auc = fluorescence_auc,
    biomass_proxy_auc = biomass_proxy_auc,
    stringsAsFactors = FALSE
  )
  d[is.finite(d$fluorescence_auc) & d$fluorescence_auc > 0 &
      is.finite(d$biomass_proxy_auc) & d$biomass_proxy_auc > 0, , drop = FALSE]
}

if (file.exists(auc_file) && !rebuild_auc) {
  message("Reading cached PFIbox AUC table: ", auc_file)
  auc_tab <- utils::read.delim(auc_file, check.names = FALSE)
} else {
  if (!file.exists(gfp_file)) {
    stop(
      "Missing multi-treatment PFIbox object: ", gfp_file, "\n",
      "Download the French et al. deposited data or create `gfpData.Rdata` with ",
      "the PFIbox analysis workflow. The GitHub sample contains only one untreated ",
      "condition and is not sufficient for DStressR perturbation analysis.",
      call. = FALSE
    )
  }

  load(gfp_file)
  if (!exists("gfpData") || !is.list(gfpData)) {
    stop("`gfpData.Rdata` must contain a list object named `gfpData`.", call. = FALSE)
  }
  if (!file.exists(proxy_map_file)) {
    stop(
      "Missing PFIbox proxy map: ", proxy_map_file, "\n",
      "Run the PFIbox audit/reconstruction step before rebuilding the AUC table.",
      call. = FALSE
    )
  }

  treatment_names <- setdiff(names(gfpData), "legend")
  treatment_metadata <- do.call(rbind, lapply(treatment_names, parse_treatment))
  message("Building colony-level AUC table from ", length(treatment_names), " treatments.")
  auc_parts <- lapply(treatment_names, auc_one_treatment)
  auc_parts <- Filter(Negate(is.null), auc_parts)
  if (!length(auc_parts)) {
    stop("No treatments with `raw.loess` and `constit.normalized` were found.", call. = FALSE)
  }
  auc_tab <- do.call(rbind, auc_parts)

  utils::write.table(
    auc_tab,
    auc_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  rm(gfpData)
  gc()
}

if (is.finite(reporter_limit) && reporter_limit > 0) {
  reporters <- sort(unique(auc_tab$reporter))
  keep_reporters <- reporters[seq_len(min(reporter_limit, length(reporters)))]
  auc_tab <- auc_tab[auc_tab$reporter %in% keep_reporters, , drop = FALSE]
  run_tag <- paste0(run_tag, "_reporters", length(keep_reporters))
  message("Using reporter subset for smoke test: ", length(keep_reporters), " reporters.")
}

if (is.finite(perturbation_limit) && perturbation_limit > 0) {
  perturbations <- sort(setdiff(unique(auc_tab$perturbation), "Control"))
  keep_perturbations <- perturbations[seq_len(min(perturbation_limit, length(perturbations)))]
  auc_tab <- auc_tab[auc_tab$perturbation %in% c("Control", keep_perturbations), , drop = FALSE]
  run_tag <- paste0(run_tag, "_perturbations", length(keep_perturbations))
  message("Using perturbation subset for smoke test: ", length(keep_perturbations), " perturbations plus Control.")
}

fit_one <- function(label, growth_exponent) {
  message("Fitting PFIbox DStressR model: ", label)
  growth_fit <- NULL
  if (is.character(growth_exponent) && identical(growth_exponent, "estimate")) {
    growth_fit <- estimate_growth_exponents(
      auc_tab,
      reporter = "reporter",
      perturbation = "perturbation",
      lux = "fluorescence_auc",
      growth = "biomass_proxy_auc",
      controls = "Control",
      covariates = "drug",
      alpha_bounds = alpha_bounds
    )
    growth_exponent <- stats::setNames(growth_fit$alpha_shrunk, growth_fit$reporter)
  }
  assay <- prepare_assay(
    auc_tab,
    reporter = "reporter",
    perturbation = "perturbation",
    control = "Control",
    lux = "fluorescence_auc",
    growth = "biomass_proxy_auc",
    growth_exponent = growth_exponent,
    control_values = "Control",
    growth_covariates = "drug"
  )
  fit <- fit_destress(
    assay,
    empirical_bayes = TRUE,
    adjustment = "by_reporter",
    interaction = FALSE
  )
  res <- results(fit)
  res$analysis <- label
  params <- model_parameters(fit)

  prefix <- paste0("pfibox_", run_tag, "_", label)
  utils::write.table(
    res,
    file.path(out_dir, paste0(prefix, "_pair_results.tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  if (!is.null(params$growth_exponents)) {
    utils::write.table(
      params$growth_exponents,
      file.path(out_dir, paste0(prefix, "_growth_exponents.tsv")),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  }
  if (!is.null(growth_fit)) {
    utils::write.table(
      growth_fit,
      file.path(out_dir, paste0(prefix, "_growth_exponents.tsv")),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  }
  list(assay = assay, fit = fit, results = res, params = params, growth_fit = growth_fit)
}

fits <- list(
  alpha_estimated = fit_one("alpha_estimated", "estimate"),
  alpha_1 = fit_one("alpha_1", 1)
)

combined <- do.call(rbind, lapply(fits, `[[`, "results"))
combined_file <- file.path(out_dir, paste0("pfibox_", run_tag, "_alpha_comparison_pair_results.tsv"))
utils::write.table(
  combined,
  combined_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

analysis_labels <- c(
  alpha_estimated = "Estimated response exponent",
  alpha_1 = expression(alpha[a] == 1)
)

pvalue_panel <- combined[
  is.finite(combined$specific_pvalue) &
    is.finite(combined$specific_padj_by_reporter),
  ,
  drop = FALSE
]
pvalue_panel$analysis_label <- factor(
  pvalue_panel$analysis,
  levels = names(analysis_labels),
  labels = as.character(analysis_labels)
)
p_pvalue <- ggplot2::ggplot(pvalue_panel, ggplot2::aes(x = specific_pvalue)) +
  ggplot2::geom_histogram(
    bins = 50,
    fill = "#4E79A7",
    color = "white",
    linewidth = 0.2
  ) +
  ggplot2::facet_wrap(ggplot2::vars(analysis_label), nrow = 1) +
  ggplot2::labs(
    x = "Specific-effect p-value",
    y = "Number of reporter-perturbation pairs"
  ) +
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(
    strip.background = ggplot2::element_rect(fill = "#111827", color = NA),
    strip.text = ggplot2::element_text(color = "white", face = "bold")
  )
ggplot2::ggsave(
  file.path(out_dir, paste0("pfibox_", run_tag, "_specific_pvalue_histograms.png")),
  p_pvalue,
  width = 7.4,
  height = 3.2,
  dpi = 300,
  bg = "white"
)
ggplot2::ggsave(
  file.path(out_dir, paste0("pfibox_", run_tag, "_specific_pvalue_histograms.pdf")),
  p_pvalue,
  width = 7.4,
  height = 3.2,
  bg = "white"
)

make_volcano <- function(tab, label, effect, pvalue, padj, xlab) {
  plot_volcano(
    tab,
    effect = effect,
    pvalue = pvalue,
    padj = padj,
    y = "pvalue",
    reporter = "reporter",
    perturbation = "perturbation",
    fdr = 0.05,
    top_n = 12,
    title = NULL,
    subtitle = "",
    xlab = xlab,
    ylab = expression(-log[10]("p-value")),
    call_labels = c("Not called", "FDR < 0.05"),
    max_label_chars = 32,
    point_size = 1.1,
    point_alpha = 0.75,
    hit_color = "#D95F02",
    background_color = "#CBD5E1"
  ) +
    ggplot2::labs(caption = label) +
    ggplot2::theme(
      plot.caption = ggplot2::element_text(
        hjust = 0,
        face = "bold",
        color = "#111827",
        size = 10
      ),
      legend.position = c(0.78, 0.86),
      legend.background = ggplot2::element_rect(fill = "white", color = "#E5E7EB")
    )
}

for (nm in names(fits)) {
  tab <- fits[[nm]]$results
  label <- if (nm == "alpha_estimated") "Estimated response exponent" else "alpha[a] = 1"
  v_specific <- make_volcano(
    tab,
    label = label,
    effect = "specific_effect",
    pvalue = "specific_pvalue",
    padj = "specific_padj_by_reporter",
    xlab = expression(hat(delta)[aj])
  )
  v_total <- make_volcano(
    tab,
    label = label,
    effect = "total_effect",
    pvalue = "total_pvalue",
    padj = "total_padj_by_reporter",
    xlab = expression(hat(Delta*y)[aj]^tot)
  )
  ggplot2::ggsave(
    file.path(out_dir, paste0("pfibox_", run_tag, "_", nm, "_specific_volcano.png")),
    v_specific,
    width = 5.8,
    height = 4.4,
    dpi = 300,
    bg = "white"
  )
  ggplot2::ggsave(
    file.path(out_dir, paste0("pfibox_", run_tag, "_", nm, "_total_volcano.png")),
    v_total,
    width = 5.8,
    height = 4.4,
    dpi = 300,
    bg = "white"
  )
}

specific_compare <- merge(
  fits$alpha_estimated$results[, c("reporter", "perturbation", "specific_effect", "specific_pvalue")],
  fits$alpha_1$results[, c("reporter", "perturbation", "specific_effect", "specific_pvalue")],
  by = c("reporter", "perturbation"),
  suffixes = c("_alpha_estimated", "_alpha_1")
)
p_compare <- ggplot2::ggplot(
  specific_compare,
  ggplot2::aes(x = specific_effect_alpha_1, y = specific_effect_alpha_estimated)
) +
  ggplot2::geom_abline(slope = 1, intercept = 0, color = "#6B7280", linewidth = 0.35) +
  ggplot2::geom_point(color = "#2A9D8F", alpha = 0.72, size = 1.25) +
  ggplot2::labs(
    x = expression(hat(delta)[aj]~","~alpha[a] == 1),
    y = expression(hat(delta)[aj]~","~hat(alpha)[a])
  ) +
  ggplot2::theme_classic(base_size = 11)
ggplot2::ggsave(
  file.path(out_dir, paste0("pfibox_", run_tag, "_specific_effect_alpha_comparison.png")),
  p_compare,
  width = 4.6,
  height = 4.2,
  dpi = 300,
  bg = "white"
)
ggplot2::ggsave(
  file.path(out_dir, paste0("pfibox_", run_tag, "_specific_effect_alpha_comparison.pdf")),
  p_compare,
  width = 4.6,
  height = 4.2,
  bg = "white"
)

alpha_est <- fits$alpha_estimated$growth_fit
if (is.null(alpha_est)) {
  alpha_est <- fits$alpha_estimated$params$growth_exponents
}
if (!is.null(alpha_est)) {
  p_alpha <- ggplot2::ggplot(alpha_est, ggplot2::aes(x = alpha_shrunk)) +
    ggplot2::geom_histogram(bins = 50, fill = "#2A9D8F", color = "white", linewidth = 0.2) +
    ggplot2::labs(
      x = expression(hat(alpha)[a]),
      y = "Number of reporters"
    ) +
    ggplot2::theme_classic(base_size = 11)
  ggplot2::ggsave(
    file.path(out_dir, paste0("pfibox_", run_tag, "_estimated_response_exponents.png")),
    p_alpha,
    width = 5.4,
    height = 3.8,
    dpi = 300,
    bg = "white"
  )
  ggplot2::ggsave(
    file.path(out_dir, paste0("pfibox_", run_tag, "_estimated_response_exponents.pdf")),
    p_alpha,
    width = 5.4,
    height = 3.8,
    bg = "white"
  )
}

message("Reference perturbation: Control")
message("Wrote colony-level AUC table: ", auc_file)
message("Wrote DStressR alpha comparison: ", combined_file)
