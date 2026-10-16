source(file.path("analysis", "_helpers.R"))
source(file.path("analysis", "_paper_artifacts.R"))
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

load_binsfeld_paper_data()

out_dir <- analysis_output_dir("binsfeld_evc_calibrated")
workflow_comparison_dir <- analysis_output_dir("binsfeld_workflow_comparison")

safe_neglog10 <- function(x) -log10(pmax(as.numeric(x), .Machine$double.xmin))

wt_auc <- binsfeld_reporter_auc[
  binsfeld_reporter_auc$strain == "WT" &
    binsfeld_reporter_auc$removed == "No",
]

assay <- prepare_assay(
  wt_auc,
  reporter = "promoter",
  perturbation = "compound",
  control = "Water",
  lux = "lux_auc",
  growth = "od_auc",
  growth_exponent = "estimate",
  batch = "dose_level",
  replicate = "replicate",
  growth_covariates = "replicate",
  numeric_covariates = "dose_level",
  background_reporter = "EVC",
  background_method = "huber",
  background_by = c("compound", "dose_level", "replicate")
)

fit <- fit_destress(
  assay,
  technical = c("replicate", "dose_level"),
  empirical_bayes = TRUE,
  adjustment = "by_reporter",
  interaction = FALSE
)

saveRDS(fit, file.path(out_dir, "evc_huber_fit.rds"))
res <- results(fit)
utils::write.table(res, file.path(out_dir, "evc_huber_full_results.tsv"),
                   sep = "\t", quote = FALSE, row.names = FALSE)
res$evc_huber_hit_class <- call_hits(
  res,
  fdr = 0.05,
  effect = "specific_effect",
  padj = "specific_padj_by_reporter"
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
  "reporter", "perturbation", "total_effect", "global_effect", "low_rank_effect",
  "specific_effect", "specific_se", "specific_statistic", "specific_pvalue",
  "specific_padj_by_reporter", "evc_huber_hit", "evc_huber_hit_class"
)]
names(pair_results)[names(pair_results) == "reporter"] <- "promoter"
names(pair_results)[names(pair_results) == "perturbation"] <- "compound"
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

comparison_file <- file.path(workflow_comparison_dir, "binsfeld_workflow_comparison_all_pair_comparison.tsv")
if (file.exists(comparison_file)) {
  comparison <- utils::read.delim(comparison_file, check.names = FALSE)
  evc <- pair_results[, c(
    "promoter", "compound", "specific_effect", "specific_pvalue",
    "specific_padj_by_reporter", "evc_huber_hit", "evc_huber_hit_class"
  )]
  names(evc) <- c(
    "promoter", "compound", "evc_huber_effect", "evc_huber_pvalue",
    "evc_huber_padj_by_reporter", "evc_huber_hit", "evc_huber_hit_class"
  )
  comparison <- merge(comparison, evc, by = c("promoter", "compound"), all.x = TRUE, sort = FALSE)
  comparison$evc_huber_hit[is.na(comparison$evc_huber_hit)] <- FALSE
  comparison$binsfeld_hit <- as.logical(comparison$binsfeld_hit)
  comparison$estimated_exponent_hit <- as.logical(comparison$estimated_exponent_hit)
  comparison$fixed_exponent_hit <- as.logical(comparison$fixed_exponent_hit)

  comparison$evc_huber_class <- apply(
    comparison[, c("binsfeld_hit", "estimated_exponent_hit", "evc_huber_hit")],
    1,
    function(x) {
      names <- c("Binsfeld reference", "DStressR without EVC", "DStressR with EVC")[as.logical(x)]
      if (length(names) == 0) {
        "None"
      } else {
        paste(names, collapse = " + ")
      }
    }
  )
  comparison$neglog10_binsfeld <- safe_neglog10(comparison$binsfeld_pvalue)
  comparison$neglog10_binsfeld_padj <- safe_neglog10(comparison$binsfeld_padj)
  comparison$neglog10_estimated_exponent <- safe_neglog10(comparison$estimated_exponent_pvalue)
  comparison$neglog10_evc_huber <- safe_neglog10(comparison$evc_huber_pvalue)

  utils::write.table(
    comparison,
    file.path(out_dir, "evc_huber_comparison_to_binsfeld_and_default.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  primary_dstressr_only <- comparison[
    !comparison$binsfeld_hit & (comparison$estimated_exponent_hit | comparison$evc_huber_hit),
    ,
    drop = FALSE
  ]
  literature_file <- file.path(workflow_comparison_dir, "dstressr_only_literature_support.tsv")
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
      "DStressR without EVC significant pairs",
      "DStressR with EVC significant pairs",
      "Reference and DStressR with EVC overlap",
      "DStressR workflows overlap",
      "All three overlap",
      "DStressR with EVC only vs reference/without EVC"
    ),
    count = c(
      nrow(comparison),
      sum(comparison$binsfeld_hit, na.rm = TRUE),
      sum(comparison$estimated_exponent_hit, na.rm = TRUE),
      sum(comparison$evc_huber_hit, na.rm = TRUE),
      sum(comparison$binsfeld_hit & comparison$evc_huber_hit, na.rm = TRUE),
      sum(comparison$estimated_exponent_hit & comparison$evc_huber_hit, na.rm = TRUE),
      sum(comparison$binsfeld_hit & comparison$estimated_exponent_hit & comparison$evc_huber_hit, na.rm = TRUE),
      sum(!comparison$binsfeld_hit & !comparison$estimated_exponent_hit & comparison$evc_huber_hit, na.rm = TRUE)
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
    binsfeld_only = sum(comparison$binsfeld_hit & !comparison$estimated_exponent_hit & !comparison$evc_huber_hit, na.rm = TRUE),
    estimated_exponent_only = sum(!comparison$binsfeld_hit & comparison$estimated_exponent_hit & !comparison$evc_huber_hit, na.rm = TRUE),
    evc_huber_only = sum(!comparison$binsfeld_hit & !comparison$estimated_exponent_hit & comparison$evc_huber_hit, na.rm = TRUE),
    binsfeld_estimated_exponent = sum(comparison$binsfeld_hit & comparison$estimated_exponent_hit & !comparison$evc_huber_hit, na.rm = TRUE),
    binsfeld_evc_huber = sum(comparison$binsfeld_hit & !comparison$estimated_exponent_hit & comparison$evc_huber_hit, na.rm = TRUE),
    estimated_exponent_evc_huber = sum(!comparison$binsfeld_hit & comparison$estimated_exponent_hit & comparison$evc_huber_hit, na.rm = TRUE),
    all_three = sum(comparison$binsfeld_hit & comparison$estimated_exponent_hit & comparison$evc_huber_hit, na.rm = TRUE)
  )
  set_counts <- c(
    binsfeld = sum(comparison$binsfeld_hit, na.rm = TRUE),
    estimated_exponent = sum(comparison$estimated_exponent_hit, na.rm = TRUE),
    evc_huber = sum(comparison$evc_huber_hit, na.rm = TRUE)
  )
  utils::write.table(
    data.frame(region = names(region_counts), count = as.integer(region_counts)),
    file.path(out_dir, "evc_huber_region_counts.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  write_paper_artifact(file.path(out_dir, "binsfeld_discovery.rds"),
    results = list(with_EVC = res), parameters = params,
    annotations = list(comparison = comparison, region_counts = region_counts,
                       set_counts = set_counts, summary = summary))
}
message("Wrote Binsfeld EVC analysis and discovery artifact: ", out_dir)
