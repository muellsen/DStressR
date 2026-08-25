source(file.path("analysis", "_helpers.R"))
load_destress_package()

load_binsfeld_paper_data()

out_dir <- analysis_output_dir("binsfeld_variance_diagnostics")

wt_auc <- binsfeld_reporter_auc[
  binsfeld_reporter_auc$strain == "WT" &
    binsfeld_reporter_auc$removed == "No",
]

fit_current <- function(use_ev_control) {
  assay_data <- if (isTRUE(use_ev_control)) {
    wt_auc
  } else {
    wt_auc[wt_auc$promoter != "EVC", , drop = FALSE]
  }
  args <- list(
    data = assay_data,
    reporter = "promoter",
    perturbation = "compound",
    control = "Water",
    lux = "lux_auc",
    growth = "od_auc",
    growth_exponent = "estimate",
    batch = "dose_level",
    replicate = "replicate",
    growth_covariates = "replicate",
    numeric_covariates = "dose_level"
  )
  if (isTRUE(use_ev_control)) {
    args$background_reporter <- "EVC"
    args$background_method <- "huber"
    args$background_by <- c("compound", "dose_level", "replicate")
  }
  assay <- do.call(prepare_assay, args)
  fit <- fit_destress(
    assay,
    technical = c("replicate", "dose_level"),
    empirical_bayes = TRUE,
    adjustment = "by_reporter",
    interaction = FALSE
  )
  res <- results(fit)
  calls <- call_hits(
    res,
    fdr = 0.05,
    effect = "specific_effect",
    padj = "specific_padj_by_reporter"
  )
  res$hit <- calls$hit != "Not DE"
  res$hit_class <- calls$hit
  res$pair_id <- paste(res$promoter, res$compound, sep = "__")
  res
}

current_no_ev <- fit_current(FALSE)
current_ev <- fit_current(TRUE)

old_comparison <- read.delim(
  analysis_path(
    "analysis", "outputs", "binsfeld_evc_calibrated",
    "evc_huber_comparison_to_binsfeld_and_default.tsv"
  ),
  check.names = FALSE
)
old_comparison$pair_id <- paste(old_comparison$promoter, old_comparison$compound, sep = "__")
old_no_ev_ids <- old_comparison$pair_id[old_comparison$modeled_hit]
old_ev_ids <- old_comparison$pair_id[old_comparison$evc_huber_hit]
binsfeld_ids <- old_comparison$pair_id[old_comparison$binsfeld_hit]

old_ev <- read.delim(
  analysis_path(
    "analysis", "outputs", "binsfeld_evc_calibrated",
    "evc_huber_significant_pairs.tsv"
  ),
  check.names = FALSE
)
old_ev$pair_id <- paste(old_ev$promoter, old_ev$compound, sep = "__")

comparison_table <- function(current, old_ids, label) {
  current_ids <- current$pair_id[current$hit]
  gained <- setdiff(current_ids, old_ids)
  lost <- setdiff(old_ids, current_ids)
  changed <- current[current$pair_id %in% c(gained, lost), , drop = FALSE]
  if (nrow(changed) == 0) {
    return(data.frame(
      workflow = character(),
      status = character(),
      promoter = character(),
      compound = character(),
      specific_effect = numeric(),
      specific_pvalue = numeric(),
      specific_padj_by_reporter = numeric(),
      hit_class = character(),
      stringsAsFactors = FALSE
    ))
  }
  changed$status <- ifelse(changed$pair_id %in% gained, "Gained", "Lost")
  changed$workflow <- label
  changed[, c(
    "workflow", "status", "promoter", "compound", "specific_effect",
    "specific_pvalue", "specific_padj_by_reporter", "hit_class"
  ), drop = FALSE]
}

changed <- rbind(
  comparison_table(current_no_ev, old_no_ev_ids, "DStressR without EV control"),
  comparison_table(current_ev, old_ev_ids, "DStressR with EV control")
)

summary <- data.frame(
  workflow = c("DStressR without EV control", "DStressR with EV control"),
  old_hits = c(length(old_no_ev_ids), length(old_ev_ids)),
  current_hits = c(sum(current_no_ev$hit), sum(current_ev$hit)),
  gained = c(
    sum(changed$workflow == "DStressR without EV control" & changed$status == "Gained"),
    sum(changed$workflow == "DStressR with EV control" & changed$status == "Gained")
  ),
  lost = c(
    sum(changed$workflow == "DStressR without EV control" & changed$status == "Lost"),
    sum(changed$workflow == "DStressR with EV control" & changed$status == "Lost")
  ),
  binsfeld_overlap = c(
    length(intersect(current_no_ev$pair_id[current_no_ev$hit], binsfeld_ids)),
    length(intersect(current_ev$pair_id[current_ev$hit], binsfeld_ids))
  ),
  stringsAsFactors = FALSE
)

pvalue_summary <- function(current, old_effect_col, old_p_col, old_q_col, label) {
  old <- old_comparison[, c("pair_id", old_effect_col, old_p_col, old_q_col), drop = FALSE]
  names(old) <- c("pair_id", "old_effect", "old_pvalue", "old_qvalue")
  cur <- current[, c("pair_id", "specific_effect", "specific_pvalue", "specific_padj_by_reporter"), drop = FALSE]
  names(cur) <- c("pair_id", "current_effect", "current_pvalue", "current_qvalue")
  merged <- merge(old, cur, by = "pair_id", all = FALSE, sort = FALSE)
  data.frame(
    workflow = label,
    max_abs_effect_difference = max(abs(merged$current_effect - merged$old_effect), na.rm = TRUE),
    median_log10_pvalue_change = stats::median(
      log10(pmax(merged$current_pvalue, .Machine$double.xmin)) -
        log10(pmax(merged$old_pvalue, .Machine$double.xmin)),
      na.rm = TRUE
    ),
    max_abs_log10_pvalue_change = max(abs(
      log10(pmax(merged$current_pvalue, .Machine$double.xmin)) -
        log10(pmax(merged$old_pvalue, .Machine$double.xmin))
    ), na.rm = TRUE),
    median_qvalue_difference = stats::median(merged$current_qvalue - merged$old_qvalue, na.rm = TRUE),
    max_abs_qvalue_difference = max(abs(merged$current_qvalue - merged$old_qvalue), na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

pvalue_summary_table <- rbind(
  pvalue_summary(
    current_no_ev,
    "modeled_effect",
    "modeled_pvalue",
    "modeled_padj_by_reporter",
    "DStressR without EV control"
  ),
  pvalue_summary(
    current_ev,
    "evc_huber_effect",
    "evc_huber_pvalue",
    "evc_huber_padj_by_reporter",
    "DStressR with EV control"
  )
)

write.table(
  summary,
  file.path(out_dir, "ecoli_moderation_update_hit_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  changed,
  file.path(out_dir, "ecoli_moderation_update_changed_hits.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  pvalue_summary_table,
  file.path(out_dir, "ecoli_moderation_update_pvalue_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

print(summary)
print(pvalue_summary_table)
if (nrow(changed)) {
  print(changed[order(changed$workflow, changed$status, changed$promoter, changed$compound), ], row.names = FALSE)
}
