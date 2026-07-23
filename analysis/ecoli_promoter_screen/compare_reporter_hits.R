source(file.path("analysis", "_helpers.R"))
load_destress_package()

load(analysis_path("data", "binsfeld_reporter_data.rda"))

out_dir <- analysis_output_dir("binsfeld")

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
  batch = "concentration_index",
  replicate = "replicate"
)

fit <- fit_destress(
  assay,
  preset = "model",
  technical = c("replicate", "concentration_index"),
  empirical_bayes = TRUE,
  adjustment = "by_promoter",
  interaction = FALSE
)

destress_results <- results(fit)
destress_hits <- call_hits(
  destress_results,
  fdr = 0.05,
  effect = "specific_effect",
  padj = "specific_padj_by_promoter"
)
destress_hits <- destress_hits[destress_hits$hit != "Not DE", ]

wt_z <- binsfeld_reporter_scores[
  binsfeld_reporter_scores$strain == "WT" &
    binsfeld_reporter_scores$statistic == "Z_scores",
]

author_hits <- do.call(rbind, lapply(sort(unique(wt_z$promoter)), function(promoter) {
  if (promoter == "EVC") {
    return(NULL)
  }
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
      pvalue = pvalue,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out$padj <- stats::p.adjust(out$pvalue, method = "BH")
  out$hit <- is.finite(out$padj) &
    out$padj < 0.05 &
    abs(out$mean_z) > 1 &
    out$compound != "Water"
  out
}))
author_hits <- author_hits[author_hits$hit, ]

hit_key <- function(x) paste(x$promoter, x$compound, sep = "\r")
author_key <- hit_key(author_hits)
destress_key <- hit_key(destress_hits)

comparison <- merge(
  author_hits[, c("promoter", "compound", "mean_z", "padj")],
  destress_hits[, c("promoter", "compound", "specific_effect", "specific_padj_by_promoter", "hit")],
  by = c("promoter", "compound"),
  all = TRUE
)
names(comparison)[names(comparison) == "padj"] <- "binsfeld_padj"
comparison$binsfeld_hit <- hit_key(comparison) %in% author_key
comparison$destress_hit <- hit_key(comparison) %in% destress_key

summary <- data.frame(
  analysis = c("Binsfeld_Wilcoxon_Z", "DStressR_default_modeled_response", "overlap"),
  hits = c(length(author_key), length(destress_key), sum(author_key %in% destress_key)),
  stringsAsFactors = FALSE
)

utils::write.table(
  summary,
  file.path(out_dir, "binsfeld_destress_hit_summary.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
utils::write.table(
  comparison[order(comparison$promoter, comparison$compound), ],
  file.path(out_dir, "binsfeld_destress_hit_comparison.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

print(summary)
