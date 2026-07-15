#!/usr/bin/env Rscript

# Template for regenerating the local Campylobacter DStressR standard-model
# package output used by the downstream manuscript analysis scripts.
#
# This script contains no data. It expects the proprietary Campylobacter
# workflow data to exist locally. Set DSTRESSR_DATA_ROOT to the directory that
# contains 03-hit_determination/expression_table.tsv.gz, or place that workflow
# data next to this repository as described in analysis/_helpers.R.
#
# Run from the repository root:
#   Rscript scripts/export_campy_standard_model_template.R

source(file.path("analysis", "_helpers.R"))
load_destress_package()

data_root <- analysis_data_root()
out_dir <- package_results_dir()
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

expr <- read_tsv_base(file.path(data_root, "03-hit_determination", "expression_table.tsv.gz"))

# Match the promoter set used in the current Campylobacter manuscript figures.
expr <- expr[!(expr$promoter %in% c("PCJnc20", "PCjas704")), , drop = FALSE]

# DMSO is coded as one common reference level; non-DMSO compounds retain their
# library-position code. Library plate and replicate enter as technical factors.
expr$libplate <- sub("_.*$", "", expr$srn_code)
expr$compound_model <- ifelse(expr$ProductName == "DMSO", "DMSO", expr$srn_code)

assay <- prepare_assay(
  expr,
  promoter = "promoter",
  compound = "compound_model",
  control = "DMSO",
  lux = "lux_auc_until16h",
  growth = "od_at_16h",
  growth_exponent = "estimate",
  plate = "libplate",
  replicate = "replicate"
)

fit_standard <- fit_destress(
  assay,
  technical = c("libplate", "replicate"),
  empirical_bayes = FALSE,
  interaction = FALSE,
  adjustment = "global"
)

res <- results(fit_standard)

standard_pairs <- data.frame(
  promoter = res$promoter,
  compound = res$compound,
  effect = res$specific_effect,
  pvalue = res$specific_pvalue,
  padj_global = res$specific_padj_global,
  padj_by_promoter = res$specific_padj_by_promoter,
  stringsAsFactors = FALSE
)

out_file <- file.path(out_dir, "destress_standard_pair_results.tsv")
write.table(
  standard_pairs,
  out_file,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

params <- model_parameters(fit_standard)
if (!is.null(params$growth_exponents)) {
  write.table(
    params$growth_exponents,
    file.path(out_dir, "growth_parameters.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
}

message("Wrote DStressR standard pair results to: ", out_file)
message("Rows: ", nrow(standard_pairs))
