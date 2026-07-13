#!/usr/bin/env Rscript

source(file.path("analysis", "_helpers.R"))

root <- analysis_data_root()
out_dir <- file.path(getwd(), "analysis", "outputs", "venn")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

rejected_file <- file.path(out_dir, "bh_rejected_pair_membership_rejected_only.tsv")
libmap_file <- file.path(root, "00-import", "Campylobacter", "LibMap.txt")
effects_file <- file.path(
  getwd(),
  "analysis",
  "outputs",
  "eb_moderated_variance",
  "workflow_vs_destress_eb_promoter_compound_pvalues.tsv"
)

read_tsv_base <- function(path) {
  read.delim(path, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
}

if (!file.exists(rejected_file)) {
  stop("Missing rejected-only table. Run analysis/plot_bh_hit_venn.R first.", call. = FALSE)
}
if (!file.exists(effects_file)) {
  stop("Missing EB promoter-compound table. Run analysis/apply_eb_moderated_variances.R first.", call. = FALSE)
}

rejected <- read_tsv_base(rejected_file)
effects <- read_tsv_base(effects_file)
libmap <- read_tsv_base(libmap_file)

libmap$libplate <- paste0("lp", libmap[["Library plate"]])
libmap$srn_code <- paste(libmap$libplate, libmap[["Well"]], sep = "_")
libmap$ProductName <- ifelse(
  is.na(libmap$ProductName) | libmap$ProductName == "NA" | libmap$ProductName == "",
  libmap[["Catalog Number"]],
  libmap$ProductName
)

effects_key <- effects[, c(
  "promoter",
  "srn_code",
  "log2FC.polished",
  "destress_global_effect",
  "destress_specific_effect",
  "destress_eb_effect_centered",
  "workflow_padj_by_promoter",
  "destress_gaussian_padj_by_promoter",
  "destress_eb_padj_by_promoter"
)]

out <- merge(
  rejected[, c(
    "promoter",
    "srn_code",
    "median_polish",
    "destress_gaussian",
    "destress_eb",
    "venn_region"
  )],
  effects_key,
  by = c("promoter", "srn_code"),
  all.x = TRUE,
  sort = FALSE
)
out <- merge(
  out,
  libmap[, c("srn_code", "ProductName", "Catalog Number", "Target")],
  by = "srn_code",
  all.x = TRUE,
  sort = FALSE
)

out$methods_rejected <- gsub(";", " + ", out$venn_region, fixed = TRUE)
out <- out[, c(
  "promoter",
  "srn_code",
  "ProductName",
  "Catalog Number",
  "Target",
  "methods_rejected",
  "median_polish",
  "destress_gaussian",
  "destress_eb",
  "workflow_padj_by_promoter",
  "destress_gaussian_padj_by_promoter",
  "destress_eb_padj_by_promoter",
  "log2FC.polished",
  "destress_global_effect",
  "destress_specific_effect",
  "destress_eb_effect_centered"
)]

out <- out[order(out$promoter, out$ProductName, out$srn_code), ]

write.table(
  out,
  file.path(out_dir, "bh_rejected_compound_promoter_pairs_explicit.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  out,
  file.path(out_dir, "bh_rejected_compound_promoter_pairs_explicit.csv"),
  row.names = FALSE,
  quote = TRUE
)

message("Wrote explicit rejected pair list to: ", out_dir)
message("Rows: ", nrow(out))
