#!/usr/bin/env Rscript

source(file.path("analysis", "_helpers.R"))

out_dir <- comparison_results_dir("default_median_pair_results")
adjustment <- comparison_adjustment()
methods <- c("destress_moderated", "median_polish")

pair_table <- merge_package_pair_results(methods)
pair_table <- add_hit_columns(pair_table, methods, fdr = 0.05, adjustment = adjustment)

libmap <- read_tsv_base(libmap_path())
libmap$library_plate <- paste0("lp", libmap[["Library plate"]])
libmap$compound <- paste(libmap$library_plate, libmap[["Well"]], sep = "_")
libmap$product_name <- ifelse(
  is.na(libmap$ProductName) | libmap$ProductName == "NA" | libmap$ProductName == "",
  libmap[["Catalog Number"]],
  libmap$ProductName
)
libmap$catalog_number <- libmap[["Catalog Number"]]
libmap$target <- libmap$Target
libmap <- libmap[, c("compound", "library_plate", "Well", "catalog_number", "product_name", "target"), drop = FALSE]
names(libmap)[names(libmap) == "Well"] <- "well"

out <- merge(pair_table, libmap, by = "compound", all.x = TRUE, sort = FALSE)

out$default_destress_hit_global <- out$destress_moderated_padj_global < 0.05
out$default_destress_hit_by_promoter <- out$destress_moderated_padj_by_promoter < 0.05
out$median_polish_hit_global <- out$median_polish_padj_global < 0.05
out$median_polish_hit_by_promoter <- out$median_polish_padj_by_promoter < 0.05
out$hit_pairing <- "not_significant"
out$hit_pairing[out$default_destress_hit_global & !out$median_polish_hit_global] <- "default_destress_only"
out$hit_pairing[!out$default_destress_hit_global & out$median_polish_hit_global] <- "median_polish_only"
out$hit_pairing[out$default_destress_hit_global & out$median_polish_hit_global] <- "both"

clean <- out[, c(
  "pair_id",
  "promoter",
  "compound",
  "library_plate",
  "well",
  "catalog_number",
  "product_name",
  "target",
  "hit_pairing",
  "default_destress_hit_global",
  "median_polish_hit_global",
  "default_destress_hit_by_promoter",
  "median_polish_hit_by_promoter",
  "destress_moderated_effect",
  "destress_moderated_pvalue",
  "destress_moderated_padj_global",
  "destress_moderated_padj_by_promoter",
  "median_polish_effect",
  "median_polish_pvalue",
  "median_polish_padj_global",
  "median_polish_padj_by_promoter"
), drop = FALSE]

names(clean) <- c(
  "pair_id",
  "promoter",
  "compound_id",
  "library_plate",
  "well",
  "catalog_number",
  "compound_name",
  "target",
  "hit_pairing_global_fdr_0_05",
  "default_destress_hit_global_fdr_0_05",
  "median_polish_hit_global_fdr_0_05",
  "default_destress_hit_promoter_fdr_0_05",
  "median_polish_hit_promoter_fdr_0_05",
  "default_destress_effect",
  "default_destress_pvalue",
  "default_destress_padj_global",
  "default_destress_padj_by_promoter",
  "median_polish_effect",
  "median_polish_pvalue",
  "median_polish_padj_global",
  "median_polish_padj_by_promoter"
)

clean <- clean[order(
  clean$hit_pairing_global_fdr_0_05 == "not_significant",
  clean$promoter,
  clean$compound_id
), , drop = FALSE]

out_file <- file.path(out_dir, "default_destress_median_polish_pair_results.tsv")
write.table(clean, out_file, sep = "\t", row.names = FALSE, quote = FALSE)

significant <- clean[
  clean$default_destress_hit_global_fdr_0_05 | clean$median_polish_hit_global_fdr_0_05,
  ,
  drop = FALSE
]
significant <- significant[order(
  significant$hit_pairing_global_fdr_0_05,
  significant$promoter,
  significant$compound_id
), , drop = FALSE]

significant_file <- file.path(out_dir, "default_destress_median_polish_significant_union.tsv")
write.table(significant, significant_file, sep = "\t", row.names = FALSE, quote = FALSE)

summary <- as.data.frame(table(clean$hit_pairing_global_fdr_0_05), stringsAsFactors = FALSE)
names(summary) <- c("hit_pairing_global_fdr_0_05", "n_pairs")
summary <- summary[order(summary$hit_pairing_global_fdr_0_05), , drop = FALSE]
write.table(
  summary,
  file.path(out_dir, "default_destress_median_polish_pairing_summary.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

message("Wrote clean default DStressR vs median-polish pair results to: ", out_file)
message("Rows: ", nrow(clean))
message("Wrote significant union to: ", significant_file)
message("Significant union rows: ", nrow(significant))
print(summary)
