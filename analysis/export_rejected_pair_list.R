#!/usr/bin/env Rscript

source(file.path("analysis", "_helpers.R"))

methods <- c("median_polish", "destress_standard", "destress_moderated")
out_dir <- comparison_results_dir("hit_overlap")
comparison_file <- comparison_results_dir("pair_level", "pair_level_pvalue_comparison.tsv")
membership_file <- comparison_results_dir("pair_level", "pair_level_hit_membership.tsv")

if (!file.exists(comparison_file) || !file.exists(membership_file)) {
  stop(
    "Missing pair-level comparison outputs.",
    "\nRun analysis/compare_pair_level_pvalues.R first.",
    call. = FALSE
  )
}

comparison <- read_tsv_base(comparison_file)
membership <- read_tsv_base(membership_file)
membership <- membership[membership$hit_region != "not_hit", , drop = FALSE]

out <- merge(
  membership,
  comparison,
  by = c("promoter", "compound", "pair_id"),
  all.x = TRUE,
  sort = FALSE
)

libmap_file <- libmap_path()
if (file.exists(libmap_file)) {
  libmap <- read_tsv_base(libmap_file)
  libmap$libplate <- paste0("lp", libmap[["Library plate"]])
  libmap$compound <- paste(libmap$libplate, libmap[["Well"]], sep = "_")
  libmap$ProductName <- ifelse(
    is.na(libmap$ProductName) | libmap$ProductName == "NA" | libmap$ProductName == "",
    libmap[["Catalog Number"]],
    libmap$ProductName
  )
  keep <- intersect(c("compound", "ProductName", "Catalog Number", "Target"), names(libmap))
  out <- merge(out, libmap[, keep, drop = FALSE], by = "compound", all.x = TRUE, sort = FALSE)
}

out$methods_rejected <- gsub(";", " + ", out$hit_region, fixed = TRUE)
method_cols <- unlist(lapply(methods, function(method) {
  c(method, paste0(method, "_effect"), paste0(method, "_pvalue"), paste0(method, "_padj"))
}))
front_cols <- intersect(
  c("promoter", "compound", "ProductName", "Catalog Number", "Target", "methods_rejected"),
  names(out)
)
out <- out[, c(front_cols, intersect(method_cols, names(out))), drop = FALSE]
out <- out[order(out$promoter, out$compound), , drop = FALSE]

write.table(out, file.path(out_dir, "differential_pair_list.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
utils::write.csv(out, file.path(out_dir, "differential_pair_list.csv"),
                 row.names = FALSE, quote = TRUE)

message("Wrote explicit differential pair list to: ", out_dir)
message("Rows: ", nrow(out))
