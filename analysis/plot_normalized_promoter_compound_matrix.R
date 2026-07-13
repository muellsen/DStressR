#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

if (requireNamespace("DStressR", quietly = TRUE)) {
  library(DStressR)
} else {
  devtools::load_all(".", quiet = TRUE)
}

out_dir <- file.path(getwd(), "analysis", "outputs", "normalized_matrix")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

result_file <- file.path(
  getwd(),
  "analysis",
  "outputs",
  "growth_exponent",
  "workflow_vs_destress_eb_estimated_growth_alpha_promoter_compound_pvalues.tsv"
)
libmap_file <- "/Users/cmueller/Documents/GitHub/campylobacter_stressregnet/workflow/data/00-import/Campylobacter/LibMap.txt"

read_tsv_base <- function(path) {
  read.delim(path, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
}

tab <- read_tsv_base(result_file)

if (file.exists(libmap_file)) {
  libmap <- read_tsv_base(libmap_file)
  libmap$libplate <- paste0("lp", libmap[["Library plate"]])
  libmap$srn_code <- paste(libmap$libplate, libmap[["Well"]], sep = "_")
  libmap$compound_label <- libmap$ProductName
  missing_name <- is.na(libmap$compound_label) |
    libmap$compound_label == "NA" |
    !nzchar(libmap$compound_label)
  libmap$compound_label[missing_name] <- libmap[["Catalog Number"]][missing_name]
  missing_name <- is.na(libmap$compound_label) |
    libmap$compound_label == "NA" |
    !nzchar(libmap$compound_label)
  libmap$compound_label[missing_name] <- libmap$srn_code[missing_name]
  tab <- merge(
    tab,
    libmap[, c("srn_code", "compound_label")],
    by = "srn_code",
    all.x = TRUE,
    sort = FALSE
  )
} else {
  tab$compound_label <- tab$srn_code
}

tab$compound_label[is.na(tab$compound_label) | !nzchar(tab$compound_label)] <-
  tab$srn_code[is.na(tab$compound_label) | !nzchar(tab$compound_label)]
tab$compound_display <- paste0(tab$compound_label, " [", tab$srn_code, "]")

matrix_long <- tab[, c(
  "promoter",
  "srn_code",
  "compound_label",
  "compound_display",
  "log2FC",
  "destress_global_effect",
  "destress_specific_effect",
  "destress_eb_effect_centered",
  "destress_eb_pvalue",
  "estimated_alpha_eb_padj_by_promoter"
)]
write.table(
  matrix_long,
  file.path(out_dir, "normalized_promoter_compound_matrix_long.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

p_full <- plot_response_heatmap(
  matrix_long,
  value = "destress_eb_effect_centered",
  promoter = "promoter",
  compound = "srn_code",
  compound_label = "compound_label",
  show_compound_ids = TRUE,
  top_n_compounds = Inf,
  title = "Full normalized promoter-by-compound matrix",
  subtitle = "All compounds; color scale clipped at the 98th percentile of absolute effects"
)
effect_mat <- attr(p_full, "response_matrix")

p_log2fc <- plot_response_heatmap(
  matrix_long,
  value = "log2FC",
  promoter = "promoter",
  compound = "srn_code",
  compound_label = "compound_label",
  show_compound_ids = TRUE,
  top_n_compounds = Inf
)
log2fc_mat <- attr(p_log2fc, "response_matrix")

write.table(
  data.frame(promoter = rownames(effect_mat), effect_mat, check.names = FALSE),
  file.path(out_dir, "normalized_promoter_by_compound_matrix_destress_eb_effect.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
write.table(
  data.frame(promoter = rownames(log2fc_mat), log2fc_mat, check.names = FALSE),
  file.path(out_dir, "normalized_promoter_by_compound_matrix_log2fc.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

p_top <- plot_response_heatmap(
  matrix_long,
  value = "destress_eb_effect_centered",
  promoter = "promoter",
  compound = "srn_code",
  compound_label = "compound_label",
  show_compound_ids = TRUE,
  top_n_compounds = 160,
  title = "Normalized promoter-by-compound matrix",
  subtitle = "Top 160 compounds by mean absolute DStressR EB effect; rows and columns hierarchically clustered"
)

ggsave(
  file.path(out_dir, "normalized_promoter_compound_matrix_top160_heatmap.png"),
  p_top,
  width = 12,
  height = 6.5,
  dpi = 300,
  bg = "white"
)
ggsave(
  file.path(out_dir, "normalized_promoter_compound_matrix_top160_heatmap.pdf"),
  p_top,
  width = 12,
  height = 6.5,
  bg = "white"
)

ggsave(
  file.path(out_dir, "normalized_promoter_compound_matrix_full_heatmap.png"),
  p_full,
  width = 14,
  height = 6.5,
  dpi = 300,
  bg = "white"
)
ggsave(
  file.path(out_dir, "normalized_promoter_compound_matrix_full_heatmap.pdf"),
  p_full,
  width = 14,
  height = 6.5,
  bg = "white"
)

summary_df <- data.frame(
  quantity = c(
    "n_promoters",
    "n_compounds",
    "matrix_entries",
    "effect_min",
    "effect_median",
    "effect_max",
    "effect_abs_98pct"
  ),
  value = c(
    nrow(effect_mat),
    ncol(effect_mat),
    length(effect_mat),
    min(effect_mat, na.rm = TRUE),
    stats::median(effect_mat, na.rm = TRUE),
    max(effect_mat, na.rm = TRUE),
    attr(p_full, "color_limit")
  )
)
write.table(
  summary_df,
  file.path(out_dir, "normalized_promoter_compound_matrix_summary.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

message("Wrote normalized matrices and heatmaps to: ", out_dir)
print(summary_df, row.names = FALSE)
