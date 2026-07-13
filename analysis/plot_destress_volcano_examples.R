#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

if (requireNamespace("DStressR", quietly = TRUE)) {
  library(DStressR)
} else {
  devtools::load_all(".", quiet = TRUE)
}

out_dir <- file.path(getwd(), "analysis", "outputs", "volcano")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

estimated_alpha_file <- file.path(
  getwd(),
  "analysis",
  "outputs",
  "growth_exponent",
  "workflow_vs_destress_eb_estimated_growth_alpha_promoter_compound_pvalues.tsv"
)

if (!file.exists(estimated_alpha_file)) {
  stop("Missing input table: ", estimated_alpha_file, call. = FALSE)
}

tab <- read.delim(estimated_alpha_file, sep = "\t", check.names = FALSE)

libmap_file <- "/Users/cmueller/Documents/GitHub/campylobacter_stressregnet/workflow/data/00-import/Campylobacter/LibMap.txt"
if (file.exists(libmap_file)) {
  libmap <- read.delim(libmap_file, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
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

p <- plot_volcano(
  tab,
  effect = "destress_eb_effect_centered",
  padj = "estimated_alpha_eb_padj_by_promoter",
  promoter = "promoter",
  compound = "srn_code",
  compound_label = "compound_label",
  fdr = 0.05,
  lfc = 0,
  top_n = 14,
  top_promoters = 6,
  label_by = "pair",
  title = "DStressR EB volcano plot",
  subtitle = "Estimated growth exponent model; colored promoters are the top hit-rich promoter groups"
)

ggsave(
  file.path(out_dir, "destress_eb_estimated_alpha_volcano.png"),
  p,
  width = 9.5,
  height = 7,
  dpi = 300,
  bg = "white"
)
ggsave(
  file.path(out_dir, "destress_eb_estimated_alpha_volcano.pdf"),
  p,
  width = 9.5,
  height = 7,
  bg = "white"
)

message("Wrote volcano plots to: ", out_dir)
