#!/usr/bin/env Rscript

source(file.path("analysis", "_helpers.R"))

suppressPackageStartupMessages({
  library(ggplot2)
})

out_dir <- file.path(getwd(), "analysis", "outputs", "venn")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

input_file <- file.path(
  getwd(),
  "analysis",
  "outputs",
  "eb_moderated_variance",
  "workflow_vs_destress_eb_promoter_compound_pvalues.tsv"
)

if (!file.exists(input_file)) {
  stop(
    "Missing EB promoter-compound table: ", input_file,
    "\nRun analysis/apply_eb_moderated_variances.R first.",
    call. = FALSE
  )
}

safe_has_package <- function(pkg) {
  requireNamespace(pkg, quietly = TRUE)
}

tab <- read.delim(input_file, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
tab$pair_id <- paste(tab$promoter, tab$srn_code, sep = "__")

membership <- data.frame(
  promoter = tab$promoter,
  srn_code = tab$srn_code,
  pair_id = tab$pair_id,
  median_polish = tab$workflow_padj_by_promoter < 0.05,
  destress_gaussian = tab$destress_gaussian_padj_by_promoter < 0.05,
  destress_eb = tab$destress_eb_padj_by_promoter < 0.05,
  workflow_padj_by_promoter = tab$workflow_padj_by_promoter,
  destress_gaussian_padj_by_promoter = tab$destress_gaussian_padj_by_promoter,
  destress_eb_padj_by_promoter = tab$destress_eb_padj_by_promoter
)
membership$venn_region <- paste(
  ifelse(membership$median_polish, "Median-polish", ""),
  ifelse(membership$destress_gaussian, "DStressR Gaussian", ""),
  ifelse(membership$destress_eb, "DStressR EB", ""),
  sep = ";"
)
membership$venn_region <- gsub(";+", ";", membership$venn_region)
membership$venn_region <- gsub("^;+|;+$", "", membership$venn_region)
membership$venn_region[membership$venn_region == ""] <- "Not rejected"

write.table(
  membership,
  file.path(out_dir, "bh_rejected_pair_membership.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

rejected_membership <- membership[membership$venn_region != "Not rejected", ]
write.table(
  rejected_membership,
  file.path(out_dir, "bh_rejected_pair_membership_rejected_only.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

sets <- list(
  `Median-polish` = membership$pair_id[membership$median_polish],
  `DStressR Gaussian` = membership$pair_id[membership$destress_gaussian],
  `DStressR EB` = membership$pair_id[membership$destress_eb]
)

region_summary <- as.data.frame(table(membership$venn_region), stringsAsFactors = FALSE)
names(region_summary) <- c("region", "n")
region_summary <- region_summary[order(-region_summary$n, region_summary$region), ]
write.table(
  region_summary,
  file.path(out_dir, "bh_rejected_pair_region_counts.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

if (safe_has_package("ggvenn")) {
  p <- ggvenn::ggvenn(
    sets,
    fill_color = c("#64748b", "#2563eb", "#16a34a"),
    stroke_size = 0.5,
    set_name_size = 4,
    text_size = 4
  ) +
    labs(
      title = "BH-rejected promoter-compound pairs",
      subtitle = "BH adjusted p < 0.05 within promoter"
    ) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
} else {
  stop("Package ggvenn is not installed; cannot draw Venn diagram.", call. = FALSE)
}

ggsave(
  file.path(out_dir, "bh_rejected_pair_venn_three_methods.png"),
  p,
  width = 7,
  height = 6,
  dpi = 300,
  bg = "white"
)
ggsave(
  file.path(out_dir, "bh_rejected_pair_venn_three_methods.pdf"),
  p,
  width = 7,
  height = 6,
  bg = "white"
)

print(region_summary)
message("Wrote Venn outputs to: ", out_dir)
