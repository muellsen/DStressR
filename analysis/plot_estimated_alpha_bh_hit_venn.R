#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

out_dir <- file.path(getwd(), "analysis", "outputs", "growth_exponent", "venn")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

standard_eb_file <- file.path(
  getwd(),
  "analysis",
  "outputs",
  "eb_moderated_variance",
  "workflow_vs_destress_eb_promoter_compound_pvalues.tsv"
)
estimated_alpha_file <- file.path(
  getwd(),
  "analysis",
  "outputs",
  "growth_exponent",
  "workflow_vs_destress_eb_estimated_growth_alpha_promoter_compound_pvalues.tsv"
)

for (path in c(standard_eb_file, estimated_alpha_file)) {
  if (!file.exists(path)) {
    stop("Missing required input: ", path, call. = FALSE)
  }
}

safe_has_package <- function(pkg) {
  requireNamespace(pkg, quietly = TRUE)
}

standard <- read.delim(standard_eb_file, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
estimated <- read.delim(estimated_alpha_file, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)

standard$pair_id <- paste(standard$promoter, standard$srn_code, sep = "__")
estimated$pair_id <- paste(estimated$promoter, estimated$srn_code, sep = "__")

membership <- merge(
  standard[, c(
    "promoter",
    "srn_code",
    "pair_id",
    "workflow_padj_by_promoter",
    "destress_eb_padj_by_promoter"
  )],
  estimated[, c("pair_id", "estimated_alpha_eb_padj_by_promoter")],
  by = "pair_id",
  all = FALSE,
  sort = FALSE
)

membership$median_polish <- membership$workflow_padj_by_promoter < 0.05
membership$destress_eb_fixed_alpha <- membership$destress_eb_padj_by_promoter < 0.05
membership$destress_eb_estimated_alpha <- membership$estimated_alpha_eb_padj_by_promoter < 0.05

membership$venn_region <- paste(
  ifelse(membership$median_polish, "Median-polish", ""),
  ifelse(membership$destress_eb_fixed_alpha, "DStressR EB alpha=1", ""),
  ifelse(membership$destress_eb_estimated_alpha, "DStressR EB estimated alpha_g", ""),
  sep = ";"
)
membership$venn_region <- gsub(";+", ";", membership$venn_region)
membership$venn_region <- gsub("^;+|;+$", "", membership$venn_region)
membership$venn_region[membership$venn_region == ""] <- "Not rejected"

write.table(
  membership,
  file.path(out_dir, "bh_rejected_pair_membership_medianpolish_standardeb_estimatedalpha.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

rejected_membership <- membership[membership$venn_region != "Not rejected", ]
write.table(
  rejected_membership,
  file.path(out_dir, "bh_rejected_pair_membership_rejected_only_medianpolish_standardeb_estimatedalpha.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

sets <- list(
  `Median-polish` = membership$pair_id[membership$median_polish],
  `DStressR EB alpha=1` = membership$pair_id[membership$destress_eb_fixed_alpha],
  `DStressR EB estimated alpha_g` = membership$pair_id[membership$destress_eb_estimated_alpha]
)

region_summary <- as.data.frame(table(membership$venn_region), stringsAsFactors = FALSE)
names(region_summary) <- c("region", "n")
region_summary <- region_summary[order(-region_summary$n, region_summary$region), ]
write.table(
  region_summary,
  file.path(out_dir, "bh_rejected_pair_region_counts_medianpolish_standardeb_estimatedalpha.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

if (safe_has_package("ggvenn")) {
  p <- ggvenn::ggvenn(
    sets,
    fill_color = c("#64748b", "#16a34a", "#dc2626"),
    stroke_size = 0.5,
    set_name_size = 3.5,
    text_size = 4
  ) +
    labs(
      title = "BH-rejected promoter-compound pairs",
      subtitle = "Median-polish vs fixed-alpha DStressR EB vs estimated-alpha DStressR EB"
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
  file.path(out_dir, "bh_rejected_pair_venn_medianpolish_standardeb_estimatedalpha.png"),
  p,
  width = 8,
  height = 6.5,
  dpi = 300,
  bg = "white"
)
ggsave(
  file.path(out_dir, "bh_rejected_pair_venn_medianpolish_standardeb_estimatedalpha.pdf"),
  p,
  width = 8,
  height = 6.5,
  bg = "white"
)

print(region_summary)
message("Wrote estimated-alpha Venn outputs to: ", out_dir)
