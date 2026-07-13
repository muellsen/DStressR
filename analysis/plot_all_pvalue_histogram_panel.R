#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

root <- "/Users/cmueller/Documents/GitHub/campylobacter_stressregnet/workflow/data"
out_dir <- file.path(getwd(), "analysis", "outputs", "pvalue_histograms")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

eb_file <- file.path(
  getwd(),
  "analysis",
  "outputs",
  "eb_moderated_variance",
  "workflow_vs_destress_eb_replicate_pvalues.tsv"
)
libmap_file <- file.path(root, "00-import", "Campylobacter", "LibMap.txt")

if (!file.exists(eb_file)) {
  stop(
    "Missing EB comparison table: ", eb_file,
    "\nRun analysis/apply_eb_moderated_variances.R first.",
    call. = FALSE
  )
}

read_tsv_base <- function(path) {
  read.delim(path, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
}

libmap <- read_tsv_base(libmap_file)
libmap$libplate <- paste0("lp", libmap[["Library plate"]])
libmap$srn_code <- paste(libmap$libplate, libmap[["Well"]], sep = "_")
libmap$ProductName <- ifelse(
  is.na(libmap$ProductName) | libmap$ProductName == "NA" | libmap$ProductName == "",
  libmap[["Catalog Number"]],
  libmap$ProductName
)
dmso_srn_codes <- libmap$srn_code[libmap$ProductName == "DMSO"]
dmso_noisy_srn_codes <- libmap$srn_code[libmap$ProductName == "DMSO noisy"]

comparison <- read_tsv_base(eb_file)
comparison <- comparison[
  is.finite(comparison$pvalue) &
    is.finite(comparison$destress_pvalue) &
    is.finite(comparison$destress_eb_pvalue) &
    !(comparison$srn_code %in% c(dmso_srn_codes, dmso_noisy_srn_codes)),
]

comparison$promoter_replicate <- paste(comparison$promoter, comparison$replicate, sep = "_")

panel_order <- unique(comparison[order(comparison$promoter, comparison$replicate), "promoter_replicate"])
comparison$promoter_replicate <- factor(comparison$promoter_replicate, levels = panel_order)

long_df <- rbind(
  data.frame(
    promoter_replicate = comparison$promoter_replicate,
    pvalue = comparison$pvalue,
    method = "Median-polish",
    stringsAsFactors = FALSE
  ),
  data.frame(
    promoter_replicate = comparison$promoter_replicate,
    pvalue = comparison$destress_pvalue,
    method = "DStressR Gaussian",
    stringsAsFactors = FALSE
  ),
  data.frame(
    promoter_replicate = comparison$promoter_replicate,
    pvalue = comparison$destress_eb_pvalue,
    method = "DStressR EB",
    stringsAsFactors = FALSE
  )
)
long_df$method <- factor(
  long_df$method,
  levels = c("Median-polish", "DStressR Gaussian", "DStressR EB")
)

panel_summary <- do.call(
  rbind,
  lapply(split(long_df, list(long_df$promoter_replicate, long_df$method), drop = TRUE), function(d) {
    data.frame(
      promoter_replicate = as.character(d$promoter_replicate[1]),
      method = as.character(d$method[1]),
      n = nrow(d),
      p_lt_0.05 = mean(d$pvalue < 0.05, na.rm = TRUE),
      median_p = stats::median(d$pvalue, na.rm = TRUE)
    )
  })
)
write.table(
  panel_summary,
  file.path(out_dir, "all_promoter_replicate_pvalue_histogram_summary.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

p <- ggplot(long_df, aes(pvalue)) +
  geom_histogram(
    breaks = seq(0, 1, by = 0.025),
    fill = "#334155",
    color = "white",
    linewidth = 0.08
  ) +
  geom_vline(xintercept = 0.05, color = "#b91c1c", linewidth = 0.25) +
  facet_grid(promoter_replicate ~ method, scales = "free_y") +
  coord_cartesian(xlim = c(0, 1)) +
  theme_bw(base_size = 7) +
  theme(
    strip.text.y = element_text(angle = 0, hjust = 0, size = 6),
    strip.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 5),
    axis.text.x = element_text(size = 5),
    panel.grid.minor = element_blank(),
    panel.spacing.y = unit(0.05, "lines")
  ) +
  labs(
    title = "All promoter-replicate p-value histograms",
    subtitle = "Columns: original median-polish workflow, DStressR Gaussian, DStressR EB moderated promoter variance",
    x = "p-value",
    y = "Count"
  )

height <- max(14, length(panel_order) * 0.38)

ggsave(
  file.path(out_dir, "all_promoter_replicate_pvalue_histograms_three_methods.png"),
  p,
  width = 11,
  height = height,
  dpi = 300,
  limitsize = FALSE
)
ggsave(
  file.path(out_dir, "all_promoter_replicate_pvalue_histograms_three_methods.pdf"),
  p,
  width = 11,
  height = height,
  limitsize = FALSE
)

message("Wrote all p-value histogram panel to: ", out_dir)
message("Panel rows: ", length(panel_order))
