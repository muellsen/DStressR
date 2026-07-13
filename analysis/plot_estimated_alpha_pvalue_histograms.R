#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

root <- "/Users/cmueller/Documents/GitHub/campylobacter_stressregnet/workflow/data"
out_dir <- file.path(getwd(), "analysis", "outputs", "growth_exponent", "pvalue_histograms")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

standard_eb_file <- file.path(
  getwd(),
  "analysis",
  "outputs",
  "eb_moderated_variance",
  "workflow_vs_destress_eb_replicate_pvalues.tsv"
)
estimated_alpha_file <- file.path(
  getwd(),
  "analysis",
  "outputs",
  "growth_exponent",
  "workflow_vs_destress_eb_estimated_growth_alpha_replicate_pvalues.tsv"
)
libmap_file <- file.path(root, "00-import", "Campylobacter", "LibMap.txt")

for (path in c(standard_eb_file, estimated_alpha_file, libmap_file)) {
  if (!file.exists(path)) {
    stop("Missing required input: ", path, call. = FALSE)
  }
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
dmso_srn_codes <- libmap$srn_code[libmap$ProductName %in% c("DMSO", "DMSO noisy")]

standard <- read_tsv_base(standard_eb_file)
estimated <- read_tsv_base(estimated_alpha_file)

key_cols <- c("promoter_libplate_replicate", "promoter", "libplate", "replicate", "srn_code")
comparison <- merge(
  standard[, c(key_cols, "pvalue", "destress_eb_pvalue")],
  estimated[, c(key_cols, "destress_eb_pvalue", "growth_alpha")],
  by = key_cols,
  all = FALSE,
  suffixes = c("_standard_eb", "_estimated_alpha_eb"),
  sort = FALSE
)

comparison <- comparison[
  is.finite(comparison$pvalue) &
    is.finite(comparison$destress_eb_pvalue_standard_eb) &
    is.finite(comparison$destress_eb_pvalue_estimated_alpha_eb) &
    !(comparison$srn_code %in% dmso_srn_codes),
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
    pvalue = comparison$destress_eb_pvalue_standard_eb,
    method = "DStressR EB (alpha = 1)",
    stringsAsFactors = FALSE
  ),
  data.frame(
    promoter_replicate = comparison$promoter_replicate,
    pvalue = comparison$destress_eb_pvalue_estimated_alpha_eb,
    method = "DStressR EB (estimated alpha_g)",
    stringsAsFactors = FALSE
  )
)
long_df$method <- factor(
  long_df$method,
  levels = c("Median-polish", "DStressR EB (alpha = 1)", "DStressR EB (estimated alpha_g)")
)

panel_summary <- do.call(
  rbind,
  lapply(split(long_df, list(long_df$promoter_replicate, long_df$method), drop = TRUE), function(d) {
    data.frame(
      promoter_replicate = as.character(d$promoter_replicate[1]),
      method = as.character(d$method[1]),
      n = nrow(d),
      p_lt_0.05 = mean(d$pvalue < 0.05, na.rm = TRUE),
      median_p = stats::median(d$pvalue, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
)
write.table(
  panel_summary,
  file.path(out_dir, "all_promoter_replicate_pvalue_histogram_summary_estimated_alpha.tsv"),
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
    subtitle = "Columns: original median-polish workflow, fixed-alpha DStressR EB, estimated-alpha DStressR EB",
    x = "p-value",
    y = "Count"
  )

height <- max(14, length(panel_order) * 0.38)

ggsave(
  file.path(out_dir, "all_promoter_replicate_pvalue_histograms_medianpolish_standardeb_estimatedalpha.png"),
  p,
  width = 12,
  height = height,
  dpi = 300,
  limitsize = FALSE
)
ggsave(
  file.path(out_dir, "all_promoter_replicate_pvalue_histograms_medianpolish_standardeb_estimatedalpha.pdf"),
  p,
  width = 12,
  height = height,
  limitsize = FALSE
)

method_summary <- aggregate(
  pvalue ~ method,
  long_df,
  function(x) c(n = length(x), p_lt_0.05 = mean(x < 0.05), median_p = stats::median(x))
)
method_summary <- do.call(data.frame, method_summary)
names(method_summary) <- c("method", "n", "p_lt_0.05", "median_p")
print(method_summary)
message("Wrote estimated-alpha p-value histogram panel to: ", out_dir)
message("Panel rows: ", length(panel_order))
