#!/usr/bin/env Rscript

source(file.path("analysis", "_helpers.R"))

suppressPackageStartupMessages({
  library(ggplot2)
})

out_dir <- comparison_results_dir("legacy_workflow_all_pvalues")

workflow_pvalue_path <- file.path(
  analysis_data_root(),
  "03-hit_determination",
  "expression_df.pvalues.tsv.gz"
)

if (!file.exists(workflow_pvalue_path)) {
  stop("Missing legacy workflow p-value export: ", workflow_pvalue_path, call. = FALSE)
}

libmap <- read_tsv_base(libmap_path())
libmap$libplate <- paste0("lp", libmap[["Library plate"]])
libmap$compound <- paste(libmap$libplate, libmap[["Well"]], sep = "_")
libmap$ProductName <- ifelse(
  is.na(libmap$ProductName) | libmap$ProductName == "NA" | libmap$ProductName == "",
  libmap[["Catalog Number"]],
  libmap$ProductName
)
dmso_srn_codes <- libmap$compound[libmap$ProductName == "DMSO"]
dmso_noisy_srn_codes <- libmap$compound[libmap$ProductName == "DMSO noisy"]

workflow <- read_tsv_base(workflow_pvalue_path)
workflow <- workflow[
  !(workflow$srn_code %in% c(dmso_srn_codes, dmso_noisy_srn_codes)) &
    !(workflow$promoter %in% c("PCJnc20", "PCjas704")) &
    is.finite(workflow$pvalue),
  c("promoter", "libplate", "replicate", "srn_code", "pvalue"),
  drop = FALSE
]
workflow$promoter <- factor(workflow$promoter, levels = sort(unique(workflow$promoter)))

summary <- do.call(rbind, lapply(split(workflow, workflow$promoter), function(x) {
  data.frame(
    promoter = as.character(x$promoter[1]),
    n = nrow(x),
    p_lt_0.05 = sum(x$pvalue < 0.05, na.rm = TRUE),
    expected_p_lt_0.05 = 0.05 * nrow(x),
    median_pvalue = stats::median(x$pvalue, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))

p_by_promoter <- ggplot(workflow, aes(pvalue)) +
  geom_histogram(binwidth = 0.05, boundary = 0, fill = "#2563eb", color = "white", linewidth = 0.12) +
  facet_wrap(~ promoter, ncol = 6) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  theme_bw(base_size = 9) +
  theme(
    strip.background = element_rect(fill = "#f8fafc", color = "#cbd5e1"),
    strip.text = element_text(face = "bold", size = 8),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "Legacy workflow raw p-values by promoter",
    subtitle = "All replicate/libplate p-values; no max-p promoter-compound collapse",
    x = "Raw p-value",
    y = "Legacy workflow rows"
  )

p_all <- ggplot(workflow, aes(pvalue)) +
  geom_histogram(binwidth = 0.025, boundary = 0, fill = "#2563eb", color = "white", linewidth = 0.12) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank()) +
  labs(
    title = "Legacy workflow raw p-values",
    subtitle = "All replicate/libplate p-values pooled across promoters",
    x = "Raw p-value",
    y = "Legacy workflow rows"
  )

write.table(
  summary,
  file.path(out_dir, "legacy_workflow_all_raw_pvalues_by_promoter_summary.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

ggsave(
  file.path(out_dir, "legacy_workflow_all_raw_pvalues_by_promoter.png"),
  p_by_promoter,
  width = 14,
  height = 10,
  dpi = 300,
  bg = "white"
)
ggsave(
  file.path(out_dir, "legacy_workflow_all_raw_pvalues_by_promoter.pdf"),
  p_by_promoter,
  width = 14,
  height = 10,
  bg = "white"
)
ggsave(
  file.path(out_dir, "legacy_workflow_all_raw_pvalues_pooled.png"),
  p_all,
  width = 6.5,
  height = 4.5,
  dpi = 300,
  bg = "white"
)
ggsave(
  file.path(out_dir, "legacy_workflow_all_raw_pvalues_pooled.pdf"),
  p_all,
  width = 6.5,
  height = 4.5,
  bg = "white"
)

message("Wrote legacy workflow all-p-value plots to: ", out_dir)
message("Rows: ", nrow(workflow))
print(summary[order(-summary$p_lt_0.05), ])
