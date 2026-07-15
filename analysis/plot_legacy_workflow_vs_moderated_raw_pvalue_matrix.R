#!/usr/bin/env Rscript

source(file.path("analysis", "_helpers.R"))

suppressPackageStartupMessages({
  library(ggplot2)
})

out_dir <- comparison_results_dir("promoter_pvalue_matrix")

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
    !(workflow$promoter %in% c("PCJnc20", "PCjas704")),
  c("promoter", "srn_code", "pvalue"),
  drop = FALSE
]

legacy_pair <- stats::aggregate(
  pvalue ~ promoter + srn_code,
  workflow,
  max,
  na.rm = TRUE
)
names(legacy_pair) <- c("promoter", "compound", "pvalue")
legacy_pair$method <- "Legacy workflow"

moderated <- read_package_pair_results("destress_moderated")
moderated <- data.frame(
  promoter = moderated$promoter,
  compound = moderated$compound,
  pvalue = moderated$destress_moderated_pvalue,
  method = "DStressR moderated model",
  stringsAsFactors = FALSE
)

common_keys <- intersect(
  paste(legacy_pair$promoter, legacy_pair$compound, sep = "\r"),
  paste(moderated$promoter, moderated$compound, sep = "\r")
)
legacy_pair <- legacy_pair[paste(legacy_pair$promoter, legacy_pair$compound, sep = "\r") %in% common_keys, ]
moderated <- moderated[paste(moderated$promoter, moderated$compound, sep = "\r") %in% common_keys, ]

plot_data <- rbind(legacy_pair, moderated)
plot_data <- plot_data[is.finite(plot_data$pvalue), , drop = FALSE]
plot_data$method <- factor(plot_data$method, levels = c("Legacy workflow", "DStressR moderated model"))
plot_data$promoter <- factor(plot_data$promoter, levels = sort(unique(plot_data$promoter)))

summary <- do.call(rbind, lapply(split(plot_data, list(plot_data$method, plot_data$promoter), drop = TRUE), function(x) {
  data.frame(
    method = as.character(x$method[1]),
    promoter = as.character(x$promoter[1]),
    n = nrow(x),
    p_lt_0.05 = sum(x$pvalue < 0.05, na.rm = TRUE),
    median_pvalue = stats::median(x$pvalue, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))

p <- ggplot(plot_data, aes(pvalue)) +
  geom_histogram(binwidth = 0.05, boundary = 0, fill = "#2563eb", color = "white", linewidth = 0.12) +
  facet_grid(promoter ~ method) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  theme_bw(base_size = 8) +
  theme(
    strip.background = element_rect(fill = "#f8fafc", color = "#cbd5e1"),
    strip.text.x = element_text(face = "bold", size = 9),
    strip.text.y = element_text(face = "bold", size = 7, angle = 0),
    panel.grid.minor = element_blank(),
    panel.spacing.y = unit(0.08, "lines"),
    panel.spacing.x = unit(0.35, "lines")
  ) +
  labs(
    title = "Raw p-value histograms by promoter",
    subtitle = "Legacy workflow median-polish p-values vs DStressR moderated model",
    x = "Raw p-value",
    y = "Promoter-compound pairs"
  )

write.table(
  summary,
  file.path(out_dir, "legacy_workflow_vs_destress_moderated_raw_pvalue_matrix_summary.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

ggsave(
  file.path(out_dir, "legacy_workflow_vs_destress_moderated_raw_pvalue_matrix.png"),
  p,
  width = 8.5,
  height = 18,
  dpi = 300,
  bg = "white"
)
ggsave(
  file.path(out_dir, "legacy_workflow_vs_destress_moderated_raw_pvalue_matrix.pdf"),
  p,
  width = 8.5,
  height = 18,
  bg = "white"
)

message("Wrote legacy workflow vs moderated p-value matrix to: ", out_dir)
message("Common promoter-compound pairs: ", length(common_keys))
print(summary[order(summary$promoter, summary$method), ])
