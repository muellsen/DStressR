#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

out_dir <- file.path(getwd(), "analysis", "outputs", "normalized_matrix", "effect_histograms")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

compound_file <- file.path(
  getwd(),
  "analysis",
  "outputs",
  "normalized_matrix",
  "normalized_promoter_compound_matrix_long.tsv"
)
replicate_file <- file.path(
  getwd(),
  "analysis",
  "outputs",
  "growth_exponent",
  "workflow_vs_destress_eb_estimated_growth_alpha_replicate_pvalues.tsv"
)
libmap_file <- "/Users/cmueller/Documents/GitHub/campylobacter_stressregnet/workflow/data/00-import/Campylobacter/LibMap.txt"

compound <- read.delim(compound_file, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
replicate <- read.delim(replicate_file, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
libmap <- read.delim(libmap_file, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
libmap$libplate <- paste0("lp", libmap[["Library plate"]])
libmap$srn_code <- paste(libmap$libplate, libmap[["Well"]], sep = "_")
dmso_label <- ifelse(
  is.na(libmap$ProductName) | libmap$ProductName == "NA" | !nzchar(libmap$ProductName),
  libmap[["Catalog Number"]],
  libmap$ProductName
)
dmso_srn_codes <- libmap$srn_code[dmso_label %in% c("DMSO", "DMSO noisy")]

compound_df <- data.frame(
  promoter = compound$promoter,
  source = "Compounds",
  effect = as.numeric(compound$destress_eb_effect_centered),
  stringsAsFactors = FALSE
)
dmso_df <- replicate[replicate$srn_code %in% dmso_srn_codes, , drop = FALSE]
if (nrow(dmso_df) == 0) {
  stop("No DMSO rows found in replicate-level output. Check library map labels and srn_code construction.", call. = FALSE)
}
dmso_df <- data.frame(
  promoter = dmso_df$promoter,
  source = "DMSO controls",
  effect = as.numeric(dmso_df$destress_eb_effect_centered),
  stringsAsFactors = FALSE
)
plot_df <- rbind(compound_df, dmso_df)
plot_df <- plot_df[is.finite(plot_df$effect), , drop = FALSE]
plot_df$source <- factor(plot_df$source, levels = c("Compounds", "DMSO controls"))
plot_df$promoter <- factor(plot_df$promoter, levels = sort(unique(plot_df$promoter)))

effect_limit <- stats::quantile(abs(plot_df$effect), 0.995, na.rm = TRUE)
if (!is.finite(effect_limit) || effect_limit <= 0) {
  effect_limit <- max(abs(plot_df$effect), na.rm = TRUE)
}

summary_df <- do.call(
  rbind,
  lapply(split(plot_df, list(plot_df$promoter, plot_df$source), drop = TRUE), function(d) {
    data.frame(
      promoter = as.character(d$promoter[1]),
      source = as.character(d$source[1]),
      n = nrow(d),
      mean = mean(d$effect, na.rm = TRUE),
      median = stats::median(d$effect, na.rm = TRUE),
      sd = stats::sd(d$effect, na.rm = TRUE),
      mad = stats::mad(d$effect, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
)
write.table(
  summary_df,
  file.path(out_dir, "effect_histogram_compounds_vs_dmso_summary.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

make_panel <- function(d, title, ncol = 4) {
  ggplot(d, aes(effect)) +
    geom_histogram(bins = 45, fill = "#4E79A7", color = "white", linewidth = 0.12) +
    geom_vline(xintercept = 0, color = "#303030", linewidth = 0.28) +
    facet_grid(promoter ~ source, scales = "free_y") +
    coord_cartesian(xlim = c(-effect_limit, effect_limit)) +
    theme_light(base_size = 8) +
    theme(
      panel.grid.minor = element_blank(),
      strip.text.y = element_text(angle = 0, size = 7),
      strip.text.x = element_text(size = 8),
      plot.title.position = "plot"
    ) +
    labs(
      title = title,
      subtitle = "Left: adjusted effects over compounds. Right: corresponding DMSO control effects from replicate-level model output.",
      x = "Centered DStressR EB effect",
      y = "Count"
    )
}

p_all <- make_panel(
  plot_df,
  "Compound-effect and DMSO-control distributions per promoter"
)
ggsave(
  file.path(out_dir, "effect_histograms_by_promoter_with_dmso_all.png"),
  p_all,
  width = 10,
  height = 18,
  dpi = 300,
  bg = "white"
)
ggsave(
  file.path(out_dir, "effect_histograms_by_promoter_with_dmso_all.pdf"),
  p_all,
  width = 10,
  height = 18,
  bg = "white"
)

promoters <- levels(plot_df$promoter)
page_id <- ceiling(seq_along(promoters) / 9)
page_map <- split(promoters, page_id)
for (i in seq_along(page_map)) {
  page_df <- plot_df[plot_df$promoter %in% page_map[[i]], , drop = FALSE]
  page_df$promoter <- factor(as.character(page_df$promoter), levels = page_map[[i]])
  p_page <- make_panel(
    page_df,
    paste0("Compound-effect and DMSO-control distributions per promoter, page ", i)
  )
  ggsave(
    file.path(out_dir, sprintf("effect_histograms_by_promoter_with_dmso_page%02d.png", i)),
    p_page,
    width = 10,
    height = 8.5,
    dpi = 300,
    bg = "white"
  )
  ggsave(
    file.path(out_dir, sprintf("effect_histograms_by_promoter_with_dmso_page%02d.pdf", i)),
    p_page,
    width = 10,
    height = 8.5,
    bg = "white"
  )
}

message("Wrote compound-vs-DMSO effect histograms to: ", out_dir)
print(utils::head(summary_df, 12), row.names = FALSE)
