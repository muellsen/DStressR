#!/usr/bin/env Rscript

source(file.path("analysis", "_helpers.R"))

methods <- c("median_polish", "destress_standard", "destress_moderated")
out_dir <- comparison_results_dir("hit_overlap")
membership_file <- comparison_results_dir("pair_level", "pair_level_hit_membership.tsv")

if (!file.exists(membership_file)) {
  stop(
    "Missing pair-level hit membership table: ", membership_file,
    "\nRun analysis/compare_pair_level_pvalues.R first.",
    call. = FALSE
  )
}
if (!requireNamespace("ggvenn", quietly = TRUE)) {
  stop("Package ggvenn is required for this plot.", call. = FALSE)
}

membership <- read_tsv_base(membership_file)
for (method in methods) {
  membership[[method]] <- as.logical(membership[[method]])
}

sets <- stats::setNames(
  lapply(methods, function(method) membership$pair_id[membership[[method]]]),
  vapply(methods, method_label, character(1))
)

region_counts <- as.data.frame(table(membership$hit_region), stringsAsFactors = FALSE)
names(region_counts) <- c("hit_region", "n")
region_counts <- region_counts[order(-region_counts$n, region_counts$hit_region), ]
write.table(region_counts, file.path(out_dir, "hit_overlap_region_counts.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

p <- ggvenn::ggvenn(
  sets,
  fill_color = c("#64748b", "#2563eb", "#16a34a"),
  stroke_size = 0.5,
  set_name_size = 4,
  text_size = 4
) +
  ggplot2::labs(
    title = "Differential promoter-compound pairs",
    subtitle = "Package outputs; adjusted p < 0.05"
  ) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    plot.background = ggplot2::element_rect(fill = "white", color = NA),
    panel.background = ggplot2::element_rect(fill = "white", color = NA)
  )

ggplot2::ggsave(file.path(out_dir, "hit_overlap_venn_three_methods.png"),
                p, width = 7, height = 6, dpi = 300, bg = "white")
ggplot2::ggsave(file.path(out_dir, "hit_overlap_venn_three_methods.pdf"),
                p, width = 7, height = 6, bg = "white")

print(region_counts)
message("Wrote hit-overlap outputs to: ", out_dir)
