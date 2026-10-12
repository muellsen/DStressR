#!/usr/bin/env Rscript

# Map the French/PFIbox Mfuzz time-course clusters from Supplementary Table S3
# to DStressR reporter clusters computed from the fitted total-effect matrix.

source(file.path("analysis", "_helpers.R"))

if (!requireNamespace("readxl", quietly = TRUE)) {
  stop("Package `readxl` is required to read the French Table S3 Excel file.", call. = FALSE)
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package `ggplot2` is required for cluster-mapping plots.", call. = FALSE)
}
if (!requireNamespace("scales", quietly = TRUE)) {
  stop("Package `scales` is required for cluster-mapping plots.", call. = FALSE)
}

ggplot2 <- asNamespace("ggplot2")
scales <- asNamespace("scales")

source_file <- analysis_path(
  "data-raw", "french_pfibox", "source",
  "french_pfibox_table_s3_clusterdata_mmc2.xlsx"
)
out_dir <- analysis_output_dir("french_pfibox", "promoter_activity_subset")
dstressr_cluster_file <- file.path(
  out_dir,
  "pfibox_pa_model_clustered_reporter_clusters_reporters2000_perturbations100_T2-18h_treatment_replicate_scalable_all_reporters_ward.tsv"
)

if (!file.exists(source_file)) {
  stop("Missing French Table S3 Excel file: ", source_file, call. = FALSE)
}
if (!file.exists(dstressr_cluster_file)) {
  stop("Missing DStressR reporter-cluster file: ", dstressr_cluster_file, call. = FALSE)
}

raw <- readxl::read_excel(source_file, sheet = "ClusterData", col_names = FALSE)
if (nrow(raw) < 4 || ncol(raw) < 14) {
  stop("Unexpected French Table S3 layout.", call. = FALSE)
}

cluster_tab <- raw[-seq_len(2), seq_len(14)]
names(cluster_tab) <- c(
  "promoter",
  "french_major_cluster",
  paste0("french_membership_", seq_len(12))
)
cluster_tab <- as.data.frame(cluster_tab, stringsAsFactors = FALSE)
cluster_tab <- cluster_tab[!is.na(cluster_tab$promoter) & nzchar(cluster_tab$promoter), , drop = FALSE]
cluster_tab$promoter <- as.character(cluster_tab$promoter)
cluster_tab$french_major_cluster <- paste0("F", as.integer(cluster_tab$french_major_cluster))
for (col in paste0("french_membership_", seq_len(12))) {
  cluster_tab[[col]] <- as.numeric(cluster_tab[[col]])
}
cluster_tab$french_max_membership <- do.call(pmax, c(cluster_tab[paste0("french_membership_", seq_len(12))], na.rm = TRUE))

utils::write.table(
  cluster_tab,
  file.path(out_dir, "french_pfibox_table_s3_mfuzz_clusters_clean.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

dstressr <- read_tsv_base(dstressr_cluster_file)
required <- c("reporter", "reporter_cluster", "max_abs_total_effect")
missing <- setdiff(required, names(dstressr))
if (length(missing)) {
  stop("DStressR cluster file is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
}

map <- merge(
  dstressr[, required, drop = FALSE],
  cluster_tab,
  by.x = "reporter",
  by.y = "promoter",
  all.x = TRUE,
  sort = FALSE
)

utils::write.table(
  map,
  file.path(out_dir, "pfibox_destress_ward_clusters_mapped_to_french_mfuzz.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

mapped <- map[!is.na(map$french_major_cluster), , drop = FALSE]
contingency <- as.data.frame.matrix(table(mapped$reporter_cluster, mapped$french_major_cluster))
contingency$dstressr_reporter_cluster <- rownames(contingency)
contingency <- contingency[, c("dstressr_reporter_cluster", setdiff(names(contingency), "dstressr_reporter_cluster")), drop = FALSE]
utils::write.table(
  contingency,
  file.path(out_dir, "pfibox_destress_vs_french_mfuzz_cluster_contingency.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

long_counts <- as.data.frame(table(mapped$reporter_cluster, mapped$french_major_cluster), stringsAsFactors = FALSE)
names(long_counts) <- c("dstressr_reporter_cluster", "french_major_cluster", "n")
row_totals <- stats::aggregate(n ~ dstressr_reporter_cluster, long_counts, sum)
names(row_totals)[2] <- "dstressr_cluster_n"
long_counts <- merge(long_counts, row_totals, by = "dstressr_reporter_cluster", all.x = TRUE, sort = FALSE)
long_counts$proportion <- ifelse(long_counts$dstressr_cluster_n > 0, long_counts$n / long_counts$dstressr_cluster_n, NA_real_)

cluster_sizes <- stats::aggregate(
  reporter ~ reporter_cluster,
  dstressr,
  length
)
names(cluster_sizes) <- c("dstressr_reporter_cluster", "n_destressr_reporters")
french_sizes <- stats::aggregate(
  reporter ~ french_major_cluster,
  mapped,
  length
)
names(french_sizes) <- c("french_major_cluster", "n_mapped_reporters")

utils::write.table(
  long_counts[order(long_counts$dstressr_reporter_cluster, long_counts$french_major_cluster), ],
  file.path(out_dir, "pfibox_destress_vs_french_mfuzz_cluster_counts_long.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

p <- ggplot2$ggplot(
  long_counts,
  ggplot2$aes(french_major_cluster, dstressr_reporter_cluster, fill = proportion)
) +
  ggplot2$geom_tile(color = "white", linewidth = 0.35) +
  ggplot2$geom_text(
    ggplot2$aes(label = ifelse(n > 0, n, "")),
    size = 2.6,
    color = "#111827"
  ) +
  ggplot2$scale_fill_gradient(
    low = "white",
    high = "#3f6fb5",
    labels = scales$percent_format(accuracy = 1),
    name = "Within-cluster\nfraction"
  ) +
  ggplot2$theme_light(base_size = 10) +
  ggplot2$theme(
    panel.grid = ggplot2$element_blank(),
    axis.text = ggplot2$element_text(color = "#4B5563"),
    axis.title = ggplot2$element_text(color = "#111827"),
    legend.position = "right",
    legend.title = ggplot2$element_text(size = 8),
    legend.text = ggplot2$element_text(size = 7)
  ) +
  ggplot2$labs(
    x = "French Mfuzz time-course cluster",
    y = "DStressR total-effect reporter cluster"
  )

base <- file.path(out_dir, "pfibox_destress_vs_french_mfuzz_cluster_overlap")
ggplot2$ggsave(paste0(base, ".png"), p, width = 7.2, height = 4.8, dpi = 300, bg = "white")
ggplot2$ggsave(paste0(base, ".pdf"), p, width = 7.2, height = 4.8, bg = "white")

message("French Mfuzz clusters: ", nrow(cluster_tab), " promoters")
message("DStressR reporter clusters: ", nrow(dstressr), " reporters")
message("Mapped reporters: ", nrow(mapped))
message("Wrote overlap plot: ", paste0(base, ".png"))
