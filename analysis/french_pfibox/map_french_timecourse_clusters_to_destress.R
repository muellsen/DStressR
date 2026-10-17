#!/usr/bin/env Rscript

# Map the French/PFIbox Mfuzz time-course clusters from Supplementary Table S3
# to DStressR reporter clusters computed from the fitted total-effect matrix.

source(file.path("analysis", "_helpers.R"))
source(file.path("analysis", "_paper_artifacts.R"))

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
out_dir <- Sys.getenv("PFIBOX_ANALYSIS_DIR", analysis_output_dir("french_pfibox", "promoter_activity_analysis"))
run_tag <- Sys.getenv("PFIBOX_RUN_TAG", "reporters2000_perturbations100_T2-18h_treatment_replicate_scalable")
suffix <- Sys.getenv("PFIBOX_CLUSTER_OUTPUT_SUFFIX", "total_effect_ward")
suffix <- if (nzchar(suffix)) paste0("_", suffix) else ""
dstressr_cluster_file <- file.path(out_dir, paste0(
  "pfibox_pa_model_clustered_reporter_clusters_", run_tag, suffix, ".tsv"))

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
  file.path(out_dir, "pfibox_effect_clusters_mapped_to_french_mfuzz.tsv"),
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

write_paper_artifact(file.path(out_dir, paste0("pfibox_paper_mapping_", run_tag, suffix, ".rds")),
  diagnostics = list(long_counts = long_counts),
  annotations = list(mapping = map, reporter_clusters = dstressr),
  provenance = list(run_tag = run_tag, cluster_file = dstressr_cluster_file))

counts <- table(mapped$reporter_cluster, mapped$french_major_cluster)
c2 <- function(x) x * (x - 1) / 2
expected <- sum(c2(rowSums(counts))) * sum(c2(colSums(counts))) / c2(sum(counts))
maximum <- (sum(c2(rowSums(counts))) + sum(c2(colSums(counts)))) / 2
audit <- data.frame(clustered_reporters = nrow(map), matched_reporters = nrow(mapped),
  unmatched_reporters = nrow(map) - nrow(mapped),
  adjusted_rand_index = (sum(c2(counts)) - expected) / (maximum - expected))
write.table(audit, file.path(out_dir, "figure6_cluster_mapping_audit.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE)
