source(file.path("analysis", "_helpers.R"))
source(file.path("analysis", "_paper_artifacts.R"))
load_destress_package()
out_dir <- Sys.getenv("PFIBOX_ANALYSIS_DIR", analysis_output_dir("french_pfibox", "promoter_activity_analysis"))
run_tag <- Sys.getenv("PFIBOX_RUN_TAG", "reporters2000_perturbations100_T2-18h_treatment_replicate_scalable")
result_file <- file.path(out_dir, paste0("pfibox_pa_destress_pair_results_", run_tag, ".tsv"))
metadata_file <- file.path(out_dir, paste0("pfibox_pa_run_metadata_", run_tag, ".tsv"))
res <- read_tsv_base(result_file)
metadata <- read_tsv_base(metadata_file)
top_n_reporters <- as.integer(Sys.getenv("PFIBOX_CLUSTER_TOP_REPORTERS", unset = "2000"))
n_reporter_clusters <- as.integer(Sys.getenv("PFIBOX_REPORTER_CLUSTERS", unset = "8"))
n_perturbation_clusters <- as.integer(Sys.getenv("PFIBOX_PERTURBATION_CLUSTERS", unset = "6"))
fdr <- as.numeric(Sys.getenv("PFIBOX_FDR", unset = "0.05"))
cluster_linkage <- Sys.getenv("PFIBOX_CLUSTER_LINKAGE", unset = "ward.D2")
if (!cluster_linkage %in% c("complete", "ward.D2", "average", "single")) {
  stop("PFIBOX_CLUSTER_LINKAGE must be one of: complete, ward.D2, average, single.", call. = FALSE)
}
cluster_distance <- Sys.getenv("PFIBOX_CLUSTER_DISTANCE", unset = "euclidean")
if (!cluster_distance %in% c("euclidean", "correlation", "cosine")) {
  stop("PFIBOX_CLUSTER_DISTANCE must be one of: euclidean, correlation, cosine.", call. = FALSE)
}
cluster_output_suffix <- Sys.getenv("PFIBOX_CLUSTER_OUTPUT_SUFFIX", unset = "total_effect_ward")
if (nzchar(cluster_output_suffix)) {
  cluster_output_suffix <- paste0("_", gsub("[^A-Za-z0-9_-]+", "_", cluster_output_suffix))
}

make_effect_matrix <- function(tab, value) {
  d <- tab[, c("reporter", "perturbation", value), drop = FALSE]
  names(d) <- c("reporter", "perturbation", "effect")
  d$effect <- as.numeric(d$effect)
  d <- d[is.finite(d$effect), , drop = FALSE]
  d <- stats::aggregate(effect ~ reporter + perturbation, d, mean, na.rm = TRUE)
  reporters <- sort(unique(d$reporter))
  perturbations <- sort(unique(d$perturbation))
  mat <- matrix(
    NA_real_,
    nrow = length(reporters),
    ncol = length(perturbations),
    dimnames = list(reporters, perturbations)
  )
  idx <- cbind(match(d$reporter, reporters), match(d$perturbation, perturbations))
  mat[idx] <- d$effect
  mat
}

impute_by_column_mean <- function(mat) {
  out <- mat
  col_means <- colMeans(out, na.rm = TRUE)
  global_mean <- mean(out, na.rm = TRUE)
  col_means[!is.finite(col_means)] <- global_mean
  if (!is.finite(global_mean)) {
    global_mean <- 0
  }
  for (j in seq_len(ncol(out))) {
    miss <- !is.finite(out[, j])
    if (any(miss)) {
      out[miss, j] <- col_means[j]
    }
  }
  out[!is.finite(out)] <- global_mean
  out
}

matrix_distance <- function(x, method = c("euclidean", "correlation", "cosine")) {
  method <- match.arg(method)
  if (method == "euclidean") {
    return(stats::dist(x))
  }
  if (method == "correlation") {
    cmat <- stats::cor(t(x), use = "pairwise.complete.obs")
    cmat[!is.finite(cmat)] <- 0
    diag(cmat) <- 1
    d <- 1 - cmat
    d[d < 0] <- 0
    return(stats::as.dist(d))
  }
  norm <- sqrt(rowSums(x^2))
  sim <- x %*% t(x)
  denom <- tcrossprod(norm)
  sim <- sim / denom
  sim[!is.finite(sim)] <- 0
  diag(sim) <- 1
  d <- 1 - sim
  d[d < 0] <- 0
  stats::as.dist(d)
}

mat_all <- make_effect_matrix(res, "total_effect")
available <- rowSums(is.finite(mat_all))
mat_all <- mat_all[available >= 5, , drop = FALSE]
row_score <- apply(abs(mat_all), 1, max, na.rm = TRUE)
top_n_reporters <- min(top_n_reporters, nrow(mat_all))
selected_reporters <- names(sort(row_score, decreasing = TRUE))[seq_len(top_n_reporters)]
mat <- mat_all[selected_reporters, , drop = FALSE]

cluster_mat <- impute_by_column_mean(mat)
if (cluster_distance == "correlation") {
  variable <- apply(cluster_mat, 1, stats::sd) > 0
  if (any(!variable)) {
    utils::write.table(data.frame(reporter = rownames(mat)[!variable],
      reason = "Constant profile; correlation distance undefined"),
      file.path(out_dir, "correlation_clustering_excluded_profiles.tsv"),
      sep = "\t", quote = FALSE, row.names = FALSE)
    mat <- mat[variable, , drop = FALSE]
    cluster_mat <- cluster_mat[variable, , drop = FALSE]
  }
}
row_hc <- stats::hclust(matrix_distance(cluster_mat, cluster_distance), method = cluster_linkage)
col_hc <- stats::hclust(matrix_distance(t(cluster_mat), cluster_distance), method = cluster_linkage)
row_order <- rownames(mat)[row_hc$order]
col_order <- colnames(mat)[col_hc$order]
row_cluster <- stats::cutree(row_hc, k = min(n_reporter_clusters, nrow(mat)))[rownames(mat)]
col_cluster <- stats::cutree(col_hc, k = min(n_perturbation_clusters, ncol(mat)))[colnames(mat)]

ordered_mat <- mat[row_order, col_order, drop = FALSE]
plot_long <- as.data.frame(as.table(ordered_mat), stringsAsFactors = FALSE)
names(plot_long) <- c("reporter", "perturbation", "effect")
plot_long$reporter <- factor(plot_long$reporter, levels = rev(row_order))
plot_long$perturbation <- factor(plot_long$perturbation, levels = col_order)
block_rows <- data.frame(reporter = rownames(mat),
                        reporter_cluster = paste0("R", row_cluster[rownames(mat)]))
block_cols <- data.frame(perturbation = colnames(mat),
                        perturbation_cluster = paste0("C", col_cluster[colnames(mat)]))
block_long <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
names(block_long) <- c("reporter", "perturbation", "effect")
block_long <- merge(block_long, block_rows, by = "reporter", sort = FALSE)
block_long <- merge(block_long, block_cols, by = "perturbation", sort = FALSE)
block_summary <- stats::aggregate(
  effect ~ reporter_cluster + perturbation_cluster,
  block_long,
  mean,
  na.rm = TRUE
)
block_summary$reporter_cluster <- factor(
  block_summary$reporter_cluster,
  levels = paste0("R", sort(unique(row_cluster)))
)
block_summary$perturbation_cluster <- factor(
  block_summary$perturbation_cluster,
  levels = paste0("C", sort(unique(col_cluster)))
)

reporter_clusters <- data.frame(
  reporter = names(row_cluster),
  reporter_cluster = paste0("R", row_cluster),
  row_order = match(names(row_cluster), row_order),
  max_abs_total_effect = row_score[names(row_cluster)],
  stringsAsFactors = FALSE
)
perturbation_clusters <- data.frame(
  perturbation = names(col_cluster),
  perturbation_cluster = paste0("C", col_cluster),
  column_order = match(names(col_cluster), col_order),
  mean_abs_total_effect = colMeans(abs(mat), na.rm = TRUE)[names(col_cluster)],
  stringsAsFactors = FALSE
)
utils::write.table(
  reporter_clusters[order(reporter_clusters$reporter_cluster, reporter_clusters$row_order), ],
  file.path(out_dir, paste0("pfibox_pa_model_clustered_reporter_clusters_", run_tag, cluster_output_suffix, ".tsv")),
  sep = "\t", quote = FALSE, row.names = FALSE
)
utils::write.table(
  perturbation_clusters[order(perturbation_clusters$perturbation_cluster, perturbation_clusters$column_order), ],
  file.path(out_dir, paste0("pfibox_pa_model_clustered_perturbation_clusters_", run_tag, cluster_output_suffix, ".tsv")),
  sep = "\t", quote = FALSE, row.names = FALSE
)
utils::write.table(
  block_summary[order(block_summary$reporter_cluster, block_summary$perturbation_cluster), ],
  file.path(out_dir, paste0("pfibox_pa_model_clustered_block_summary_", run_tag, cluster_output_suffix, ".tsv")),
  sep = "\t", quote = FALSE, row.names = FALSE
)

write_paper_artifact(file.path(out_dir, paste0("pfibox_paper_diagnostics_", run_tag, cluster_output_suffix, ".rds")),
  results = list(default = res), annotations = list(metadata = metadata,
    reporter_clusters = reporter_clusters, perturbation_clusters = perturbation_clusters),
  diagnostics = list(clustering = list(mat = mat, ordered_mat = ordered_mat, row_score = row_score,
    row_hc = row_hc, col_hc = col_hc, row_order = row_order, col_order = col_order,
    row_cluster = row_cluster, col_cluster = col_cluster, block_summary = block_summary)),
  provenance = list(run_tag = run_tag, linkage = cluster_linkage, distance = cluster_distance,
                    top_n_reporters = top_n_reporters, cluster_output_suffix = cluster_output_suffix))
