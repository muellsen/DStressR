#!/usr/bin/env Rscript

# Paper-style EDA and discovery figures for the PFIbox promoter-activity
# DStressR analysis. The clustered figure displays the fitted DStressR total
# effect matrix. The discovery figure uses the Binsfeld-style raw p-value
# histogram and volcano grammar for centered specific-effect tests.

source(file.path("analysis", "_helpers.R"))
load_destress_package()

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package `ggplot2` is required for PFIbox figures.", call. = FALSE)
}
if (!requireNamespace("gridExtra", quietly = TRUE)) {
  stop("Package `gridExtra` is required for PFIbox figures.", call. = FALSE)
}
if (!requireNamespace("scales", quietly = TRUE)) {
  stop("Package `scales` is required for PFIbox figures.", call. = FALSE)
}

ggplot2 <- asNamespace("ggplot2")
gridExtra <- asNamespace("gridExtra")
scales <- asNamespace("scales")

out_dir <- analysis_output_dir("french_pfibox", "promoter_activity_subset")
run_tag <- Sys.getenv("PFIBOX_RUN_TAG", unset = "")
if (!nzchar(run_tag)) {
  pattern <- "^pfibox_pa_destress_pair_results_.*treatment_replicate_scalable[.]tsv$"
  candidates <- list.files(out_dir, pattern = pattern, full.names = TRUE)
  if (!length(candidates)) {
    candidates <- list.files(out_dir, pattern = "^pfibox_pa_destress_pair_results_.*[.]tsv$", full.names = TRUE)
  }
  if (!length(candidates)) {
    stop("No PFIbox DStressR result table found in ", out_dir, call. = FALSE)
  }
  candidates <- candidates[order(file.info(candidates)$mtime, decreasing = TRUE)]
  result_file <- candidates[1]
  run_tag <- sub("^pfibox_pa_destress_pair_results_", "", basename(result_file))
  run_tag <- sub("[.]tsv$", "", run_tag)
} else {
  result_file <- file.path(out_dir, paste0("pfibox_pa_destress_pair_results_", run_tag, ".tsv"))
}

metadata_file <- file.path(out_dir, paste0("pfibox_pa_run_metadata_", run_tag, ".tsv"))
if (!file.exists(result_file)) {
  stop("Missing PFIbox DStressR result table: ", result_file, call. = FALSE)
}
if (!file.exists(metadata_file)) {
  stop("Missing PFIbox run metadata: ", metadata_file, call. = FALSE)
}

res <- read_tsv_base(result_file)
metadata <- read_tsv_base(metadata_file)

required <- c(
  "reporter", "perturbation", "total_effect", "specific_effect",
  "specific_pvalue", "specific_padj_by_reporter"
)
missing <- setdiff(required, names(res))
if (length(missing)) {
  stop("Missing required result columns: ", paste(missing, collapse = ", "), call. = FALSE)
}

top_n_reporters <- as.integer(Sys.getenv("PFIBOX_CLUSTER_TOP_REPORTERS", unset = "450"))
n_reporter_clusters <- as.integer(Sys.getenv("PFIBOX_REPORTER_CLUSTERS", unset = "8"))
n_perturbation_clusters <- as.integer(Sys.getenv("PFIBOX_PERTURBATION_CLUSTERS", unset = "6"))
fdr <- as.numeric(Sys.getenv("PFIBOX_FDR", unset = "0.05"))
cluster_linkage <- Sys.getenv("PFIBOX_CLUSTER_LINKAGE", unset = "complete")
if (!cluster_linkage %in% c("complete", "ward.D2", "average", "single")) {
  stop("PFIBOX_CLUSTER_LINKAGE must be one of: complete, ward.D2, average, single.", call. = FALSE)
}
cluster_distance <- Sys.getenv("PFIBOX_CLUSTER_DISTANCE", unset = "euclidean")
if (!cluster_distance %in% c("euclidean", "correlation", "cosine")) {
  stop("PFIBOX_CLUSTER_DISTANCE must be one of: euclidean, correlation, cosine.", call. = FALSE)
}
cluster_output_suffix <- Sys.getenv("PFIBOX_CLUSTER_OUTPUT_SUFFIX", unset = "")
if (nzchar(cluster_output_suffix)) {
  cluster_output_suffix <- paste0("_", gsub("[^A-Za-z0-9_-]+", "_", cluster_output_suffix))
}

paper_theme <- function(base_size = 10) {
  ggplot2$theme_light(base_size = base_size) +
    ggplot2$theme(
      panel.grid.minor = ggplot2$element_blank(),
      panel.grid.major = ggplot2$element_line(linewidth = 0.24, color = "#E5E7EB"),
      axis.title = ggplot2$element_text(color = "#111827"),
      axis.text = ggplot2$element_text(color = "#4B5563"),
      strip.background = ggplot2$element_rect(fill = "white", color = "#D1D5DB", linewidth = 0.25),
      strip.text = ggplot2$element_text(face = "bold", color = "#111827"),
      legend.background = ggplot2$element_blank(),
      legend.key = ggplot2$element_blank(),
      plot.title = ggplot2$element_blank(),
      plot.subtitle = ggplot2$element_blank()
    )
}

panel_label <- function(label, plot) {
  gridExtra$arrangeGrob(
    grid::textGrob(
      label,
      x = grid::unit(0, "npc"),
      y = grid::unit(0.20, "npc"),
      hjust = 0,
      vjust = 0,
      gp = grid::gpar(fontsize = 15, fontface = "bold", col = "#111827")
    ),
    plot,
    ncol = 1,
    heights = c(0.055, 1)
  )
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

clip_limit <- stats::quantile(abs(plot_long$effect), 0.985, na.rm = TRUE, names = FALSE)
if (!is.finite(clip_limit) || clip_limit <= 0) {
  clip_limit <- max(abs(plot_long$effect), na.rm = TRUE)
}

row_clusters_ordered <- row_cluster[row_order]
col_clusters_ordered <- col_cluster[col_order]
row_breaks <- which(row_clusters_ordered[-1] != row_clusters_ordered[-length(row_clusters_ordered)]) + 0.5
col_breaks <- which(col_clusters_ordered[-1] != col_clusters_ordered[-length(col_clusters_ordered)]) + 0.5
plot_long$row_index <- match(as.character(plot_long$reporter), rev(row_order))
plot_long$col_index <- match(as.character(plot_long$perturbation), col_order)

p_cluster <- ggplot2$ggplot(plot_long, ggplot2$aes(perturbation, reporter, fill = effect)) +
  ggplot2$geom_tile(width = 1, height = 1, color = NA) +
  ggplot2$geom_hline(yintercept = length(row_order) - row_breaks + 1, color = "white", linewidth = 0.18) +
  ggplot2$geom_vline(xintercept = col_breaks, color = "white", linewidth = 0.22) +
  ggplot2$scale_fill_gradient2(
    low = "#3f6fb5",
    mid = "white",
    high = "#c93f3f",
    midpoint = 0,
    limits = c(-clip_limit, clip_limit),
    oob = scales$squish,
    na.value = "#F3F4F6",
    name = expression(hat(Delta * r)[aj]^"tot")
  ) +
  ggplot2$labs(x = NULL, y = "Reporters clustered") +
  ggplot2$theme_minimal(base_size = 10) +
  ggplot2$theme(
    panel.grid = ggplot2$element_blank(),
    axis.text.x = ggplot2$element_text(angle = 58, hjust = 1, vjust = 1, size = 7.6),
    axis.text.y = ggplot2$element_blank(),
    axis.ticks.y = ggplot2$element_blank(),
    axis.title.y = ggplot2$element_text(size = 9, color = "#111827"),
    legend.position = "right",
    legend.title = ggplot2$element_text(size = 9),
    legend.text = ggplot2$element_text(size = 8),
    plot.margin = ggplot2$margin(0, 6, 0, 4)
  ) +
  ggplot2$guides(
    fill = ggplot2$guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barheight = grid::unit(30, "mm"),
      barwidth = grid::unit(4, "mm")
    )
  )

block_rows <- data.frame(
  reporter = rownames(mat),
  reporter_cluster = paste0("R", row_cluster[rownames(mat)]),
  stringsAsFactors = FALSE
)
block_cols <- data.frame(
  perturbation = colnames(mat),
  perturbation_cluster = paste0("C", col_cluster[colnames(mat)]),
  stringsAsFactors = FALSE
)
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

p_blocks <- ggplot2$ggplot(block_summary, ggplot2$aes(perturbation_cluster, reporter_cluster, fill = effect)) +
  ggplot2$geom_tile(color = "white", linewidth = 0.45) +
  ggplot2$scale_fill_gradient2(
    low = "#3f6fb5",
    mid = "white",
    high = "#c93f3f",
    midpoint = 0,
    limits = c(-clip_limit, clip_limit),
    oob = scales$squish,
    na.value = "#F3F4F6",
    name = expression(hat(Delta * r)[aj]^"tot")
  ) +
  ggplot2$labs(x = "Perturbation clusters", y = "Reporter clusters") +
  paper_theme(10) +
  ggplot2$theme(
    legend.position = "right",
    legend.title = ggplot2$element_text(size = 9),
    legend.text = ggplot2$element_text(size = 8),
    axis.text = ggplot2$element_text(size = 9)
  ) +
  ggplot2$guides(fill = "none")

cluster_figure <- gridExtra$arrangeGrob(
  panel_label("a", p_cluster),
  panel_label("b", p_blocks),
  ncol = 1,
  heights = c(1.35, 0.72)
)

cluster_base <- file.path(out_dir, paste0("pfibox_pa_model_clustered_total_effect_", run_tag, cluster_output_suffix))
ggplot2$ggsave(paste0(cluster_base, ".png"), cluster_figure, width = 10.4, height = 9.2, dpi = 300, bg = "white")
ggplot2$ggsave(paste0(cluster_base, ".pdf"), cluster_figure, width = 10.4, height = 9.2, bg = "white")

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

discovery <- res
discovery$specific_pvalue <- as.numeric(discovery$specific_pvalue)
discovery$specific_padj_by_reporter <- as.numeric(discovery$specific_padj_by_reporter)
discovery$specific_effect <- as.numeric(discovery$specific_effect)
discovery$hit <- is.finite(discovery$specific_padj_by_reporter) &
  discovery$specific_padj_by_reporter < fdr
discovery$status <- ifelse(discovery$hit, "Hit", "Not called")
discovery$status <- factor(discovery$status, levels = c("Not called", "Hit"))
discovery$neglog10_pvalue <- -log10(pmax(discovery$specific_pvalue, .Machine$double.xmin))
volcano_y_cap <- as.numeric(Sys.getenv("PFIBOX_VOLCANO_Y_CAP", unset = "50"))
if (!is.finite(volcano_y_cap) || volcano_y_cap <= 0) {
  volcano_y_cap <- 50
}
volcano_cap_pad <- min(0.7, volcano_y_cap * 0.02)
discovery$neglog10_pvalue_plot <- pmin(discovery$neglog10_pvalue, volcano_y_cap - volcano_cap_pad)
discovery$volcano_capped <- is.finite(discovery$neglog10_pvalue) &
  discovery$neglog10_pvalue > volcano_y_cap
discovery$label <- paste(discovery$reporter, discovery$perturbation, sep = "-")
discovery$label_hjust <- ifelse(discovery$specific_effect > 0, 1.04, -0.04)
discovery$drug <- sub("-[0-9]+$", "", discovery$perturbation)
discovery$method <- "DStressR"

top_labels <- discovery[
  discovery$hit &
    is.finite(discovery$specific_effect) &
    is.finite(discovery$neglog10_pvalue),
  ,
  drop = FALSE
]
top_labels <- top_labels[
  order(-top_labels$neglog10_pvalue, -abs(top_labels$specific_effect)),
  ,
  drop = FALSE
]
top_labels <- utils::head(top_labels, as.integer(Sys.getenv("PFIBOX_VOLCANO_LABELS", unset = "18")))
top_labels$label_display <- ifelse(
  top_labels$volcano_capped,
  paste0(top_labels$label, "\n", sprintf("%.1f", top_labels$neglog10_pvalue)),
  top_labels$label
)

drug_palette <- c(
  Ampicillin = "#0072B2",
  Azithromycin = "#D55E00",
  Cefmetazole = "#009E73",
  Chloramphenicol = "#CC79A7",
  Ciprofloxacin = "#E69F00",
  Meropenem = "#56B4E9",
  Metronidazole = "#6A3D9A",
  Minocycline = "#A6761D",
  NalidixicAcid = "#1B9E77",
  Neomycin = "#E7298A",
  Penicillin = "#7570B3",
  PolymyxinB = "#66A61E",
  Sulfamethoxazole = "#E6AB02",
  Trimethoprim = "#A6CEE3",
  TrimSulfa4 = "#B15928"
)
missing_drugs <- setdiff(sort(unique(discovery$drug)), names(drug_palette))
if (length(missing_drugs)) {
  drug_palette <- c(
    drug_palette,
    stats::setNames(grDevices::hcl.colors(length(missing_drugs), palette = "Dark 3"), missing_drugs)
  )
}

p_hist <- ggplot2$ggplot(
  discovery[is.finite(discovery$specific_pvalue), , drop = FALSE],
  ggplot2$aes(specific_pvalue)
) +
  ggplot2$geom_histogram(
    bins = 40,
    fill = "#d1d5db",
    color = "white",
    linewidth = 0.2
  ) +
  ggplot2$geom_histogram(
    data = discovery[
      discovery$hit & is.finite(discovery$specific_pvalue),
      ,
      drop = FALSE
    ],
    bins = 40,
    fill = "#2B6CB0",
    color = "white",
    linewidth = 0.18,
    alpha = 0.88
  ) +
  ggplot2$facet_wrap(ggplot2$vars(method), nrow = 1, scales = "free_y") +
  ggplot2$scale_x_continuous(breaks = seq(0, 1, by = 0.5)) +
  ggplot2$coord_cartesian(xlim = c(0, 1)) +
  paper_theme(11) +
  ggplot2$theme(
    strip.background = ggplot2$element_rect(fill = "#b3b3b3", color = "#b3b3b3"),
    strip.text = ggplot2$element_text(face = "bold", color = "white", size = 10),
    plot.margin = ggplot2$margin(4, 4, 4, 4)
  ) +
  ggplot2$labs(
    x = "Raw p-value",
    y = "Pairs"
  )

xmax <- max(abs(discovery$specific_effect[is.finite(discovery$specific_effect)]), na.rm = TRUE)
xmax <- ceiling(xmax * 1.10 * 10) / 10

p_volcano <- ggplot2$ggplot(
  discovery[is.finite(discovery$specific_effect) & is.finite(discovery$neglog10_pvalue_plot), , drop = FALSE],
  ggplot2$aes(specific_effect, neglog10_pvalue_plot)
) +
  ggplot2$geom_hline(yintercept = -log10(fdr), color = "#9ca3af", linewidth = 0.3, linetype = "dashed") +
  ggplot2$geom_hline(yintercept = volcano_y_cap, color = "#9ca3af", linewidth = 0.24, linetype = "dotted") +
  ggplot2$geom_vline(xintercept = 0, color = "#9ca3af", linewidth = 0.3) +
  ggplot2$geom_point(color = "#d1d5db", alpha = 0.42, size = 0.9) +
  ggplot2$geom_point(
    data = discovery[
      discovery$hit & is.finite(discovery$specific_effect) & is.finite(discovery$neglog10_pvalue_plot),
      ,
      drop = FALSE
    ],
    ggplot2$aes(color = drug),
    alpha = 0.95,
    size = 1.45
  ) +
  ggplot2$facet_wrap(ggplot2$vars(method), nrow = 1) +
  ggplot2$scale_color_manual(values = drug_palette, drop = FALSE, name = "Antibiotic") +
  ggplot2$scale_x_continuous(limits = c(-xmax, xmax), expand = ggplot2$expansion(mult = c(0.04, 0.1))) +
  ggplot2$scale_y_continuous(
    limits = c(0, volcano_y_cap),
    breaks = pretty(c(0, volcano_y_cap), n = 5),
    expand = ggplot2$expansion(mult = c(0.02, 0.04))
  ) +
  paper_theme(11) +
  ggplot2$theme(
    strip.background = ggplot2$element_rect(fill = "#b3b3b3", color = "#b3b3b3"),
    strip.text = ggplot2$element_text(face = "bold", color = "white", size = 10),
    legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.direction = "horizontal",
    legend.background = ggplot2$element_rect(fill = grDevices::adjustcolor("white", alpha.f = 0.86), color = "#D1D5DB", linewidth = 0.2),
    legend.key.size = grid::unit(3.5, "mm"),
    legend.text = ggplot2$element_text(size = 6.2),
    legend.title = ggplot2$element_text(size = 6.8, face = "bold"),
    plot.margin = ggplot2$margin(4, 4, 4, 4)
  ) +
  ggplot2$guides(
    color = ggplot2$guide_legend(
      ncol = 2,
      byrow = TRUE,
      title.position = "top",
      title.hjust = 0,
      override.aes = list(size = 1.8, alpha = 1)
    )
  ) +
  ggplot2$labs(
    x = expression(hat(delta)[aj]),
    y = expression(-log[10] * " raw p-value")
  )

if (nrow(top_labels) > 0 && requireNamespace("ggrepel", quietly = TRUE)) {
  p_volcano <- p_volcano +
    ggrepel::geom_text_repel(
      data = top_labels,
      ggplot2$aes(label = label_display, color = drug),
      size = 2.15,
      min.segment.length = 0,
      segment.color = "#6b7280",
      segment.size = 0.16,
      segment.alpha = 0.72,
      box.padding = 0.16,
      point.padding = 0.10,
      force = 1.8,
      force_pull = 0.12,
      max.overlaps = Inf,
      max.time = 3,
      max.iter = 10000,
      seed = 3261,
      show.legend = FALSE
    )
} else if (nrow(top_labels) > 0) {
  p_volcano <- p_volcano +
    ggplot2$geom_text(
      data = top_labels,
      ggplot2$aes(label = label_display, hjust = label_hjust, color = drug),
      size = 2.15,
      vjust = -0.55,
      check_overlap = TRUE,
      show.legend = FALSE
    )
}

discovery_figure <- gridExtra$arrangeGrob(
  panel_label("a", p_hist),
  panel_label("b", p_volcano),
  ncol = 1,
  heights = c(0.72, 1.28)
)

discovery_base <- file.path(out_dir, paste0("pfibox_pa_binsfeld_style_discovery_", run_tag))
ggplot2$ggsave(paste0(discovery_base, ".png"), discovery_figure, width = 9.4, height = 8.4, dpi = 300, bg = "white")
ggplot2$ggsave(paste0(discovery_base, ".pdf"), discovery_figure, width = 9.4, height = 8.4, bg = "white")

utils::write.table(
  top_labels[order(top_labels$specific_pvalue), ],
  file.path(out_dir, paste0("pfibox_pa_binsfeld_style_discovery_top_labels_", run_tag, ".tsv")),
  sep = "\t", quote = FALSE, row.names = FALSE
)

message("Using PFIbox DStressR run tag: ", run_tag)
message("Clustered matrix: ", nrow(mat), " reporters x ", ncol(mat), " treatments")
message("Within-reporter FDR < ", fdr, " hits: ", sum(discovery$hit, na.rm = TRUE))
message("Wrote clustered figure: ", paste0(cluster_base, ".png"))
message("Wrote discovery figure: ", paste0(discovery_base, ".png"))
