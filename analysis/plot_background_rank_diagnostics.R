source(file.path("analysis", "_helpers.R"))
load_destress_package()

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package `ggplot2` is required for rank diagnostics.", call. = FALSE)
}
if (!requireNamespace("gridExtra", quietly = TRUE)) {
  stop("Package `gridExtra` is required for rank diagnostics.", call. = FALSE)
}

ggplot2 <- asNamespace("ggplot2")
gridExtra <- asNamespace("gridExtra")

out_dir <- comparison_results_dir("background_rank")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
pkg_dir <- package_results_dir()

data_root <- analysis_data_root()
expr_path <- file.path(data_root, "03-hit_determination", "expression_table.tsv.gz")
if (!file.exists(expr_path)) {
  stop("Campylobacter expression table not found: ", expr_path, call. = FALSE)
}

expr <- read_tsv_base(expr_path)
expr <- expr[!(expr$promoter %in% c("PCJnc20", "PCjas704")), , drop = FALSE]
expr$libplate <- sub("_.*$", "", expr$srn_code)
expr$compound_model <- ifelse(expr$ProductName == "DMSO", "DMSO", expr$srn_code)

libmap <- read_tsv_base(libmap_path())
libmap$libplate <- paste0("lp", libmap[["Library plate"]])
libmap$compound <- paste(libmap$libplate, libmap[["Well"]], sep = "_")
libmap$ProductName <- ifelse(
  is.na(libmap$ProductName) | libmap$ProductName == "NA" | libmap$ProductName == "",
  libmap[["Catalog Number"]],
  libmap$ProductName
)
compound_lookup <- unique(libmap[, c("compound", "ProductName", "Catalog Number", "Target"), drop = FALSE])

clean_compound_label <- function(x, max_chars = 34) {
  x <- iconv(x, from = "", to = "ASCII//TRANSLIT")
  x[is.na(x) | x == "" | x == "NA"] <- NA_character_
  too_long <- !is.na(x) & nchar(x) > max_chars
  x[too_long] <- paste0(substr(x[too_long], 1, max_chars - 3), "...")
  x
}

assay <- prepare_assay(
  expr,
  promoter = "promoter",
  compound = "compound_model",
  control = "DMSO",
  lux = "lux_auc_until16h",
  growth = "od_at_16h",
  growth_exponent = "estimate",
  plate = "libplate",
  replicate = "replicate"
)

export_rank_results <- function(rank) {
  fit <- fit_destress(
    assay,
    technical = c("libplate", "replicate"),
    empirical_bayes = TRUE,
    interaction = FALSE,
    adjustment = "global",
    background_rank = rank
  )
  res <- results(fit)
  pairs <- data.frame(
    promoter = res$promoter,
    compound = res$compound,
    effect = res$specific_effect,
    low_rank_effect = res$low_rank_effect,
    global_effect = res$global_effect,
    pvalue = res$specific_pvalue,
    padj_global = res$specific_padj_global,
    padj_by_promoter = res$specific_padj_by_promoter,
    stringsAsFactors = FALSE
  )
  write.table(
    pairs,
    file.path(pkg_dir, paste0("destress_moderated_rank", rank, "_pair_results.tsv")),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  pairs$rank <- rank
  pairs
}

rank_results <- do.call(rbind, lapply(0:2, export_rank_results))
rank_results$hit <- is.finite(rank_results$padj_global) & rank_results$padj_global < 0.05
rank_results$direction <- ifelse(rank_results$effect < 0, "Down-regulated", "Up-regulated")

rank_summary <- do.call(rbind, lapply(split(rank_results, rank_results$rank), function(d) {
  hits <- d[d$hit, , drop = FALSE]
  data.frame(
    rank = unique(d$rank),
    significant_pairs = nrow(hits),
    significant_compounds = length(unique(hits$compound)),
    max_promoters_per_compound = if (nrow(hits)) max(table(hits$compound)) else 0,
    positive_hits = sum(hits$effect > 0),
    negative_hits = sum(hits$effect < 0),
    stringsAsFactors = FALSE
  )
}))
rank_summary <- rank_summary[order(rank_summary$rank), , drop = FALSE]
write.table(
  rank_summary,
  file.path(out_dir, "background_rank_hit_summary.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

rank0 <- rank_results[rank_results$rank == 0, , drop = FALSE]
diag <- background_rank_diagnostics(
  rank0,
  effect = "effect",
  rank_max = 8,
  permutations = 200,
  seed = 1
)
write.table(
  diag,
  file.path(out_dir, "background_rank_scree_diagnostics.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

effect_matrix <- effect_matrix_from_table(rank0, effect = "effect")
observed_matrix <- is.finite(effect_matrix)
column_mean <- colMeans(effect_matrix, na.rm = TRUE)
column_mean[!is.finite(column_mean)] <- 0
centered_matrix <- sweep(effect_matrix, 2, column_mean, "-")
centered_matrix[!observed_matrix] <- 0
sv <- svd(centered_matrix, nu = 2, nv = 2)
promoter_loadings <- data.frame(
  promoter = rep(rownames(centered_matrix), 2),
  component = rep(c("Component 1", "Component 2"), each = nrow(centered_matrix)),
  loading = c(sv$u[, 1] * sv$d[1], sv$u[, 2] * sv$d[2]),
  stringsAsFactors = FALSE
)
promoter_loadings$component <- factor(
  promoter_loadings$component,
  levels = c("Component 1", "Component 2")
)
promoter_loadings <- promoter_loadings[order(
  promoter_loadings$component,
  promoter_loadings$loading
), , drop = FALSE]
write.table(
  promoter_loadings,
  file.path(out_dir, "background_rank_promoter_loadings.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

compound_loadings <- data.frame(
  compound = rep(colnames(centered_matrix), 2),
  component = rep(c("Component 1", "Component 2"), each = ncol(centered_matrix)),
  loading = c(sv$v[, 1] * sv$d[1], sv$v[, 2] * sv$d[2]),
  stringsAsFactors = FALSE
)
compound_loadings <- merge(compound_loadings, compound_lookup, by = "compound", all.x = TRUE, sort = FALSE)
compound_loadings$label <- ifelse(
  is.na(compound_loadings$ProductName) | compound_loadings$ProductName == "" | compound_loadings$ProductName == "NA",
  compound_loadings$compound,
  compound_loadings$ProductName
)
compound_loadings$label <- clean_compound_label(compound_loadings$label)
compound_loadings$label[is.na(compound_loadings$label)] <- compound_loadings$compound[is.na(compound_loadings$label)]
compound_loadings$component <- factor(
  compound_loadings$component,
  levels = c("Component 1", "Component 2")
)
write.table(
  compound_loadings[order(compound_loadings$component, -abs(compound_loadings$loading)), ],
  file.path(out_dir, "background_rank_compound_loadings.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

overlap <- Reduce(
  merge,
  lapply(0:2, function(rank) {
    d <- rank_results[rank_results$rank == rank, c("promoter", "compound", "hit"), drop = FALSE]
    names(d)[3] <- paste0("rank_", rank)
    d
  })
)
overlap$pair_id <- paste(overlap$promoter, overlap$compound, sep = "__")
write.table(
  overlap,
  file.path(out_dir, "background_rank_hit_membership.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

compound_summary <- do.call(rbind, lapply(split(rank_results, rank_results$rank), function(d) {
  hits <- d[d$hit, , drop = FALSE]
  if (!nrow(hits)) {
    return(data.frame())
  }
  tab <- as.data.frame(table(hits$compound), stringsAsFactors = FALSE)
  names(tab) <- c("compound", "n_promoters")
  tab$rank <- unique(d$rank)
  tab
}))
compound_summary <- compound_summary[order(compound_summary$rank, -compound_summary$n_promoters), ]
write.table(
  compound_summary,
  file.path(out_dir, "background_rank_compound_breadth.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

diag_plot <- diag[diag$component <= 6, , drop = FALSE]
scree <- ggplot2$ggplot(diag_plot, ggplot2$aes(component, observed)) +
  ggplot2$geom_ribbon(
    ggplot2$aes(ymin = null_median, ymax = null_q95),
    fill = "#bfdbfe",
    alpha = 0.55
  ) +
  ggplot2$geom_line(color = "#1d4ed8", linewidth = 0.7) +
  ggplot2$geom_point(color = "#1d4ed8", size = 2) +
  ggplot2$geom_line(ggplot2$aes(y = null_q95), color = "#64748b", linetype = "dashed", linewidth = 0.45) +
  ggplot2$scale_x_continuous(breaks = diag_plot$component) +
  ggplot2$labs(
    title = "Signed low-rank diagnostic",
    subtitle = "Observed singular values versus within-compound permutation reference",
    x = "Component",
    y = "Singular value"
  ) +
  ggplot2$theme_minimal(base_size = 10) +
  ggplot2$theme(plot.title = ggplot2$element_text(face = "bold"))

summary_long <- reshape(
  rank_summary[, c("rank", "significant_pairs", "significant_compounds", "max_promoters_per_compound")],
  varying = list(c("significant_pairs", "significant_compounds", "max_promoters_per_compound")),
  v.names = "value",
  timevar = "metric",
  times = c("Significant pairs", "Significant compounds", "Max promoters per compound"),
  direction = "long"
)
summary_long$metric <- factor(
  summary_long$metric,
  levels = c("Significant pairs", "Significant compounds", "Max promoters per compound")
)
summary_long$rank_label <- paste0("k=", summary_long$rank)

hits <- ggplot2$ggplot(summary_long, ggplot2$aes(rank_label, value, fill = metric)) +
  ggplot2$geom_col(position = "dodge", width = 0.68) +
  ggplot2$geom_text(
    ggplot2$aes(label = value),
    position = ggplot2$position_dodge(width = 0.68),
    vjust = -0.25,
    size = 2.7
  ) +
  ggplot2$scale_fill_manual(values = c("#1d4ed8", "#16a34a", "#f97316")) +
  ggplot2$labs(
    title = "Hit sensitivity",
    subtitle = "Default moderated model; global BH FDR < 0.05 after rank-k background subtraction",
    x = "Background rank",
    y = "Count",
    fill = NULL
  ) +
  ggplot2$theme_minimal(base_size = 10) +
  ggplot2$theme(
    plot.title = ggplot2$element_text(face = "bold"),
    legend.position = "bottom"
  )

promoter_loadings$promoter <- factor(
  promoter_loadings$promoter,
  levels = unique(promoter_loadings$promoter[order(promoter_loadings$loading)])
)

loadings <- ggplot2$ggplot(
  promoter_loadings,
  ggplot2$aes(loading, promoter, fill = loading > 0)
) +
  ggplot2$geom_col(width = 0.72) +
  ggplot2$geom_vline(xintercept = 0, color = "#334155", linewidth = 0.35) +
  ggplot2$facet_wrap(~ component, scales = "free_x") +
  ggplot2$scale_fill_manual(values = c("TRUE" = "#1d4ed8", "FALSE" = "#f97316")) +
  ggplot2$labs(
    title = "Promoter loadings",
    subtitle = "Signed promoter scores for the first two low-rank components",
    x = "Signed loading",
    y = NULL
  ) +
  ggplot2$theme_minimal(base_size = 9) +
  ggplot2$theme(
    plot.title = ggplot2$element_text(face = "bold"),
    legend.position = "none",
    axis.text.y = ggplot2$element_text(size = 6.2)
  )

pan_stressor_ids <- c("lp7_K22", "lp1_A8", "lp7_C18")
top_compound_loadings <- do.call(rbind, lapply(split(compound_loadings, compound_loadings$component), function(d) {
  priority <- d[d$compound %in% pan_stressor_ids, , drop = FALSE]
  strongest <- d[order(-abs(d$loading)), , drop = FALSE]
  out <- rbind(priority, strongest)
  out <- out[!duplicated(out$compound), , drop = FALSE]
  utils::head(out, 14)
}))
top_compound_loadings$label <- factor(
  top_compound_loadings$label,
  levels = unique(top_compound_loadings$label[order(top_compound_loadings$loading)])
)

compound_loading_plot <- ggplot2$ggplot(
  top_compound_loadings,
  ggplot2$aes(loading, label, fill = loading > 0)
) +
  ggplot2$geom_col(width = 0.72) +
  ggplot2$geom_vline(xintercept = 0, color = "#334155", linewidth = 0.35) +
  ggplot2$facet_wrap(~ component, scales = "free_x") +
  ggplot2$scale_fill_manual(values = c("TRUE" = "#1d4ed8", "FALSE" = "#f97316")) +
  ggplot2$labs(
    title = "Compound loadings",
    subtitle = "Top signed compound scores; pan-stressor candidates included",
    x = "Signed loading",
    y = NULL
  ) +
  ggplot2$theme_minimal(base_size = 9) +
  ggplot2$theme(
    plot.title = ggplot2$element_text(face = "bold"),
    legend.position = "none",
    axis.text.y = ggplot2$element_text(size = 6.2)
  )

combined <- gridExtra$grid.arrange(
  scree,
  hits,
  loadings,
  compound_loading_plot,
  ncol = 1,
  heights = c(1.0, 0.85, 1.15, 1.15)
)

ggplot2$ggsave(
  file.path(out_dir, "background_rank_diagnostics.png"),
  combined,
  width = 8.8,
  height = 13.2,
  dpi = 300
)
ggplot2$ggsave(
  file.path(out_dir, "background_rank_diagnostics.pdf"),
  combined,
  width = 8.8,
  height = 13.2
)

message("Wrote background-rank diagnostics to: ", out_dir)
