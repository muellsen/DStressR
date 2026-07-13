#!/usr/bin/env Rscript

source(file.path("analysis", "_helpers.R"))

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
})

root <- analysis_data_root()
import_root <- file.path(root, "00-import", "Campylobacter")
out_dir <- file.path(getwd(), "analysis", "outputs", "growth_exponent", "educational_examples")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

examples <- data.frame(
  promoter = c("PCj0164c2", "PCJnc110"),
  srn_code = c("lp4_K10", "lp4_K10"),
  stringsAsFactors = FALSE
)

read_tsv_base <- function(path) {
  read.delim(path, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
}

wide_curve_to_long <- function(path, wells, value_name) {
  d <- read_tsv_base(path)
  d$.time <- as.POSIXct(d$Time, format = "%m/%d/%Y %H:%M:%S")
  if (any(!is.finite(as.numeric(d$.time)))) {
    d$.time <- as.POSIXct(d$Time, format = "%m/%d/%y %H:%M")
  }
  d$time_h <- as.numeric(d$.time - d$.time[1], units = "hours")
  wells <- intersect(wells, names(d))
  if (length(wells) == 0) {
    stop("No requested wells found in ", path, call. = FALSE)
  }
  out <- do.call(
    rbind,
    lapply(wells, function(w) {
      data.frame(
        Time = d$Time,
        time_h = d$time_h,
        well = w,
        value = as.numeric(d[[w]]),
        stringsAsFactors = FALSE
      )
    })
  )
  names(out)[names(out) == "value"] <- value_name
  out
}

resolve_raw_path <- function(file_origin, assay) {
  batch <- dirname(file_origin)
  file <- basename(file_origin)
  if (assay == "lux") {
    file <- sub("_OD\\.tsv\\.gz$", "_LUX.tsv.gz", file)
    file.path(import_root, "lux_data", batch, file)
  } else {
    file.path(import_root, "od_data", batch, file)
  }
}

safe_log2 <- function(x) log2(pmax(x, .Machine$double.eps))

libmap <- read_tsv_base(file.path(root, "00-import", "Campylobacter", "LibMap.txt"))
libmap$libplate <- paste0("lp", libmap[["Library plate"]])
libmap$srn_code <- paste(libmap$libplate, libmap[["Well"]], sep = "_")
libmap$ProductName <- ifelse(
  is.na(libmap$ProductName) | libmap$ProductName == "NA" | libmap$ProductName == "",
  libmap[["Catalog Number"]],
  libmap$ProductName
)

metadata <- read_tsv_base(file.path(root, "02-lux_expression", "complete_metadata.tsv.gz"))
expr <- read_tsv_base(file.path(root, "02-lux_expression", "expression_values.tsv.gz"))
expr <- merge(
  expr,
  libmap[, c("srn_code", "ProductName", "Catalog Number")],
  by = "srn_code",
  all.x = TRUE,
  sort = FALSE
)

alpha <- read_tsv_base(file.path(getwd(), "analysis", "outputs", "growth_exponent", "libplate_alpha", "growth_alpha_pooled_vs_adjusted_promoter_level.tsv"))

curve_panels <- list()
response_panels <- list()
summary_rows <- list()

for (ii in seq_len(nrow(examples))) {
  ex <- examples[ii, ]
  selected_meta <- metadata[
    metadata$promoter == ex$promoter &
      metadata$srn_code == ex$srn_code &
      metadata$filter.status == "keep",
  ]
  selected_meta <- selected_meta[order(selected_meta$replicate), ]
  selected_meta <- selected_meta[!duplicated(selected_meta$replicate), ]
  selected_meta <- head(selected_meta, 2)
  if (nrow(selected_meta) < 2) {
    stop("Need at least two kept replicate rows for ", ex$promoter, " / ", ex$srn_code, call. = FALSE)
  }

  ex_alpha <- alpha[alpha$promoter == ex$promoter, ]
  selected_expr <- expr[expr$experiment_id %in% selected_meta$experiment_id, ]
  dmso_codes <- libmap$srn_code[
    libmap$libplate == selected_meta$libplate[1] &
      libmap$ProductName == "DMSO"
  ]
  dmso_expr <- expr[
    expr$promoter == ex$promoter &
      expr$libplate == selected_meta$libplate[1] &
      expr$replicate %in% selected_meta$replicate &
      expr$srn_code %in% dmso_codes,
  ]

  y_df <- rbind(
    data.frame(
      experiment_id = selected_expr$experiment_id,
      replicate = selected_expr$replicate,
      srn_code = selected_expr$srn_code,
      well = selected_expr$well,
      group = "Selected compound",
      LUX.AUC_16 = selected_expr$LUX.AUC_16,
      od_16h.measured = selected_expr$od_16h.measured,
      stringsAsFactors = FALSE
    ),
    data.frame(
      experiment_id = dmso_expr$experiment_id,
      replicate = dmso_expr$replicate,
      srn_code = dmso_expr$srn_code,
      well = dmso_expr$well,
      group = "Plate DMSO",
      LUX.AUC_16 = dmso_expr$LUX.AUC_16,
      od_16h.measured = dmso_expr$od_16h.measured,
      stringsAsFactors = FALSE
    )
  )
  y_long <- rbind(
    data.frame(
      promoter = ex$promoter,
      compound = ex$srn_code,
      replicate = y_df$replicate,
      srn_code = y_df$srn_code,
      group = y_df$group,
      method = "Fixed alpha = 1",
      Y = safe_log2(y_df$LUX.AUC_16) - safe_log2(y_df$od_16h.measured),
      stringsAsFactors = FALSE
    ),
    data.frame(
      promoter = ex$promoter,
      compound = ex$srn_code,
      replicate = y_df$replicate,
      srn_code = y_df$srn_code,
      group = y_df$group,
      method = "Old pooled alpha_g",
      Y = safe_log2(y_df$LUX.AUC_16) - ex_alpha$pooled_alpha_shrunk * safe_log2(y_df$od_16h.measured),
      stringsAsFactors = FALSE
    ),
    data.frame(
      promoter = ex$promoter,
      compound = ex$srn_code,
      replicate = y_df$replicate,
      srn_code = y_df$srn_code,
      group = y_df$group,
      method = "Adjusted alpha_g",
      Y = safe_log2(y_df$LUX.AUC_16) - ex_alpha$adjusted_alpha_shrunk * safe_log2(y_df$od_16h.measured),
      stringsAsFactors = FALSE
    )
  )
  y_long$method <- factor(y_long$method, levels = c("Fixed alpha = 1", "Old pooled alpha_g", "Adjusted alpha_g"))
  response_panels[[ii]] <- y_long

  selected_product <- unique(selected_expr$ProductName)
  selected_product <- selected_product[is.finite(match(selected_product, selected_product)) & nzchar(selected_product)]
  if (length(selected_product) == 0) selected_product <- ex$srn_code
  summary_rows[[ii]] <- data.frame(
    promoter = ex$promoter,
    compound = ex$srn_code,
    product = selected_product[1],
    libplate = selected_meta$libplate[1],
    replicates = paste(selected_meta$replicate, collapse = ","),
    pooled_alpha = ex_alpha$pooled_alpha_shrunk,
    adjusted_alpha = ex_alpha$adjusted_alpha_shrunk,
    stringsAsFactors = FALSE
  )

  for (jj in seq_len(nrow(selected_meta))) {
    m <- selected_meta[jj, ]
    replicate_dmso <- metadata[
      metadata$promoter == ex$promoter &
        metadata$libplate == m$libplate &
        metadata$replicate == m$replicate &
        metadata$srn_code %in% dmso_codes &
        metadata$filter.status == "keep",
    ]
    wells <- unique(c(m$well, replicate_dmso$well))
    od_path <- resolve_raw_path(m$file_origin, "od")
    lux_path <- resolve_raw_path(m$file_origin, "lux")
    if (!file.exists(od_path) || !file.exists(lux_path)) {
      stop("Missing raw curve files for ", m$experiment_id, call. = FALSE)
    }
    od_long <- wide_curve_to_long(od_path, wells, "value")
    lux_long <- wide_curve_to_long(lux_path, wells, "value")
    od_long$assay <- "Growth OD"
    lux_long$assay <- "Luminescence"
    names(od_long)[names(od_long) == "value"] <- "measurement"
    names(lux_long)[names(lux_long) == "value"] <- "measurement"
    curves <- rbind(od_long, lux_long)
    curves$promoter <- ex$promoter
    curves$compound <- ex$srn_code
    curves$product <- selected_product[1]
    curves$libplate <- m$libplate
    curves$replicate <- m$replicate
    curves$curve_group <- ifelse(curves$well == m$well, "Selected compound", "Plate DMSO")
    curves$curve_id <- paste(curves$replicate, curves$well, curves$curve_group, sep = "_")
    curve_panels[[length(curve_panels) + 1]] <- curves
  }
}

curve_df <- do.call(rbind, curve_panels)
response_df <- do.call(rbind, response_panels)
summary_df <- do.call(rbind, summary_rows)
write.table(
  summary_df,
  file.path(out_dir, "educational_example_summary.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

curve_df$example <- paste(curve_df$promoter, curve_df$compound, curve_df$product, sep = " / ")
curve_df$panel <- paste(curve_df$promoter, curve_df$assay, sep = " - ")
curve_df$replicate <- factor(curve_df$replicate, levels = unique(curve_df$replicate))
curve_df$curve_group <- factor(curve_df$curve_group, levels = c("Plate DMSO", "Selected compound"))

p_curves <- ggplot(curve_df, aes(time_h, measurement, group = curve_id)) +
  geom_line(
    data = curve_df[curve_df$curve_group == "Plate DMSO", ],
    color = "#94a3b8",
    alpha = 0.35,
    linewidth = 0.35
  ) +
  geom_line(
    data = curve_df[curve_df$curve_group == "Selected compound", ],
    aes(color = replicate),
    linewidth = 0.9
  ) +
  facet_grid(panel ~ replicate, scales = "free_y") +
  theme_bw(base_size = 9) +
  theme(
    strip.text.y = element_text(angle = 0, hjust = 0, size = 7),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "Raw growth and Lux curves for selected examples",
    subtitle = "Selected compound: lp4_K10, Miglustat (hydrochloride). Gray curves are same-promoter, same-libplate DMSO wells.",
    x = "Hours since first measurement",
    y = "Raw measurement",
    color = "Replicate"
  )

ggsave(
  file.path(out_dir, "raw_lux_growth_curves_with_plate_dmso.png"),
  p_curves,
  width = 12,
  height = 8.5,
  dpi = 300,
  bg = "white"
)
ggsave(
  file.path(out_dir, "raw_lux_growth_curves_with_plate_dmso.pdf"),
  p_curves,
  width = 12,
  height = 8.5,
  bg = "white"
)

response_df$example <- paste(response_df$promoter, response_df$compound, sep = " / ")
response_df$point_alpha <- ifelse(response_df$group == "Selected compound", 1, 0.35)
response_df$point_size <- ifelse(response_df$group == "Selected compound", 2.4, 1.25)

p_y <- ggplot(response_df, aes(method, Y)) +
  geom_point(
    data = response_df[response_df$group == "Plate DMSO", ],
    color = "#94a3b8",
    alpha = 0.45,
    size = 1.2,
    position = position_jitter(width = 0.12, height = 0, seed = 2)
  ) +
  geom_point(
    data = response_df[response_df$group == "Selected compound", ],
    aes(color = replicate),
    size = 2.5,
    position = position_jitter(width = 0.04, height = 0, seed = 3)
  ) +
  stat_summary(
    data = response_df[response_df$group == "Plate DMSO", ],
    fun = mean,
    geom = "crossbar",
    width = 0.45,
    color = "#475569",
    linewidth = 0.25
  ) +
  facet_wrap(~ example, scales = "free_y", ncol = 1) +
  theme_bw(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 20, hjust = 1),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "How growth normalization changes raw activity estimates",
    subtitle = "Gray points are same-promoter, same-libplate DMSO wells; colored points are selected compound replicates",
    x = NULL,
    y = expression(Y == log[2](LUX.AUC) - alpha[g] %.% log[2](OD[16])),
    color = "Replicate"
  )

ggsave(
  file.path(out_dir, "raw_activity_estimates_fixed_pooled_adjusted_alpha.png"),
  p_y,
  width = 10,
  height = 6.5,
  dpi = 300,
  bg = "white"
)
ggsave(
  file.path(out_dir, "raw_activity_estimates_fixed_pooled_adjusted_alpha.pdf"),
  p_y,
  width = 10,
  height = 6.5,
  bg = "white"
)

fixed_compound <- read_tsv_base(file.path(
  getwd(),
  "analysis",
  "outputs",
  "eb_moderated_variance",
  "workflow_vs_destress_eb_promoter_compound_pvalues.tsv"
))
fixed_replicate <- read_tsv_base(file.path(
  getwd(),
  "analysis",
  "outputs",
  "eb_moderated_variance",
  "workflow_vs_destress_eb_replicate_pvalues.tsv"
))
adjusted_compound <- read_tsv_base(file.path(
  getwd(),
  "analysis",
  "outputs",
  "growth_exponent",
  "workflow_vs_destress_eb_estimated_growth_alpha_promoter_compound_pvalues.tsv"
))
adjusted_replicate <- read_tsv_base(file.path(
  getwd(),
  "analysis",
  "outputs",
  "growth_exponent",
  "workflow_vs_destress_eb_estimated_growth_alpha_replicate_pvalues.tsv"
))

example_key <- paste(examples$promoter, examples$srn_code, sep = "\r")
fixed_example <- fixed_compound[paste(fixed_compound$promoter, fixed_compound$srn_code, sep = "\r") %in% example_key, ]
adjusted_example <- adjusted_compound[
  paste(adjusted_compound$promoter, adjusted_compound$srn_code, sep = "\r") %in% example_key,
]
fixed_replicate_example <- fixed_replicate[
  paste(fixed_replicate$promoter, fixed_replicate$srn_code, sep = "\r") %in% example_key &
    fixed_replicate$replicate %in% c("r1", "r2"),
]
adjusted_replicate_example <- adjusted_replicate[
  paste(adjusted_replicate$promoter, adjusted_replicate$srn_code, sep = "\r") %in% example_key &
    adjusted_replicate$replicate %in% c("r1", "r2"),
]

downstream_df <- rbind(
  data.frame(
    promoter = fixed_example$promoter,
    compound = fixed_example$srn_code,
    replicate = fixed_example$replicate,
    method = "Median polish",
    response = "polished log2FC",
    effect = fixed_example$log2FC.polished,
    pvalue = fixed_example$pvalue,
    padj_by_promoter = fixed_example$workflow_padj_by_promoter,
    stringsAsFactors = FALSE
  ),
  data.frame(
    promoter = fixed_example$promoter,
    compound = fixed_example$srn_code,
    replicate = fixed_example$replicate,
    method = "DStressR EB, alpha = 1",
    response = "log2(LUX) - log2(OD)",
    effect = fixed_example$destress_eb_effect_centered,
    pvalue = fixed_example$destress_eb_pvalue,
    padj_by_promoter = fixed_example$destress_eb_padj_by_promoter,
    stringsAsFactors = FALSE
  ),
  data.frame(
    promoter = adjusted_example$promoter,
    compound = adjusted_example$srn_code,
    replicate = adjusted_example$replicate,
    method = "DStressR EB, adjusted alpha_g",
    response = "log2(LUX) - alpha_g log2(OD)",
    effect = adjusted_example$destress_eb_effect_centered,
    pvalue = adjusted_example$destress_eb_pvalue,
    padj_by_promoter = adjusted_example$estimated_alpha_eb_padj_by_promoter,
    stringsAsFactors = FALSE
  )
)
downstream_df$example <- paste(downstream_df$promoter, downstream_df$compound, sep = " / ")
downstream_df$method <- factor(
  downstream_df$method,
  levels = c("Median polish", "DStressR EB, alpha = 1", "DStressR EB, adjusted alpha_g")
)
downstream_df$neglog10p <- -log10(pmax(downstream_df$pvalue, .Machine$double.xmin))
downstream_df$neglog10padj <- -log10(pmax(downstream_df$padj_by_promoter, .Machine$double.xmin))
downstream_df$bh_rejected <- downstream_df$padj_by_promoter < 0.05

write.table(
  downstream_df,
  file.path(out_dir, "educational_example_downstream_pvalues.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

replicate_downstream_df <- rbind(
  data.frame(
    promoter = fixed_replicate_example$promoter,
    compound = fixed_replicate_example$srn_code,
    replicate = fixed_replicate_example$replicate,
    method = "Median polish",
    effect = fixed_replicate_example$log2FC.polished,
    pvalue = fixed_replicate_example$pvalue,
    stringsAsFactors = FALSE
  ),
  data.frame(
    promoter = fixed_replicate_example$promoter,
    compound = fixed_replicate_example$srn_code,
    replicate = fixed_replicate_example$replicate,
    method = "DStressR EB, alpha = 1",
    effect = fixed_replicate_example$destress_eb_effect_centered,
    pvalue = fixed_replicate_example$destress_eb_pvalue,
    stringsAsFactors = FALSE
  ),
  data.frame(
    promoter = adjusted_replicate_example$promoter,
    compound = adjusted_replicate_example$srn_code,
    replicate = adjusted_replicate_example$replicate,
    method = "DStressR EB, adjusted alpha_g",
    effect = adjusted_replicate_example$destress_eb_effect_centered,
    pvalue = adjusted_replicate_example$destress_eb_pvalue,
    stringsAsFactors = FALSE
  )
)
replicate_downstream_df$example <- paste(
  replicate_downstream_df$promoter,
  replicate_downstream_df$compound,
  replicate_downstream_df$replicate,
  sep = " / "
)
replicate_downstream_df$method <- factor(
  replicate_downstream_df$method,
  levels = c("Median polish", "DStressR EB, alpha = 1", "DStressR EB, adjusted alpha_g")
)
replicate_downstream_df$neglog10p <- -log10(pmax(replicate_downstream_df$pvalue, .Machine$double.xmin))

write.table(
  replicate_downstream_df,
  file.path(out_dir, "educational_example_replicate_downstream_pvalues.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

p_downstream_effect <- ggplot(downstream_df, aes(method, effect, fill = method)) +
  geom_col(width = 0.68, color = "white", linewidth = 0.25) +
  geom_hline(yintercept = 0, color = "#334155", linewidth = 0.25) +
  facet_wrap(~ example, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c("#64748b", "#2563eb", "#dc2626"), guide = "none") +
  theme_bw(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 18, hjust = 1),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "Downstream effect estimate after normalization and centering",
    subtitle = "DStressR effects are EB-centered promoter-compound effects after subtracting the compound-wide effect",
    x = NULL,
    y = "Effect used for test"
  )

p_downstream_p <- ggplot(downstream_df, aes(method, neglog10p, fill = method)) +
  geom_col(width = 0.68, color = "white", linewidth = 0.25) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, color = "#334155", linewidth = 0.25) +
  facet_wrap(~ example, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c("#64748b", "#2563eb", "#dc2626"), guide = "none") +
  theme_bw(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 18, hjust = 1),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "Raw p-values for the same example pairs",
    subtitle = "Dashed line marks p = 0.05 before BH correction",
    x = NULL,
    y = expression(-log[10](p))
  )

png(
  file.path(out_dir, "downstream_effects_and_pvalues_for_examples.png"),
  width = 13,
  height = 5.5,
  units = "in",
  res = 300,
  bg = "white"
)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))
print(p_downstream_effect, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p_downstream_p, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
dev.off()

pdf(
  file.path(out_dir, "downstream_effects_and_pvalues_for_examples.pdf"),
  width = 13,
  height = 5.5
)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))
print(p_downstream_effect, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(p_downstream_p, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
dev.off()

p_replicate <- ggplot(replicate_downstream_df, aes(method, neglog10p, fill = method)) +
  geom_col(width = 0.68, color = "white", linewidth = 0.25) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, color = "#334155", linewidth = 0.25) +
  facet_wrap(~ example, ncol = 2, scales = "free_y") +
  scale_fill_manual(values = c("#64748b", "#2563eb", "#dc2626"), guide = "none") +
  theme_bw(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 18, hjust = 1),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "Replicate-level raw p-values before promoter-wise BH aggregation",
    subtitle = "Dashed line marks p = 0.05; final pair-level tables keep the more conservative replicate-level p-value",
    x = NULL,
    y = expression(-log[10](p))
  )

ggsave(
  file.path(out_dir, "replicate_level_downstream_pvalues_for_examples.png"),
  p_replicate,
  width = 13,
  height = 7,
  dpi = 300,
  bg = "white"
)
ggsave(
  file.path(out_dir, "replicate_level_downstream_pvalues_for_examples.pdf"),
  p_replicate,
  width = 13,
  height = 7,
  bg = "white"
)

message("Wrote educational growth-normalization examples to: ", out_dir)
print(summary_df, row.names = FALSE)
