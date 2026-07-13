#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
})

source(file.path(getwd(), "R", "growth.R"))

root <- "/Users/cmueller/Documents/GitHub/campylobacter_stressregnet/workflow/data"
out_dir <- file.path(getwd(), "analysis", "outputs", "growth_exponent", "libplate_alpha")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

expression_file <- file.path(root, "02-lux_expression", "expression_values.tsv.gz")
libmap_file <- file.path(root, "00-import", "Campylobacter", "LibMap.txt")

read_tsv_base <- function(path) {
  read.delim(path, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
}

libmap <- read_tsv_base(libmap_file)
libmap$libplate <- paste0("lp", libmap[["Library plate"]])
libmap$srn_code <- paste(libmap$libplate, libmap[["Well"]], sep = "_")
libmap$ProductName <- ifelse(
  is.na(libmap$ProductName) | libmap$ProductName == "NA" | libmap$ProductName == "",
  libmap[["Catalog Number"]],
  libmap$ProductName
)
dmso_srn_codes <- libmap$srn_code[libmap$ProductName == "DMSO"]

expr <- read_tsv_base(expression_file)
expr <- merge(
  expr,
  libmap[, c("srn_code", "ProductName", "Catalog Number")],
  by = "srn_code",
  all.x = TRUE,
  sort = FALSE
)
expr <- expr[
  !(expr[["Catalog Number"]] %in% "DMSO noisy") &
    !(expr$promoter %in% c("PCJnc20", "PCjas704")) &
    is.finite(expr$LUX.AUC_16) &
    is.finite(expr$od_16h.measured) &
    expr$LUX.AUC_16 > 0 &
    expr$od_16h.measured > 0,
]

pooled <- estimate_growth_exponents(
  expr,
  promoter = "promoter",
  compound = "srn_code",
  lux = "LUX.AUC_16",
  growth = "od_16h.measured",
  controls = dmso_srn_codes,
  min_control_n = 20,
  shrink = TRUE,
  alpha_bounds = c(-2, 3)
)
pooled <- pooled[, c("promoter", "control_n", "alpha_raw", "alpha_raw_se", "alpha_shrunk", "alpha_shrunk_se")]
names(pooled) <- c(
  "promoter",
  "pooled_control_n",
  "pooled_alpha_raw",
  "pooled_alpha_raw_se",
  "pooled_alpha_shrunk",
  "pooled_alpha_shrunk_se"
)

adjusted <- estimate_growth_exponents(
  expr,
  promoter = "promoter",
  compound = "srn_code",
  lux = "LUX.AUC_16",
  growth = "od_16h.measured",
  covariates = c("libplate", "replicate"),
  controls = dmso_srn_codes,
  min_control_n = 20,
  shrink = TRUE,
  alpha_bounds = c(-2, 3)
)
adjusted <- adjusted[, c("promoter", "control_n", "alpha_raw", "alpha_raw_se", "alpha_shrunk", "alpha_shrunk_se")]
names(adjusted) <- c(
  "promoter",
  "adjusted_control_n",
  "adjusted_alpha_raw",
  "adjusted_alpha_raw_se",
  "adjusted_alpha_shrunk",
  "adjusted_alpha_shrunk_se"
)

promoter_level <- merge(pooled, adjusted, by = "promoter", all = TRUE, sort = FALSE)
promoter_level$delta_adjusted_minus_pooled <- promoter_level$adjusted_alpha_shrunk -
  promoter_level$pooled_alpha_shrunk
promoter_level <- promoter_level[order(-abs(promoter_level$delta_adjusted_minus_pooled)), ]
write.table(
  promoter_level,
  file.path(out_dir, "growth_alpha_pooled_vs_adjusted_promoter_level.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

per_libplate <- do.call(
  rbind,
  lapply(split(expr, expr$libplate), function(d) {
    est <- estimate_growth_exponents(
      d,
      promoter = "promoter",
      compound = "srn_code",
      lux = "LUX.AUC_16",
      growth = "od_16h.measured",
      covariates = "replicate",
      controls = dmso_srn_codes,
      min_control_n = 8,
      shrink = TRUE,
      alpha_bounds = c(-2, 3)
    )
    est$libplate <- d$libplate[1]
    est
  })
)
rownames(per_libplate) <- NULL
per_libplate$alpha_shrunk_lower <- per_libplate$alpha_shrunk - 1.96 * per_libplate$alpha_shrunk_se
per_libplate$alpha_shrunk_upper <- per_libplate$alpha_shrunk + 1.96 * per_libplate$alpha_shrunk_se
write.table(
  per_libplate,
  file.path(out_dir, "growth_alpha_by_promoter_libplate.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

plot_df <- merge(
  per_libplate[, c("promoter", "libplate", "control_n", "alpha_shrunk", "alpha_shrunk_se")],
  promoter_level[, c("promoter", "pooled_alpha_shrunk", "adjusted_alpha_shrunk")],
  by = "promoter",
  all.x = TRUE,
  sort = FALSE
)
order_df <- promoter_level[order(-promoter_level$adjusted_alpha_shrunk), ]
plot_df$promoter <- factor(plot_df$promoter, levels = order_df$promoter)
promoter_level$promoter <- factor(promoter_level$promoter, levels = order_df$promoter)

scatter_df <- promoter_level
scatter_df$label <- ifelse(
  abs(scatter_df$delta_adjusted_minus_pooled) >=
    stats::quantile(abs(scatter_df$delta_adjusted_minus_pooled), 0.85, na.rm = TRUE),
  as.character(scatter_df$promoter),
  ""
)
p_scatter <- ggplot(scatter_df, aes(pooled_alpha_shrunk, adjusted_alpha_shrunk)) +
  geom_abline(slope = 1, intercept = 0, color = "#525252", linewidth = 0.35) +
  geom_hline(yintercept = 1, color = "#b91c1c", linewidth = 0.35, linetype = "dashed") +
  geom_vline(xintercept = 1, color = "#b91c1c", linewidth = 0.35, linetype = "dashed") +
  geom_point(color = "#0f766e", size = 2.2, alpha = 0.85) +
  ggrepel::geom_text_repel(aes(label = label), min.segment.length = 0, size = 3, max.overlaps = Inf) +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(size = 13),
    plot.subtitle = element_text(size = 9)
  ) +
  coord_equal() +
  labs(
    title = "Pooled vs covariate-adjusted promoter growth exponents",
    subtitle = "Pooled ignores libplate/replicate; adjusted includes both in the DMSO slope model",
    x = expression(paste("Pooled ", alpha[g])),
    y = expression(paste("Adjusted ", alpha[g]))
  )
ggsave(
  file.path(out_dir, "growth_alpha_pooled_vs_adjusted_scatter.png"),
  p_scatter,
  width = 7.5,
  height = 6,
  dpi = 300,
  bg = "white"
)
ggsave(
  file.path(out_dir, "growth_alpha_pooled_vs_adjusted_scatter.pdf"),
  p_scatter,
  width = 6.5,
  height = 6,
  bg = "white"
)

summary_points <- rbind(
  data.frame(
    promoter = promoter_level$promoter,
    alpha = promoter_level$pooled_alpha_shrunk,
    method = "Pooled",
    stringsAsFactors = FALSE
  ),
  data.frame(
    promoter = promoter_level$promoter,
    alpha = promoter_level$adjusted_alpha_shrunk,
    method = "Adjusted",
    stringsAsFactors = FALSE
  )
)
summary_points$promoter <- factor(summary_points$promoter, levels = order_df$promoter)
summary_points$method <- factor(summary_points$method, levels = c("Pooled", "Adjusted"))

p_panel <- ggplot() +
  geom_hline(yintercept = 1, color = "#b91c1c", linewidth = 0.35, linetype = "dashed") +
  geom_point(
    data = plot_df,
    aes(promoter, alpha_shrunk, color = libplate),
    position = position_jitter(width = 0.15, height = 0, seed = 1),
    size = 1.9,
    alpha = 0.78
  ) +
  geom_point(
    data = summary_points,
    aes(promoter, alpha, shape = method),
    color = "black",
    size = 2.4,
    stroke = 0.8
  ) +
  scale_shape_manual(values = c(Pooled = 4, Adjusted = 18)) +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(
    title = "Promoter growth exponents by library plate",
    subtitle = "Colored points are promoter-by-libplate estimates; black x = old pooled estimate; black diamond = adjusted default",
    x = "Promoter, ordered by adjusted alpha_g",
    y = expression(alpha[g]),
    color = "Library plate",
    shape = "Promoter-level estimate"
  )
ggsave(
  file.path(out_dir, "growth_alpha_pooled_adjusted_libplate_panel.png"),
  p_panel,
  width = 12,
  height = 6,
  dpi = 300,
  bg = "white"
)
ggsave(
  file.path(out_dir, "growth_alpha_pooled_adjusted_libplate_panel.pdf"),
  p_panel,
  width = 12,
  height = 6,
  bg = "white"
)

p_heatmap <- ggplot(plot_df, aes(libplate, promoter, fill = alpha_shrunk)) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_gradient2(
    low = "#2563eb",
    mid = "white",
    high = "#b91c1c",
    midpoint = stats::median(plot_df$alpha_shrunk, na.rm = TRUE),
    na.value = "#f3f4f6"
  ) +
  theme_bw(base_size = 9) +
  theme(panel.grid = element_blank()) +
  labs(
    title = "Promoter-by-libplate growth exponent estimates",
    subtitle = "Each tile is alpha_g estimated from DMSO wells within one library plate, adjusted for replicate",
    x = "Library plate",
    y = "Promoter",
    fill = expression(alpha[g])
  )
ggsave(
  file.path(out_dir, "growth_alpha_libplate_heatmap.png"),
  p_heatmap,
  width = 7.5,
  height = 8,
  dpi = 300,
  bg = "white"
)
ggsave(
  file.path(out_dir, "growth_alpha_libplate_heatmap.pdf"),
  p_heatmap,
  width = 7.5,
  height = 8,
  bg = "white"
)

message("Wrote pooled-vs-libplate alpha diagnostics to: ", out_dir)
print(head(promoter_level, 10), row.names = FALSE)
