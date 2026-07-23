source(file.path("analysis", "_helpers.R"))
load_destress_package()

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package `ggplot2` is required for Binsfeld modeling-step plots.", call. = FALSE)
}

ggplot2 <- asNamespace("ggplot2")

load(analysis_path("data", "binsfeld_reporter_data.rda"))
out_dir <- analysis_output_dir("binsfeld_modeling_steps")
three_method_dir <- analysis_output_dir("binsfeld_three_method")

wt_auc <- binsfeld_reporter_auc[
  binsfeld_reporter_auc$strain == "WT" &
    binsfeld_reporter_auc$removed == "No",
]

prepare_binsfeld <- function(growth_exponent) {
  prepare_assay(
    wt_auc,
    promoter = "promoter",
    compound = "compound",
    control = "Water",
    lux = "lux_auc",
    growth = "od_auc",
    growth_exponent = growth_exponent,
    batch = "dose_level",
    replicate = "replicate"
  )
}

wt_auc$raw_log2_lux <- log2(wt_auc$lux_auc + 1e-8)
raw <- prepare_assay(
  wt_auc,
  promoter = "promoter",
  compound = "compound",
  control = "Water",
  response = "raw_log2_lux",
  batch = "dose_level",
  replicate = "replicate"
)
modeled <- prepare_binsfeld("estimate")
alpha1 <- prepare_binsfeld(1)
evc_huber <- prepare_assay(
  wt_auc,
  promoter = "promoter",
  compound = "compound",
  control = "Water",
  lux = "lux_auc",
  growth = "od_auc",
  growth_exponent = "estimate",
  batch = "dose_level",
  replicate = "replicate",
  background_promoter = "EVC",
  background_method = "huber",
  background_by = c("compound", "dose_level", "replicate")
)

growth_parameters <- attr(modeled, "destress")$growth_exponent_fit
utils::write.table(
  growth_parameters,
  file.path(out_dir, "binsfeld_growth_parameter_estimates.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

growth_long <- rbind(
  data.frame(
    promoter = growth_parameters$promoter,
    estimate = growth_parameters$alpha_raw,
    estimate_type = "Raw promoter slope",
    stringsAsFactors = FALSE
  ),
  data.frame(
    promoter = growth_parameters$promoter,
    estimate = growth_parameters$alpha_shrunk,
    estimate_type = "Shrunken exponent",
    stringsAsFactors = FALSE
  )
)
growth_long$promoter <- factor(growth_long$promoter, levels = growth_parameters$promoter)

growth_plot <- ggplot2$ggplot(
  growth_long,
  ggplot2$aes(promoter, estimate, color = estimate_type, shape = estimate_type)
) +
  ggplot2$geom_hline(yintercept = 1, color = "#9ca3af", linetype = "dashed", linewidth = 0.4) +
  ggplot2$geom_hline(
    yintercept = unique(growth_parameters$alpha_global)[1],
    color = "#111827",
    linetype = "dotted",
    linewidth = 0.5
  ) +
  ggplot2$geom_point(size = 2.8, alpha = 0.9, position = ggplot2$position_dodge(width = 0.45)) +
  ggplot2$scale_color_manual(values = c("Raw promoter slope" = "#2563eb", "Shrunken exponent" = "#dc2626")) +
  ggplot2$theme_light(base_size = 10) +
  ggplot2$theme(
    panel.grid.minor = ggplot2$element_blank(),
    plot.title = ggplot2$element_text(face = "bold"),
    axis.text.x = ggplot2$element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  ) +
  ggplot2$labs(
    title = "E. coli WT control-well growth-response exponents",
    subtitle = "Dashed line: fixed alpha=1; dotted line: global control-well slope",
    x = "Reporter promoter",
    y = expression(alpha[g]),
    color = NULL,
    shape = NULL
  )

response_matrix <- function(assay) {
  tab <- assay[, c("promoter", "compound", ".response")]
  means <- stats::aggregate(.response ~ promoter + compound, tab, mean, na.rm = TRUE)
  control <- means[means$compound == "Water", c("promoter", ".response")]
  names(control)[2] <- "control_response"
  means <- merge(means, control, by = "promoter", all.x = TRUE, sort = FALSE)
  means$delta_response <- means$.response - means$control_response
  means <- means[means$compound != "Water" & means$promoter != "EVC", ]
  means
}

modeled_matrix <- response_matrix(modeled)
alpha1_matrix <- response_matrix(alpha1)
raw_matrix <- response_matrix(raw)
evc_huber_matrix <- response_matrix(evc_huber)
names(modeled_matrix)[names(modeled_matrix) == "delta_response"] <- "modeled_response"
names(alpha1_matrix)[names(alpha1_matrix) == "delta_response"] <- "alpha1_response"
names(raw_matrix)[names(raw_matrix) == "delta_response"] <- "raw_response"
names(evc_huber_matrix)[names(evc_huber_matrix) == "delta_response"] <- "evc_huber_response"

matched <- merge(
  modeled_matrix[, c("promoter", "compound", "modeled_response")],
  alpha1_matrix[, c("promoter", "compound", "alpha1_response")],
  by = c("promoter", "compound"),
  all = FALSE,
  sort = FALSE
)
matched$difference <- matched$modeled_response - matched$alpha1_response

response_construction <- Reduce(
  function(x, y) merge(x, y, by = c("promoter", "compound"), all = FALSE, sort = FALSE),
  list(
    raw_matrix[, c("promoter", "compound", "raw_response")],
    modeled_matrix[, c("promoter", "compound", "modeled_response")],
    evc_huber_matrix[, c("promoter", "compound", "evc_huber_response")]
  )
)
response_construction$modeled_minus_raw <- response_construction$modeled_response -
  response_construction$raw_response
response_construction$evc_huber_minus_raw <- response_construction$evc_huber_response -
  response_construction$raw_response

union_file <- file.path(three_method_dir, "binsfeld_three_method_significant_union.tsv")
if (file.exists(union_file)) {
  union_pairs <- read.delim(union_file, check.names = FALSE, stringsAsFactors = FALSE)
  keep_compounds <- sort(unique(union_pairs$compound))
  matched <- matched[matched$compound %in% keep_compounds, ]
  response_construction <- response_construction[response_construction$compound %in% keep_compounds, ]
}

matrix_wide <- stats::reshape(
  matched[, c("promoter", "compound", "modeled_response")],
  idvar = "promoter",
  timevar = "compound",
  direction = "wide"
)
compound_cols <- setdiff(names(matrix_wide), "promoter")
compound_order <- sub("^modeled_response[.]", "", compound_cols)
if (length(compound_cols) > 2) {
  mat <- as.matrix(matrix_wide[, compound_cols, drop = FALSE])
  rownames(mat) <- matrix_wide$promoter
  mat[is.na(mat)] <- 0
  compound_order <- sub("^modeled_response[.]", "", colnames(mat)[stats::hclust(stats::dist(t(mat)))$order])
}
promoter_order <- c("acrABp", "marRABp", "micFp", "ompFp", "robp", "soxSp", "tolCp")
matched$promoter <- factor(matched$promoter, levels = promoter_order)
matched$compound <- factor(matched$compound, levels = compound_order)
response_construction$promoter <- factor(response_construction$promoter, levels = promoter_order)
response_construction$compound <- factor(response_construction$compound, levels = compound_order)

heat_long <- rbind(
  data.frame(
    promoter = matched$promoter,
    compound = matched$compound,
    response = matched$modeled_response,
    response_type = "Modeled response",
    stringsAsFactors = FALSE
  ),
  data.frame(
    promoter = matched$promoter,
    compound = matched$compound,
    response = matched$alpha1_response,
    response_type = "Alpha=1 response",
    stringsAsFactors = FALSE
  )
)
heat_long$response_type <- factor(heat_long$response_type, levels = c("Modeled response", "Alpha=1 response"))

response_construction_long <- rbind(
  data.frame(
    promoter = response_construction$promoter,
    compound = response_construction$compound,
    response = response_construction$raw_response,
    response_type = "Raw log2 Lux",
    stringsAsFactors = FALSE
  ),
  data.frame(
    promoter = response_construction$promoter,
    compound = response_construction$compound,
    response = response_construction$modeled_response,
    response_type = "DStressR alpha_g",
    stringsAsFactors = FALSE
  ),
  data.frame(
    promoter = response_construction$promoter,
    compound = response_construction$compound,
    response = response_construction$evc_huber_response,
    response_type = "DStressR alpha_g + EVC-Huber",
    stringsAsFactors = FALSE
  )
)
response_construction_long$response_type <- factor(
  response_construction_long$response_type,
  levels = c("Raw log2 Lux", "DStressR alpha_g", "DStressR alpha_g + EVC-Huber")
)

response_difference_long <- rbind(
  data.frame(
    promoter = response_construction$promoter,
    compound = response_construction$compound,
    difference = response_construction$modeled_minus_raw,
    comparison = "DStressR alpha_g minus raw",
    stringsAsFactors = FALSE
  ),
  data.frame(
    promoter = response_construction$promoter,
    compound = response_construction$compound,
    difference = response_construction$evc_huber_minus_raw,
    comparison = "DStressR alpha_g + EVC-Huber minus raw",
    stringsAsFactors = FALSE
  )
)
response_difference_long$comparison <- factor(
  response_difference_long$comparison,
  levels = c("DStressR alpha_g minus raw", "DStressR alpha_g + EVC-Huber minus raw")
)

limit <- max(abs(c(heat_long$response, matched$difference)), na.rm = TRUE)
fill_scale <- ggplot2$scale_fill_gradient2(
  low = "#2563eb",
  mid = "white",
  high = "#dc2626",
  midpoint = 0,
  limits = c(-limit, limit),
  oob = scales::squish
)

heat_theme <- ggplot2$theme_minimal(base_size = 8) +
  ggplot2$theme(
    panel.grid = ggplot2$element_blank(),
    plot.title = ggplot2$element_text(face = "bold", size = 11),
    axis.text.x = ggplot2$element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5.5),
    axis.text.y = ggplot2$element_text(size = 8),
    strip.text = ggplot2$element_text(face = "bold"),
    legend.position = "bottom"
  )

matched_heatmap <- ggplot2$ggplot(
  heat_long,
  ggplot2$aes(compound, promoter, fill = response)
) +
  ggplot2$geom_tile(color = "white", linewidth = 0.12) +
  ggplot2$facet_wrap(ggplot2$vars(response_type), ncol = 1) +
  fill_scale +
  heat_theme +
  ggplot2$labs(
    title = "E. coli WT response matrices under two DStressR response scales",
    subtitle = "Promoter-centered compound responses; compounds are the significant union from the three-method comparison",
    x = "Compound",
    y = "Reporter promoter",
    fill = "Response"
  )

difference_limit <- max(abs(response_difference_long$difference), na.rm = TRUE)
response_construction_heatmap <- ggplot2$ggplot(
  response_construction_long,
  ggplot2$aes(compound, promoter, fill = response)
) +
  ggplot2$geom_tile(color = "white", linewidth = 0.12) +
  ggplot2$facet_wrap(ggplot2$vars(response_type), ncol = 1) +
  ggplot2$scale_fill_gradient2(
    low = "#2563eb",
    mid = "white",
    high = "#dc2626",
    midpoint = 0,
    oob = scales::squish
  ) +
  heat_theme +
  ggplot2$labs(
    title = "E. coli WT response construction",
    subtitle = "Promoter-centered compound responses compared across raw and DStressR response scales",
    x = "Compound",
    y = "Reporter promoter",
    fill = "Response"
  )

response_difference_heatmap <- ggplot2$ggplot(
  response_difference_long,
  ggplot2$aes(compound, promoter, fill = difference)
) +
  ggplot2$geom_tile(color = "white", linewidth = 0.12) +
  ggplot2$facet_wrap(ggplot2$vars(comparison), ncol = 1) +
  ggplot2$scale_fill_gradient2(
    low = "#2563eb",
    mid = "white",
    high = "#dc2626",
    midpoint = 0,
    limits = c(-difference_limit, difference_limit),
    oob = scales::squish
  ) +
  heat_theme +
  ggplot2$labs(
    title = "Change relative to the raw unadjusted response",
    subtitle = "Differences are computed after promoter-centering each response matrix against water",
    x = "Compound",
    y = "Reporter promoter",
    fill = "Difference"
  )

difference_heatmap <- ggplot2$ggplot(
  matched,
  ggplot2$aes(compound, promoter, fill = difference)
) +
  ggplot2$geom_tile(color = "white", linewidth = 0.12) +
  ggplot2$scale_fill_gradient2(
    low = "#2563eb",
    mid = "white",
    high = "#dc2626",
    midpoint = 0,
    oob = scales::squish
  ) +
  heat_theme +
  ggplot2$labs(
    title = "Influence of modeled growth adjustment in the E. coli WT screen",
    subtitle = "Difference: modeled response minus fixed alpha=1 response",
    x = "Compound",
    y = "Reporter promoter",
    fill = "Difference"
  )

utils::write.table(
  matched[order(matched$promoter, matched$compound), ],
  file.path(out_dir, "binsfeld_modeled_alpha1_response_matrix_long.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
utils::write.table(
  response_construction[order(response_construction$promoter, response_construction$compound), ],
  file.path(out_dir, "binsfeld_raw_modeled_evc_response_matrix_long.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

ggplot2$ggsave(file.path(out_dir, "binsfeld_growth_parameter_estimates.png"), growth_plot, width = 7.2, height = 4.8, dpi = 220)
ggplot2$ggsave(file.path(out_dir, "binsfeld_growth_parameter_estimates.pdf"), growth_plot, width = 7.2, height = 4.8)
ggplot2$ggsave(file.path(out_dir, "binsfeld_raw_modeled_evc_response_heatmaps.png"), response_construction_heatmap, width = 12.5, height = 9.4, dpi = 220)
ggplot2$ggsave(file.path(out_dir, "binsfeld_raw_modeled_evc_response_heatmaps.pdf"), response_construction_heatmap, width = 12.5, height = 9.4)
ggplot2$ggsave(file.path(out_dir, "binsfeld_response_minus_raw_heatmaps.png"), response_difference_heatmap, width = 12.5, height = 6.8, dpi = 220)
ggplot2$ggsave(file.path(out_dir, "binsfeld_response_minus_raw_heatmaps.pdf"), response_difference_heatmap, width = 12.5, height = 6.8)
ggplot2$ggsave(file.path(out_dir, "binsfeld_matched_response_heatmaps.png"), matched_heatmap, width = 12.5, height = 6.8, dpi = 220)
ggplot2$ggsave(file.path(out_dir, "binsfeld_matched_response_heatmaps.pdf"), matched_heatmap, width = 12.5, height = 6.8)
ggplot2$ggsave(file.path(out_dir, "binsfeld_modeled_minus_alpha1_response_heatmap.png"), difference_heatmap, width = 12.5, height = 4.2, dpi = 220)
ggplot2$ggsave(file.path(out_dir, "binsfeld_modeled_minus_alpha1_response_heatmap.pdf"), difference_heatmap, width = 12.5, height = 4.2)

message("Wrote E. coli modeling-step figures to: ", out_dir)
