source(file.path("analysis", "_helpers.R"))
load_destress_package()

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package `ggplot2` is required for Binsfeld modeling-step plots.", call. = FALSE)
}
if (!requireNamespace("gridExtra", quietly = TRUE)) {
  stop("Package `gridExtra` is required for Binsfeld modeling-step plots.", call. = FALSE)
}

ggplot2 <- asNamespace("ggplot2")
gridExtra <- asNamespace("gridExtra")

load(analysis_path("data", "binsfeld_reporter_data.rda"))
out_dir <- analysis_output_dir("binsfeld_modeling_steps")
evc_method_dir <- analysis_output_dir("binsfeld_evc_calibrated")

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
    replicate = "replicate",
    growth_covariates = "replicate",
    numeric_covariates = "dose_level"
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
  replicate = "replicate",
  numeric_covariates = "dose_level"
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
  growth_covariates = "replicate",
  numeric_covariates = "dose_level",
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

promoter_order <- c("acrABp", "marRABp", "micFp", "ompFp", "robp", "soxSp", "tolCp")
parameter_promoter_order <- c("EVC", promoter_order)
plot_growth_parameters <- growth_parameters[
  growth_parameters$promoter %in% parameter_promoter_order,
  ,
  drop = FALSE
]

growth_long <- rbind(
  data.frame(
    promoter = plot_growth_parameters$promoter,
    estimate = plot_growth_parameters$alpha_raw,
    estimate_type = "Raw promoter slope",
    stringsAsFactors = FALSE
  ),
  data.frame(
    promoter = plot_growth_parameters$promoter,
    estimate = plot_growth_parameters$alpha_shrunk,
    estimate_type = "Shrunken exponent",
    stringsAsFactors = FALSE
  )
)
growth_long$promoter <- factor(growth_long$promoter, levels = rev(parameter_promoter_order))

growth_plot <- ggplot2$ggplot(
  growth_long,
  ggplot2$aes(estimate, promoter, color = estimate_type, shape = estimate_type)
) +
  ggplot2$geom_vline(xintercept = 1, color = "#9ca3af", linetype = "dashed", linewidth = 0.4) +
  ggplot2$geom_vline(
    xintercept = unique(growth_parameters$alpha_global)[1],
    color = "#111827",
    linetype = "dotted",
    linewidth = 0.5
  ) +
  ggplot2$geom_point(size = 2.8, alpha = 0.9, position = ggplot2$position_dodge(width = 0.45)) +
  ggplot2$scale_color_manual(values = c("Raw promoter slope" = "#2563eb", "Shrunken exponent" = "#dc2626")) +
  ggplot2$theme_light(base_size = 10) +
  ggplot2$theme(
    panel.grid.minor = ggplot2$element_blank(),
    legend.position = "bottom"
  ) +
  ggplot2$labs(
    x = expression(alpha[g]),
    y = "Reporter promoter",
    color = NULL,
    shape = NULL
  )

background_calibration <- attr(evc_huber, "destress")$background_fit
if (!is.null(background_calibration)) {
  utils::write.table(
    background_calibration,
    file.path(out_dir, "binsfeld_evc_background_calibration.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
}

evc_slope_plot <- ggplot2$ggplot(
  background_calibration,
  ggplot2$aes(slope, factor(promoter, levels = rev(parameter_promoter_order)))
) +
  ggplot2$geom_vline(xintercept = 1, color = "#9ca3af", linetype = "dashed", linewidth = 0.4) +
  ggplot2$geom_point(size = 2.8, color = "#059669", alpha = 0.9) +
  ggplot2$theme_light(base_size = 10) +
  ggplot2$theme(
    panel.grid.minor = ggplot2$element_blank()
  ) +
  ggplot2$labs(
    x = "EV calibration slope",
    y = "Reporter promoter"
  )

parameter_theme <- ggplot2$theme_light(base_size = 12) +
  ggplot2$theme(
    panel.grid.minor = ggplot2$element_blank(),
    axis.title.x = ggplot2$element_text(size = 18, margin = ggplot2$margin(t = 8)),
    axis.title.y = ggplot2$element_text(size = 16),
    axis.text = ggplot2$element_text(size = 13),
    plot.title = ggplot2$element_text(face = "bold", hjust = 0.5, size = 14),
    legend.position = "bottom"
  )

parameter_alpha_plot <- ggplot2$ggplot(
  growth_long,
  ggplot2$aes(estimate, promoter, color = estimate_type, shape = estimate_type)
) +
  ggplot2$geom_vline(
    xintercept = unique(growth_parameters$alpha_global)[1],
    color = "#111827",
    linetype = "dotted",
    linewidth = 0.45
  ) +
  ggplot2$geom_vline(xintercept = 1, color = "#9ca3af", linetype = "dashed", linewidth = 0.4) +
  ggplot2$geom_point(size = 3.0, alpha = 0.9, position = ggplot2$position_dodge(width = 0.45)) +
  ggplot2$scale_color_manual(
    values = c(
      "Raw promoter slope" = "#2563eb",
      "Shrunken exponent" = "#dc2626"
    ),
    drop = FALSE
  ) +
  parameter_theme +
  ggplot2$labs(
    title = "Growth-response exponent",
    x = expression(hat(alpha)[a]),
    y = "Reporter promoter",
    color = NULL,
    shape = NULL
  )

parameter_gamma_plot <- ggplot2$ggplot(
  background_calibration,
  ggplot2$aes(slope, factor(promoter, levels = rev(parameter_promoter_order)))
) +
  ggplot2$geom_vline(xintercept = 1, color = "#9ca3af", linetype = "dashed", linewidth = 0.4) +
  ggplot2$geom_point(size = 3.0, color = "#059669", alpha = 0.9) +
  ggplot2$scale_y_discrete(limits = rev(parameter_promoter_order), drop = FALSE) +
  parameter_theme +
  ggplot2$theme(
    axis.title.y = ggplot2$element_blank(),
    axis.text.y = ggplot2$element_blank(),
    axis.ticks.y = ggplot2$element_blank(),
    legend.position = "none"
  ) +
  ggplot2$labs(
    title = "EV calibration slope",
    x = expression(hat(gamma)[1 * a]),
    y = NULL
  )

parameter_plot <- gridExtra$arrangeGrob(
  parameter_alpha_plot,
  parameter_gamma_plot,
  ncol = 2,
  widths = c(1.05, 1)
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
response_scatter_data <- response_construction

union_file <- file.path(evc_method_dir, "evc_huber_comparison_to_binsfeld_and_default.tsv")
if (file.exists(union_file)) {
  union_pairs <- read.delim(union_file, check.names = FALSE, stringsAsFactors = FALSE)
  union_pairs <- union_pairs[
    union_pairs$binsfeld_hit | union_pairs$modeled_hit | union_pairs$evc_huber_hit,
    ,
    drop = FALSE
  ]
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
matched$promoter <- factor(matched$promoter, levels = promoter_order)
matched$compound <- factor(matched$compound, levels = compound_order)
response_construction$promoter <- factor(response_construction$promoter, levels = promoter_order)
response_construction$compound <- factor(response_construction$compound, levels = compound_order)
response_scatter_data$promoter <- factor(response_scatter_data$promoter, levels = promoter_order)

scatter_long <- rbind(
	  data.frame(
	    promoter = response_scatter_data$promoter,
	    compound = response_scatter_data$compound,
	    x = response_scatter_data$raw_response,
	    y = response_scatter_data$modeled_response,
	    comparison = "DStressR without EV control",
	    stringsAsFactors = FALSE
	  ),
	  data.frame(
	    promoter = response_scatter_data$promoter,
	    compound = response_scatter_data$compound,
	    x = response_scatter_data$raw_response,
	    y = response_scatter_data$evc_huber_response,
	    comparison = "DStressR with EV control",
	    stringsAsFactors = FALSE
	  )
	)
scatter_long$comparison <- factor(
  scatter_long$comparison,
  levels = c(
    "DStressR without EV control",
    "DStressR with EV control"
  )
)
scatter_limit <- range(c(scatter_long$x, scatter_long$y), finite = TRUE)
scatter_pad <- diff(scatter_limit) * 0.04
scatter_limit <- scatter_limit + c(-scatter_pad, scatter_pad)
promoter_colors <- c(
  acrABp = "#2563eb",
  marRABp = "#dc2626",
  micFp = "#7c3aed",
  ompFp = "#ea580c",
  robp = "#059669",
  soxSp = "#0891b2",
  tolCp = "#64748b"
)

response_scatter_plot <- ggplot2$ggplot(
  scatter_long,
  ggplot2$aes(x, y, color = promoter)
) +
  ggplot2$geom_abline(slope = 1, intercept = 0, color = "#111827", linewidth = 0.35, linetype = "dashed") +
  ggplot2$geom_point(size = 1.7, alpha = 0.72) +
  ggplot2$facet_wrap(ggplot2$vars(comparison), ncol = 2) +
  ggplot2$coord_equal(xlim = scatter_limit, ylim = scatter_limit) +
  ggplot2$scale_color_manual(values = promoter_colors, drop = FALSE) +
  ggplot2$theme_light(base_size = 10) +
	  ggplot2$theme(
	    panel.grid.minor = ggplot2$element_blank(),
	    strip.text = ggplot2$element_text(face = "bold", hjust = 0.5),
	    legend.position = "bottom"
	  ) +
	  ggplot2$labs(
	    x = "Raw log2 Lux response",
	    y = "DStressR response",
	    color = "Promoter"
	  )

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
    response_type = "DStressR without EV control",
    stringsAsFactors = FALSE
  ),
  data.frame(
    promoter = response_construction$promoter,
    compound = response_construction$compound,
    response = response_construction$evc_huber_response,
    response_type = "DStressR with EV control",
    stringsAsFactors = FALSE
  )
)
response_construction_long$response_type <- factor(
  response_construction_long$response_type,
  levels = c("Raw log2 Lux", "DStressR without EV control", "DStressR with EV control")
)

response_difference_long <- rbind(
  data.frame(
    promoter = response_construction$promoter,
    compound = response_construction$compound,
    difference = response_construction$modeled_minus_raw,
    comparison = "DStressR without EV control minus raw",
    stringsAsFactors = FALSE
  ),
  data.frame(
    promoter = response_construction$promoter,
    compound = response_construction$compound,
    difference = response_construction$evc_huber_minus_raw,
    comparison = "DStressR with EV control minus raw",
    stringsAsFactors = FALSE
  )
)
response_difference_long$comparison <- factor(
  response_difference_long$comparison,
  levels = c("DStressR without EV control minus raw", "DStressR with EV control minus raw")
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
	    axis.text.x = ggplot2$element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5.5),
	    axis.text.y = ggplot2$element_text(size = 8),
	    strip.text = ggplot2$element_text(face = "bold", hjust = 0.5),
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
utils::write.table(
  response_scatter_data[order(response_scatter_data$promoter, response_scatter_data$compound), ],
  file.path(out_dir, "binsfeld_response_scale_scatter_data.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

ggplot2$ggsave(file.path(out_dir, "binsfeld_growth_parameter_estimates.png"), growth_plot, width = 7.2, height = 4.8, dpi = 220)
ggplot2$ggsave(file.path(out_dir, "binsfeld_growth_parameter_estimates.pdf"), growth_plot, width = 7.2, height = 4.8)
ggplot2$ggsave(file.path(out_dir, "binsfeld_evc_background_calibration_slopes.png"), evc_slope_plot, width = 6.8, height = 4.4, dpi = 220)
ggplot2$ggsave(file.path(out_dir, "binsfeld_evc_background_calibration_slopes.pdf"), evc_slope_plot, width = 6.8, height = 4.4)
ggplot2$ggsave(file.path(out_dir, "binsfeld_response_model_parameters.png"), parameter_plot, width = 9.2, height = 4.8, dpi = 220)
ggplot2$ggsave(file.path(out_dir, "binsfeld_response_model_parameters.pdf"), parameter_plot, width = 9.2, height = 4.8)
ggplot2$ggsave(file.path(out_dir, "binsfeld_response_scale_scatter.png"), response_scatter_plot, width = 10.8, height = 5.4, dpi = 220)
ggplot2$ggsave(file.path(out_dir, "binsfeld_response_scale_scatter.pdf"), response_scatter_plot, width = 10.8, height = 5.4)
ggplot2$ggsave(file.path(out_dir, "binsfeld_raw_modeled_evc_response_heatmaps.png"), response_construction_heatmap, width = 12.5, height = 9.4, dpi = 220)
ggplot2$ggsave(file.path(out_dir, "binsfeld_raw_modeled_evc_response_heatmaps.pdf"), response_construction_heatmap, width = 12.5, height = 9.4)
ggplot2$ggsave(file.path(out_dir, "binsfeld_response_minus_raw_heatmaps.png"), response_difference_heatmap, width = 12.5, height = 6.8, dpi = 220)
ggplot2$ggsave(file.path(out_dir, "binsfeld_response_minus_raw_heatmaps.pdf"), response_difference_heatmap, width = 12.5, height = 6.8)
ggplot2$ggsave(file.path(out_dir, "binsfeld_matched_response_heatmaps.png"), matched_heatmap, width = 12.5, height = 6.8, dpi = 220)
ggplot2$ggsave(file.path(out_dir, "binsfeld_matched_response_heatmaps.pdf"), matched_heatmap, width = 12.5, height = 6.8)
ggplot2$ggsave(file.path(out_dir, "binsfeld_modeled_minus_alpha1_response_heatmap.png"), difference_heatmap, width = 12.5, height = 4.2, dpi = 220)
ggplot2$ggsave(file.path(out_dir, "binsfeld_modeled_minus_alpha1_response_heatmap.pdf"), difference_heatmap, width = 12.5, height = 4.2)

message("Wrote E. coli modeling-step figures to: ", out_dir)
