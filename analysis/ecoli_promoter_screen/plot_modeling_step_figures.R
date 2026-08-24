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
after_stat <- ggplot2$after_stat

panel_label <- function(label, size = 18) {
  grid::textGrob(
    label,
    x = grid::unit(0, "npc"),
    y = grid::unit(0.08, "npc"),
    hjust = 0,
    vjust = 0,
    gp = grid::gpar(fontsize = size, fontface = "bold", col = "#111827")
  )
}

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
    x = "EVC calibration slope",
    y = "Reporter promoter"
  )

parameter_panel <- rbind(
  data.frame(
    promoter = as.character(growth_long$promoter),
    estimate = growth_long$estimate,
    estimate_type = as.character(growth_long$estimate_type),
    panel = "hat(alpha)[a]",
    stringsAsFactors = FALSE
  ),
  data.frame(
    promoter = background_calibration$promoter,
    estimate = background_calibration$slope,
    estimate_type = "EVC calibration slope",
    panel = "hat(gamma)[1*a]",
    stringsAsFactors = FALSE
  )
)
parameter_panel$promoter <- factor(parameter_panel$promoter, levels = rev(parameter_promoter_order))
parameter_panel$panel <- factor(
  parameter_panel$panel,
  levels = c(
    "hat(alpha)[a]",
    "hat(gamma)[1*a]"
  )
)
parameter_panel$estimate_type <- factor(
  parameter_panel$estimate_type,
  levels = c("Raw promoter slope", "Shrunken exponent", "EVC calibration slope")
)
parameter_reference_lines <- data.frame(
  panel = factor(
    c(
      "hat(alpha)[a]",
      "hat(gamma)[1*a]"
    ),
    levels = levels(parameter_panel$panel)
  ),
  xintercept = 1,
  linetype = "Fixed value 1",
  stringsAsFactors = FALSE
)
parameter_global_line <- data.frame(
  panel = factor("hat(alpha)[a]", levels = levels(parameter_panel$panel)),
  xintercept = unique(growth_parameters$alpha_global)[1],
  stringsAsFactors = FALSE
)
parameter_legend_boxes <- data.frame(
  panel = factor(
    c("hat(alpha)[a]", "hat(alpha)[a]", "hat(gamma)[1*a]"),
    levels = levels(parameter_panel$panel)
  ),
  xmin = c(0.55, 0.55, 0.70),
  xmax = c(0.985, 0.985, 0.99),
  ymin = c(7.68, 6.68, 7.68),
  ymax = c(8.32, 7.32, 8.32),
  stringsAsFactors = FALSE
)
parameter_legend <- data.frame(
  panel = factor(
    c("hat(alpha)[a]", "hat(alpha)[a]", "hat(gamma)[1*a]"),
    levels = levels(parameter_panel$panel)
  ),
  x = c(0.58, 0.58, 0.72),
  y = factor(c("EVC", "acrABp", "EVC"), levels = rev(parameter_promoter_order)),
  label = c("Raw promoter slope", "Shrunken exponent", "EVC calibration"),
  estimate_type = factor(
    c("Raw promoter slope", "Shrunken exponent", "EVC calibration slope"),
    levels = levels(parameter_panel$estimate_type)
  ),
  stringsAsFactors = FALSE
)
parameter_legend_text <- parameter_legend
parameter_legend_text$x <- parameter_legend_text$x + c(0.045, 0.045, 0.045)
parameter_plot <- ggplot2$ggplot(
  parameter_panel,
  ggplot2$aes(estimate, promoter, shape = estimate_type)
) +
  ggplot2$geom_vline(
    data = parameter_reference_lines,
    ggplot2$aes(xintercept = xintercept),
    color = "#9ca3af",
    linetype = "dashed",
    linewidth = 0.4
  ) +
  ggplot2$geom_vline(
    data = parameter_global_line,
    ggplot2$aes(xintercept = xintercept),
    color = "#111827",
    linetype = "dotted",
    linewidth = 0.45
  ) +
  ggplot2$geom_rect(
    data = parameter_legend_boxes,
    ggplot2$aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = grDevices::adjustcolor("white", alpha.f = 0.88),
    color = "#d1d5db",
    linewidth = 0.35
  ) +
  ggplot2$geom_point(size = 3.0, alpha = 0.9) +
  ggplot2$geom_point(
    data = parameter_legend,
    ggplot2$aes(x = x, y = y, shape = estimate_type),
    size = 3.4,
    inherit.aes = FALSE
  ) +
  ggplot2$geom_text(
    data = parameter_legend_text,
    ggplot2$aes(x = x, y = y, label = label),
    hjust = 0,
    vjust = 0.5,
    size = 5.0,
    color = "#111827",
    inherit.aes = FALSE
  ) +
  ggplot2$facet_wrap(
    ggplot2$vars(panel),
    ncol = 2,
    scales = "free_x",
    labeller = ggplot2$label_parsed,
    strip.position = "bottom"
  ) +
  ggplot2$scale_y_discrete(drop = FALSE) +
  ggplot2$scale_shape_manual(
    values = c(
      "Raw promoter slope" = 1,
      "Shrunken exponent" = 17,
      "EVC calibration slope" = 15
    ),
    drop = FALSE
  ) +
  ggplot2$theme_light(base_size = 13) +
  ggplot2$theme(
    panel.grid.minor = ggplot2$element_blank(),
    axis.title.x = ggplot2$element_blank(),
    axis.title.y = ggplot2$element_text(size = 16),
    axis.text = ggplot2$element_text(size = 13),
    strip.background = ggplot2$element_blank(),
    strip.placement = "outside",
    strip.text = ggplot2$element_text(face = "bold", hjust = 0.5, size = 16, color = "#111827"),
    strip.text.x.bottom = ggplot2$element_text(
      face = "bold",
      hjust = 0.5,
      size = 16,
      color = "#111827",
      margin = ggplot2$margin(t = 8, b = 3)
    ),
    legend.position = "none"
  ) +
  ggplot2$labs(
    x = NULL,
    y = "Reporter promoter",
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
evc_huber_matrix <- response_matrix(evc_huber)
names(modeled_matrix)[names(modeled_matrix) == "delta_response"] <- "modeled_response"
names(alpha1_matrix)[names(alpha1_matrix) == "delta_response"] <- "alpha1_response"
names(evc_huber_matrix)[names(evc_huber_matrix) == "delta_response"] <- "evc_huber_response"

score_keys <- unique(wt_auc[, c("well", "drug", "promoter", "replicate", "compound")])
binsfeld_scores <- binsfeld_reporter_scores[
  binsfeld_reporter_scores$strain == "WT" &
    binsfeld_reporter_scores$statistic == "Scores",
  ,
  drop = FALSE
]
binsfeld_scores <- merge(
  binsfeld_scores[, c("well", "drug", "promoter", "replicate", "value")],
  score_keys,
  by = c("well", "drug", "promoter", "replicate"),
  all = FALSE,
  sort = FALSE
)
binsfeld_adjusted_response_matrix <- stats::aggregate(
  value ~ promoter + compound,
  binsfeld_scores[
    binsfeld_scores$compound != "Water" &
      binsfeld_scores$promoter != "EVC" &
      is.finite(binsfeld_scores$value),
    ,
    drop = FALSE
  ],
  mean
)
names(binsfeld_adjusted_response_matrix)[3] <- "binsfeld_adjusted_response"

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
    alpha1_matrix[, c("promoter", "compound", "alpha1_response")],
    modeled_matrix[, c("promoter", "compound", "modeled_response")],
    evc_huber_matrix[, c("promoter", "compound", "evc_huber_response")],
    binsfeld_adjusted_response_matrix[, c("promoter", "compound", "binsfeld_adjusted_response")]
  )
)
response_construction$modeled_minus_alpha1 <- response_construction$modeled_response -
  response_construction$alpha1_response
response_construction$evc_huber_minus_alpha1 <- response_construction$evc_huber_response -
  response_construction$alpha1_response
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
	  x = response_scatter_data$binsfeld_adjusted_response,
	  y = response_scatter_data$modeled_response,
	  comparison = "DStressR without EVC",
	  stringsAsFactors = FALSE
	),
	data.frame(
	  promoter = response_scatter_data$promoter,
	  compound = response_scatter_data$compound,
	  x = response_scatter_data$binsfeld_adjusted_response,
	  y = response_scatter_data$evc_huber_response,
	  comparison = "DStressR with EVC",
	    stringsAsFactors = FALSE
	  )
	)
scatter_long$comparison <- factor(
  scatter_long$comparison,
  levels = c(
    "DStressR without EVC",
    "DStressR with EVC"
  )
)
scatter_limit <- range(c(scatter_long$x, scatter_long$y), finite = TRUE)
scatter_pad <- diff(scatter_limit) * 0.04
scatter_limit <- scatter_limit + c(-scatter_pad, scatter_pad)
finite_x <- abs(scatter_long$x[is.finite(scatter_long$x) & scatter_long$x != 0])
pseudo_log_sigma <- if (length(finite_x) > 0) {
  max(1, stats::quantile(finite_x, probs = 0.05, names = FALSE, na.rm = TRUE))
} else {
  1
}
promoter_colors <- c(
  acrABp = "#0072B2",
  marRABp = "#D55E00",
  micFp = "#CC79A7",
  ompFp = "#E69F00",
  robp = "#009E73",
  soxSp = "#56B4E9",
  tolCp = "#6A3D9A"
)

response_scatter_plot <- ggplot2$ggplot(
  scatter_long,
  ggplot2$aes(x, y, color = promoter)
) +
  ggplot2$geom_point(size = 2.0, alpha = 0.78) +
  ggplot2$facet_wrap(ggplot2$vars(comparison), ncol = 2) +
  ggplot2$scale_color_manual(values = promoter_colors, drop = FALSE) +
  ggplot2$scale_x_continuous(
    trans = scales::pseudo_log_trans(sigma = pseudo_log_sigma),
    breaks = c(-30000, -3000, 0, 3000, 30000)
  ) +
  ggplot2$theme_light(base_size = 13) +
  ggplot2$theme(
    panel.grid.minor = ggplot2$element_blank(),
    axis.title = ggplot2$element_text(size = 16),
    axis.text = ggplot2$element_text(size = 12),
    strip.text = ggplot2$element_text(face = "bold", hjust = 0.5, size = 14),
    legend.position = c(0.985, 0.015),
    legend.justification = c(1, 0),
    legend.background = ggplot2$element_rect(fill = grDevices::adjustcolor("white", alpha.f = 0.86), color = "#d1d5db"),
    legend.key = ggplot2$element_rect(fill = grDevices::adjustcolor("white", alpha.f = 0.86), color = NA),
    legend.margin = ggplot2$margin(4, 6, 4, 6)
  ) +
  ggplot2$labs(
    x = expression("Original response " * bar(s)[aj] * " (shown on a pseudo log-scale)"),
    y = "DStressR response",
    color = "Promoter"
  )

response_modeling_figure <- gridExtra$arrangeGrob(
  gridExtra$arrangeGrob(panel_label("a"), panel_label("b"), ncol = 2),
  parameter_plot,
  gridExtra$arrangeGrob(panel_label("c"), panel_label("d"), ncol = 2),
  response_scatter_plot,
  ncol = 1,
  heights = c(0.045, 0.90, 0.045, 1.08)
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
    response = response_construction$binsfeld_adjusted_response,
    response_type = "Published adjusted response",
    stringsAsFactors = FALSE
  ),
  data.frame(
    promoter = response_construction$promoter,
    compound = response_construction$compound,
    response = response_construction$modeled_response,
    response_type = "DStressR without EVC",
    stringsAsFactors = FALSE
  ),
  data.frame(
    promoter = response_construction$promoter,
    compound = response_construction$compound,
    response = response_construction$evc_huber_response,
    response_type = "DStressR with EVC",
    stringsAsFactors = FALSE
  )
)
response_construction_long$response_type <- factor(
  response_construction_long$response_type,
  levels = c("Published adjusted response", "DStressR without EVC", "DStressR with EVC")
)

response_distribution_data <- rbind(
  data.frame(
    response = response_scatter_data$binsfeld_adjusted_response,
    response_type = "Original response",
    stringsAsFactors = FALSE
  ),
  data.frame(
    response = response_scatter_data$modeled_response,
    response_type = "DStressR without EVC",
    stringsAsFactors = FALSE
  ),
  data.frame(
    response = response_scatter_data$evc_huber_response,
    response_type = "DStressR with EVC",
    stringsAsFactors = FALSE
  )
)
response_distribution_data <- response_distribution_data[
  is.finite(response_distribution_data$response),
  ,
  drop = FALSE
]
response_distribution_data$response_type <- factor(
  response_distribution_data$response_type,
  levels = c("Original response", "DStressR without EVC", "DStressR with EVC")
)
response_distribution_summary <- stats::aggregate(
  response ~ response_type,
  response_distribution_data,
  function(z) {
    stats::quantile(z, probs = c(0.1, 0.5, 0.9), names = FALSE, na.rm = TRUE)
  }
)
response_distribution_summary <- do.call(
  data.frame,
  response_distribution_summary
)
names(response_distribution_summary) <- c("response_type", "q10", "median", "q90")
response_distribution_summary$response_type <- factor(
  response_distribution_summary$response_type,
  levels = levels(response_distribution_data$response_type)
)

distribution_theme <- ggplot2$theme_light(base_size = 11) +
  ggplot2$theme(
    panel.grid.minor = ggplot2$element_blank(),
    panel.grid.major.x = ggplot2$element_blank(),
    axis.title = ggplot2$element_text(size = 12),
    axis.text = ggplot2$element_text(size = 10, color = "#374151"),
    plot.margin = ggplot2$margin(5, 9, 5, 9),
    legend.position = "none"
  )

make_response_distribution_plot <- function(data, summary, response_type, fill, x_label, x_scale = NULL) {
  panel_data <- data[data$response_type == response_type, , drop = FALSE]
  panel_summary <- summary[summary$response_type == response_type, , drop = FALSE]
  bin_width <- diff(range(panel_data$response, finite = TRUE)) / 34
  if (!is.finite(bin_width) || bin_width <= 0) {
    bin_width <- 1
  }
  plot <- ggplot2$ggplot(panel_data, ggplot2$aes(response)) +
    ggplot2$annotate(
      "rect",
      xmin = panel_summary$q10,
      xmax = panel_summary$q90,
      ymin = -Inf,
      ymax = Inf,
      fill = grDevices::adjustcolor(fill, alpha.f = 0.10),
      color = NA
    ) +
    ggplot2$geom_histogram(
      ggplot2$aes(y = after_stat(density)),
      binwidth = bin_width,
      boundary = 0,
      fill = grDevices::adjustcolor(fill, alpha.f = 0.42),
      color = "white",
      linewidth = 0.25
    ) +
    ggplot2$geom_density(color = fill, linewidth = 0.95, adjust = 1.05) +
    ggplot2$geom_vline(xintercept = 0, color = "#9ca3af", linewidth = 0.35, linetype = "dashed") +
    ggplot2$geom_vline(xintercept = panel_summary$median, color = "#111827", linewidth = 0.45) +
    distribution_theme +
    ggplot2$labs(x = x_label, y = "Density")
  if (!is.null(x_scale)) {
    plot <- plot + x_scale
  }
  plot
}

signed_pseudolog <- function(x, sigma) {
  sign(x) * log10(1 + abs(x) / sigma)
}
original_display_data <- response_distribution_data
original_display_data$response <- signed_pseudolog(original_display_data$response, pseudo_log_sigma)
original_display_data <- original_display_data[
  original_display_data$response_type == "Original response",
  ,
  drop = FALSE
]
original_display_summary <- response_distribution_summary
original_display_summary$q10 <- signed_pseudolog(original_display_summary$q10, pseudo_log_sigma)
original_display_summary$median <- signed_pseudolog(original_display_summary$median, pseudo_log_sigma)
original_display_summary$q90 <- signed_pseudolog(original_display_summary$q90, pseudo_log_sigma)
original_display_break_values <- c(-30000, -3000, -300, 0, 300, 3000, 30000)
original_display_scale <- ggplot2$scale_x_continuous(
  breaks = signed_pseudolog(original_display_break_values, pseudo_log_sigma),
  labels = c("-3e4", "-3e3", "-300", "0", "300", "3e3", "3e4")
)
original_distribution_plot <- make_response_distribution_plot(
  original_display_data,
  original_display_summary,
  "Original response",
  "#4b5563",
  expression("Original response " * bar(s)[aj] * " (pseudo-log scale)"),
  x_scale = original_display_scale
)
modeled_distribution_plot <- make_response_distribution_plot(
  response_distribution_data,
  response_distribution_summary,
  "DStressR without EVC",
  "#2563eb",
  expression("DStressR response " * widehat(y)[aj])
)
evc_distribution_plot <- make_response_distribution_plot(
  response_distribution_data,
  response_distribution_summary,
  "DStressR with EVC",
  "#059669",
  expression("DStressR response with EVC " * widehat(y)[aj])
)
response_distribution_figure <- gridExtra$arrangeGrob(
  gridExtra$arrangeGrob(panel_label("a", size = 15), panel_label("b", size = 15), panel_label("c", size = 15), ncol = 3),
  gridExtra$arrangeGrob(original_distribution_plot, modeled_distribution_plot, evc_distribution_plot, ncol = 3),
  ncol = 1,
  heights = c(0.06, 1)
)

response_difference_long <- rbind(
  data.frame(
    promoter = response_construction$promoter,
    compound = response_construction$compound,
    difference = response_construction$modeled_minus_alpha1,
    comparison = "DStressR without EVC minus alpha=1",
    stringsAsFactors = FALSE
  ),
  data.frame(
    promoter = response_construction$promoter,
    compound = response_construction$compound,
    difference = response_construction$evc_huber_minus_alpha1,
    comparison = "DStressR with EVC minus alpha=1",
    stringsAsFactors = FALSE
  )
)
response_difference_long$comparison <- factor(
  response_difference_long$comparison,
  levels = c("DStressR without EVC minus raw", "DStressR with EVC minus raw")
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
utils::write.table(
  response_distribution_summary,
  file.path(out_dir, "binsfeld_response_distribution_summary.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

ggplot2$ggsave(file.path(out_dir, "binsfeld_growth_parameter_estimates.png"), growth_plot, width = 7.2, height = 4.8, dpi = 220)
ggplot2$ggsave(file.path(out_dir, "binsfeld_growth_parameter_estimates.pdf"), growth_plot, width = 7.2, height = 4.8)
ggplot2$ggsave(file.path(out_dir, "binsfeld_evc_background_calibration_slopes.png"), evc_slope_plot, width = 6.8, height = 4.4, dpi = 220)
ggplot2$ggsave(file.path(out_dir, "binsfeld_evc_background_calibration_slopes.pdf"), evc_slope_plot, width = 6.8, height = 4.4)
ggplot2$ggsave(file.path(out_dir, "binsfeld_response_model_parameters.png"), parameter_plot, width = 10.0, height = 4.8, dpi = 220)
ggplot2$ggsave(file.path(out_dir, "binsfeld_response_model_parameters.pdf"), parameter_plot, width = 10.0, height = 4.8)
ggplot2$ggsave(file.path(out_dir, "binsfeld_response_scale_scatter.png"), response_scatter_plot, width = 10.8, height = 5.4, dpi = 220)
ggplot2$ggsave(file.path(out_dir, "binsfeld_response_scale_scatter.pdf"), response_scatter_plot, width = 10.8, height = 5.4)
ggplot2$ggsave(file.path(out_dir, "binsfeld_response_modeling_combined.png"), response_modeling_figure, width = 10.8, height = 10.5, dpi = 220)
ggplot2$ggsave(file.path(out_dir, "binsfeld_response_modeling_combined.pdf"), response_modeling_figure, width = 10.8, height = 10.5)
ggplot2$ggsave(file.path(out_dir, "binsfeld_response_distribution_histograms.png"), response_distribution_figure, width = 11.0, height = 3.8, dpi = 300, bg = "white")
ggplot2$ggsave(file.path(out_dir, "binsfeld_response_distribution_histograms.pdf"), response_distribution_figure, width = 11.0, height = 3.8, bg = "white")
ggplot2$ggsave(file.path(out_dir, "binsfeld_raw_modeled_evc_response_heatmaps.png"), response_construction_heatmap, width = 12.5, height = 9.4, dpi = 220)
ggplot2$ggsave(file.path(out_dir, "binsfeld_raw_modeled_evc_response_heatmaps.pdf"), response_construction_heatmap, width = 12.5, height = 9.4)
ggplot2$ggsave(file.path(out_dir, "binsfeld_response_minus_raw_heatmaps.png"), response_difference_heatmap, width = 12.5, height = 6.8, dpi = 220)
ggplot2$ggsave(file.path(out_dir, "binsfeld_response_minus_raw_heatmaps.pdf"), response_difference_heatmap, width = 12.5, height = 6.8)
ggplot2$ggsave(file.path(out_dir, "binsfeld_matched_response_heatmaps.png"), matched_heatmap, width = 12.5, height = 6.8, dpi = 220)
ggplot2$ggsave(file.path(out_dir, "binsfeld_matched_response_heatmaps.pdf"), matched_heatmap, width = 12.5, height = 6.8)
ggplot2$ggsave(file.path(out_dir, "binsfeld_modeled_minus_alpha1_response_heatmap.png"), difference_heatmap, width = 12.5, height = 4.2, dpi = 220)
ggplot2$ggsave(file.path(out_dir, "binsfeld_modeled_minus_alpha1_response_heatmap.pdf"), difference_heatmap, width = 12.5, height = 4.2)

message("Wrote E. coli modeling-step figures to: ", out_dir)
