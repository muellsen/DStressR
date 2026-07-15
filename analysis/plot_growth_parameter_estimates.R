#!/usr/bin/env Rscript

source(file.path("analysis", "_helpers.R"))

suppressPackageStartupMessages({
  library(ggplot2)
})

out_dir <- comparison_results_dir("growth_parameters")
input_file <- package_results_dir("growth_parameters.tsv")

if (!file.exists(input_file)) {
  stop(
    "Missing package growth-parameter output: ", input_file,
    "\nExpected a TSV exported by the package with columns including ",
    "promoter, a_raw, a_raw_se, alpha_raw, alpha_raw_se, alpha_shrunk, alpha_shrunk_se.",
    call. = FALSE
  )
}

params <- read_tsv_base(input_file)
required <- c(
  "promoter",
  "control_n",
  "a_raw",
  "a_raw_se",
  "alpha_raw",
  "alpha_raw_se",
  "alpha_shrunk",
  "alpha_shrunk_se",
  "alpha_global"
)
missing <- setdiff(required, names(params))
if (length(missing) > 0) {
  stop("Growth-parameter output is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
}

params <- params[order(params$a_raw), , drop = FALSE]
promoter_order <- params$promoter

make_panel <- function(parameter, estimate, se, reference = NA_real_) {
  data.frame(
    promoter = params$promoter,
    parameter = parameter,
    estimate = estimate,
    se = se,
    reference = reference,
    control_n = params$control_n,
    stringsAsFactors = FALSE
  )
}

plot_data <- rbind(
  make_panel("a_raw", params$a_raw, params$a_raw_se),
  make_panel("alpha_raw", params$alpha_raw, params$alpha_raw_se, 1),
  make_panel("alpha_shrunk", params$alpha_shrunk, params$alpha_shrunk_se, 1)
)
plot_data$ci_low <- plot_data$estimate - 1.96 * plot_data$se
plot_data$ci_high <- plot_data$estimate + 1.96 * plot_data$se
plot_data$promoter <- factor(plot_data$promoter, levels = promoter_order)
plot_data$parameter <- factor(plot_data$parameter, levels = c("a_raw", "alpha_raw", "alpha_shrunk"))

reference_data <- data.frame(
  parameter = factor(c("alpha_raw", "alpha_shrunk"), levels = levels(plot_data$parameter)),
  xintercept = 1
)
global_data <- data.frame(
  parameter = factor(c("alpha_raw", "alpha_shrunk"), levels = levels(plot_data$parameter)),
  xintercept = unique(params$alpha_global)[1]
)

p <- ggplot(plot_data, aes(x = estimate, y = promoter)) +
  geom_vline(
    data = reference_data,
    aes(xintercept = xintercept),
    color = "#94a3b8",
    linewidth = 0.35,
    linetype = "dashed"
  ) +
  geom_vline(
    data = global_data,
    aes(xintercept = xintercept),
    color = "#475569",
    linewidth = 0.4,
    linetype = "dotted"
  ) +
  geom_errorbarh(
    aes(xmin = ci_low, xmax = ci_high),
    height = 0,
    linewidth = 0.35,
    color = "#64748b"
  ) +
  geom_point(aes(fill = control_n), shape = 21, size = 2.8, stroke = 0.35, color = "#0f172a") +
  facet_grid(. ~ parameter, scales = "free_x", labeller = as_labeller(c(
    a_raw = "Promoter~intercept~a[g]",
    alpha_raw = "Raw~growth~exponent~alpha[g]",
    alpha_shrunk = "Shrunken~growth~exponent~alpha[g]"
  ), default = label_parsed)) +
  scale_fill_gradient(low = "#dbeafe", high = "#1d4ed8", name = "Control wells") +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.major.y = element_line(color = "#e2e8f0", linewidth = 0.25),
    panel.grid.major.x = element_line(color = "#f1f5f9", linewidth = 0.25),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "#f8fafc", color = "#cbd5e1"),
    strip.text = element_text(face = "bold"),
    axis.text.y = element_text(size = 8),
    legend.position = "bottom",
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(8, 14, 8, 10)
  ) +
  labs(
    title = "Promoter growth-normalization parameter estimates",
    subtitle = "Package output from DMSO control wells; intervals show estimate +/- 1.96 SE",
    x = "Estimate",
    y = "Promoter"
  )

write.table(
  plot_data,
  file.path(out_dir, "growth_parameter_estimates_long.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
write.table(
  params,
  file.path(out_dir, "growth_parameter_estimates.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

ggsave(
  file.path(out_dir, "growth_parameter_estimates.png"),
  p,
  width = 12.5,
  height = 8.6,
  dpi = 300,
  bg = "white"
)
ggsave(
  file.path(out_dir, "growth_parameter_estimates.pdf"),
  p,
  width = 12.5,
  height = 8.6,
  bg = "white"
)

message("Wrote growth parameter figure to: ", out_dir)
message("Promoters: ", nrow(params))
