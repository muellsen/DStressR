#!/usr/bin/env Rscript

source(file.path("analysis", "_helpers.R"))

suppressPackageStartupMessages({
  library(ggplot2)
})

in_dir <- analysis_path("analysis", "outputs", "ecoli_interaction_sensitivity")
out_png <- file.path(in_dir, "interaction_hit_overlap_venn.png")
out_pdf <- file.path(in_dir, "interaction_hit_overlap_venn.pdf")
comparison_file <- file.path(in_dir, "interaction_comparison_to_primary_methods.tsv")

if (!file.exists(comparison_file)) {
  stop("Missing interaction comparison table: ", comparison_file, call. = FALSE)
}

comparison <- read.delim(comparison_file, check.names = FALSE, stringsAsFactors = FALSE)

make_panel <- function(title, left_label, right_label, left_col, right_col) {
  left <- comparison$pair_id[comparison[[left_col]]]
  right <- comparison$pair_id[comparison[[right_col]]]
  data.frame(
    panel = title,
    left_label = left_label,
    right_label = right_label,
    left_only = length(setdiff(left, right)),
    overlap = length(intersect(left, right)),
    right_only = length(setdiff(right, left)),
    stringsAsFactors = FALSE
  )
}

panels <- rbind(
  make_panel(
    "Reference vs interaction without EV",
    "Binsfeld reference",
    "Interaction without EV",
    "binsfeld_hit",
    "interaction_without_ev_hit"
  ),
  make_panel(
    "Reference vs interaction with EV",
    "Binsfeld reference",
    "Interaction with EV",
    "binsfeld_hit",
    "interaction_with_ev_hit"
  ),
  make_panel(
    "DStressR without EV vs interaction",
    "DStressR without EV",
    "Interaction without EV",
    "modeled_hit",
    "interaction_without_ev_hit"
  ),
  make_panel(
    "DStressR with EV vs interaction",
    "DStressR with EV",
    "Interaction with EV",
    "evc_huber_hit",
    "interaction_with_ev_hit"
  )
)
panel_levels <- panels$panel
panels$panel <- factor(panels$panel, levels = panel_levels)

theta <- seq(0, 2 * pi, length.out = 240)
circle <- rbind(
  data.frame(circle = "left", x = -0.55 + cos(theta), y = sin(theta)),
  data.frame(circle = "right", x = 0.55 + cos(theta), y = sin(theta))
)
circle <- merge(circle, panels[, "panel", drop = FALSE], all = TRUE)
circle$panel <- factor(circle$panel, levels = panel_levels)

label_positions <- rbind(
  data.frame(label_type = "left_only", x = -1.02, y = 0.02),
  data.frame(label_type = "overlap", x = 0, y = 0.02),
  data.frame(label_type = "right_only", x = 1.02, y = 0.02)
)
count_labels <- merge(panels, label_positions, all = TRUE)
count_labels$panel <- factor(count_labels$panel, levels = panel_levels)
count_labels$count <- ifelse(
  count_labels$label_type == "left_only",
  count_labels$left_only,
  ifelse(count_labels$label_type == "overlap", count_labels$overlap, count_labels$right_only)
)

method_labels <- rbind(
  data.frame(label_side = "left", x = -0.75, y = -1.18),
  data.frame(label_side = "right", x = 0.75, y = -1.18)
)
method_labels <- merge(panels, method_labels, all = TRUE)
method_labels$panel <- factor(method_labels$panel, levels = panel_levels)
method_labels$label <- ifelse(
  method_labels$label_side == "left",
  method_labels$left_label,
  method_labels$right_label
)

p <- ggplot() +
  geom_polygon(
    data = circle,
    aes(x, y, group = interaction(panel, circle), fill = circle),
    alpha = 0.35,
    color = "#334155",
    linewidth = 0.35
  ) +
  geom_text(
    data = count_labels,
    aes(x, y, label = count),
    size = 5,
    fontface = "bold"
  ) +
  geom_text(
    data = method_labels,
    aes(x, y, label = label),
    size = 2.9,
    lineheight = 0.92
  ) +
  facet_wrap(~ panel, ncol = 2) +
  scale_fill_manual(values = c(left = "#60a5fa", right = "#f87171"), guide = "none") +
  coord_equal(xlim = c(-1.75, 1.75), ylim = c(-1.35, 1.2), clip = "off") +
  theme_void(base_size = 9) +
  theme(
    strip.text = element_text(face = "bold", size = 9),
    plot.title = element_text(face = "bold"),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  labs(
    title = "Hit-overlap summary for the full interaction-model sensitivity analysis"
  )

ggsave(out_png, p, width = 7.2, height = 5.4, dpi = 300, bg = "white")
ggsave(out_pdf, p, width = 7.2, height = 5.4, bg = "white")

message("Wrote ", out_png)
