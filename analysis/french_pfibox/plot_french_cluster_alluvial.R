#!/usr/bin/env Rscript

source(file.path("analysis", "_helpers.R"))
source(file.path("analysis", "_paper_artifacts.R"))

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
})

in_dir <- file.path(
  "analysis", "outputs", "french_pfibox", "promoter_activity_analysis"
)
out_dir <- Sys.getenv("PFIBOX_FIGURE_OUTPUT_DIR", unset = in_dir)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

run_tag <- Sys.getenv("PFIBOX_RUN_TAG", "reporters2000_perturbations100_T2-18h_treatment_replicate_scalable")
suffix <- Sys.getenv("PFIBOX_CLUSTER_OUTPUT_SUFFIX", "total_effect_ward")
suffix <- if (nzchar(suffix)) paste0("_", suffix) else ""
artifact <- read_paper_artifact(file.path(in_dir,
  paste0("pfibox_paper_mapping_", run_tag, suffix, ".rds")))
mapping <- artifact$annotations$mapping
dstressr_clusters <- artifact$annotations$reporter_clusters

flow_counts <- mapping %>%
  filter(!is.na(french_major_cluster), !is.na(reporter_cluster)) %>%
  mutate(
    french_major_cluster = factor(
      french_major_cluster,
      levels = paste0("F", seq_len(12))
    ),
    reporter_cluster = factor(
      reporter_cluster,
      levels = paste0("R", seq_len(8))
    )
  ) %>%
  count(french_major_cluster, reporter_cluster, name = "n") %>%
  filter(n > 0) %>%
  arrange(reporter_cluster, french_major_cluster)

reporter_cluster_order <- dstressr_clusters %>%
  filter(!is.na(reporter_cluster), is.finite(row_order)) %>%
  group_by(reporter_cluster) %>%
  summarise(first_row = min(row_order), total_n = n(), .groups = "drop") %>%
  arrange(first_row) %>%
  mutate(
    reporter_cluster = as.character(reporter_cluster),
    display_cluster = paste0("R", row_number())
  )

if (!nrow(reporter_cluster_order)) {
  reporter_cluster_order <- tibble(
    reporter_cluster = paste0("R", seq_len(8)),
    display_cluster = paste0("R", seq_len(8)),
    total_n = as.numeric(table(mapping$reporter_cluster)[paste0("R", seq_len(8))])
  ) %>%
    filter(is.finite(total_n))
}

reporter_levels <- reporter_cluster_order$display_cluster
french_levels <- paste0("F", seq_len(12))
flow_counts <- flow_counts %>%
  left_join(
    reporter_cluster_order %>% select(reporter_cluster, display_cluster),
    by = "reporter_cluster"
  ) %>%
  mutate(
    reporter_cluster = factor(display_cluster, levels = reporter_levels),
    french_major_cluster = factor(as.character(french_major_cluster), levels = french_levels)
  ) %>%
  filter(!is.na(reporter_cluster), !is.na(french_major_cluster)) %>%
  select(-display_cluster) %>%
  arrange(reporter_cluster, french_major_cluster)

layout_axis <- function(values, gap = 8) {
  totals <- tibble(cluster = names(values), total_n = as.numeric(values)) %>%
    mutate(cluster = factor(cluster, levels = names(values))) %>%
    arrange(cluster)

  total_height <- sum(totals$total_n)
  n_gaps <- max(0, nrow(totals) - 1)
  y_top <- total_height + gap * n_gaps

  y <- y_top
  ranges <- vector("list", nrow(totals))
  for (i in seq_len(nrow(totals))) {
    y0 <- y - totals$total_n[[i]]
    ranges[[i]] <- tibble(
      cluster = as.character(totals$cluster[[i]]),
      total_n = totals$total_n[[i]],
      ymin = y0,
      ymax = y
    )
    y <- y0 - gap
  }
  bind_rows(ranges)
}

make_stack <- function(flows, side = c("left", "right"), gap = 8) {
  side <- match.arg(side)
  if (side == "left") {
    totals <- reporter_cluster_order$total_n
    names(totals) <- reporter_cluster_order$display_cluster
    blocks <- layout_axis(totals[!is.na(totals)], gap = gap)
    stack <- flows %>%
      arrange(reporter_cluster, french_major_cluster) %>%
      group_by(reporter_cluster) %>%
      mutate(
        offset = cumsum(lag(n, default = 0)),
        cluster = as.character(reporter_cluster)
      ) %>%
      ungroup()
  } else {
    totals <- tapply(flows$n, flows$french_major_cluster, sum)
    totals <- totals[french_levels]
    blocks <- layout_axis(totals[!is.na(totals)], gap = gap)
    stack <- flows %>%
      arrange(french_major_cluster, reporter_cluster) %>%
      group_by(french_major_cluster) %>%
      mutate(
        offset = cumsum(lag(n, default = 0)),
        cluster = as.character(french_major_cluster)
      ) %>%
      ungroup()
  }

  stack %>%
    left_join(blocks, by = "cluster") %>%
    mutate(
      flow_ymin = ymax - offset - n,
      flow_ymax = ymax - offset
    ) %>%
    select(french_major_cluster, reporter_cluster, n, flow_ymin, flow_ymax)
}

bezier <- function(p0, p1, p2, p3, t) {
  (1 - t)^3 * p0 + 3 * (1 - t)^2 * t * p1 +
    3 * (1 - t) * t^2 * p2 + t^3 * p3
}

make_ribbon <- function(flow_id, x0, x1, y0_low, y0_high, y1_low, y1_high) {
  t <- seq(0, 1, length.out = 60)
  cx0 <- x0 + 0.34 * (x1 - x0)
  cx1 <- x1 - 0.34 * (x1 - x0)

  top <- tibble(
    flow_id = flow_id,
    x = bezier(x0, cx0, cx1, x1, t),
    y = bezier(y0_high, y0_high, y1_high, y1_high, t)
  )
  bottom <- tibble(
    flow_id = flow_id,
    x = bezier(x1, cx1, cx0, x0, t),
    y = bezier(y1_low, y1_low, y0_low, y0_low, t)
  )
  bind_rows(top, bottom)
}

gap <- 10
left_blocks <- layout_axis(
  setNames(reporter_cluster_order$total_n, reporter_cluster_order$display_cluster),
  gap = gap
) %>%
  filter(!is.na(total_n)) %>%
  mutate(axis = "DStressR reporter clusters")

right_blocks <- layout_axis(
  tapply(flow_counts$n, flow_counts$french_major_cluster, sum)[french_levels],
  gap = gap
) %>%
  filter(!is.na(total_n)) %>%
  mutate(axis = "French time-course clusters")

left_stack <- make_stack(flow_counts, "left", gap = gap) %>%
  rename(left_ymin = flow_ymin, left_ymax = flow_ymax)
right_stack <- make_stack(flow_counts, "right", gap = gap) %>%
  rename(right_ymin = flow_ymin, right_ymax = flow_ymax)

flows <- left_stack %>%
  left_join(
    right_stack,
    by = c("french_major_cluster", "reporter_cluster", "n")
  ) %>%
  arrange(french_major_cluster, reporter_cluster) %>%
  mutate(
    flow_id = row_number(),
    french_major_cluster = factor(french_major_cluster, levels = paste0("F", seq_len(12))),
    reporter_cluster = factor(reporter_cluster, levels = paste0("R", seq_len(8)))
  )

ribbons <- bind_rows(lapply(seq_len(nrow(flows)), function(i) {
  make_ribbon(
    flow_id = flows$flow_id[[i]],
    x0 = 0.25,
    x1 = 0.75,
    y0_low = flows$left_ymin[[i]],
    y0_high = flows$left_ymax[[i]],
    y1_low = flows$right_ymin[[i]],
    y1_high = flows$right_ymax[[i]]
  )
})) %>%
  left_join(
    flows %>%
      select(flow_id, french_major_cluster, reporter_cluster, n),
    by = "flow_id"
  )

french_cols <- c(
  "#5C8D89", "#88B04B", "#C7D36F", "#F3C969",
  "#E6A157", "#D96C4A", "#C75D72", "#9F6AAE",
  "#6F7DBD", "#4D96B8", "#58A6A6", "#A7A7A7"
)
names(french_cols) <- paste0("F", seq_len(12))

block_df <- bind_rows(
  left_blocks %>% mutate(xmin = 0.14, xmax = 0.25),
  right_blocks %>% mutate(xmin = 0.75, xmax = 0.86)
) %>%
  mutate(label_y = (ymin + ymax) / 2)

reporter_cols <- grDevices::hcl.colors(length(reporter_levels), palette = "Dark 3")
names(reporter_cols) <- reporter_levels

block_cols <- c(reporter_cols, french_cols)

label_df <- block_df %>%
  mutate(
    label = paste0(cluster, " (", total_n, ")"),
    x = if_else(axis == "DStressR reporter clusters", xmin - 0.018, xmax + 0.018),
    hjust = if_else(axis == "DStressR reporter clusters", 1, 0)
  )

axis_df <- tibble(
  x = c(0.195, 0.805),
  y = max(block_df$ymax) + 28,
  label = c("DStressR reporter clusters", "French et al. time-course clusters")
)

p <- ggplot() +
  geom_polygon(
    data = ribbons,
    aes(x = x, y = y, group = flow_id, fill = reporter_cluster),
    color = NA,
    alpha = 0.45
  ) +
  geom_rect(
    data = block_df,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = cluster),
    color = "white",
    linewidth = 0.35
  ) +
  geom_text(
    data = label_df,
    aes(x = x, y = label_y, label = label, hjust = hjust),
    size = 3.0,
    color = "#1F2933"
  ) +
  geom_text(
    data = axis_df,
    aes(x = x, y = y, label = label),
    size = 3.6,
    fontface = "bold",
    color = "#1F2933"
  ) +
  scale_fill_manual(values = block_cols, guide = "none") +
  coord_cartesian(
    xlim = c(0.02, 0.98),
    ylim = c(min(block_df$ymin) - 12, max(block_df$ymax) + 42),
    expand = FALSE
  ) +
  theme_void(base_size = 11) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(10, 16, 10, 16)
  )

png_path <- file.path(out_dir, "pfibox_french_to_destress_cluster_alluvial.png")
pdf_path <- file.path(out_dir, "pfibox_french_to_destress_cluster_alluvial.pdf")
tsv_path <- file.path(out_dir, "pfibox_french_to_destress_cluster_alluvial_counts.tsv")

ggsave(png_path, p, width = 8.2, height = 7.2, dpi = 320)
ggsave(pdf_path, p, width = 8.2, height = 7.2, device = grDevices::pdf)
write_tsv(flow_counts, tsv_path)

message("Wrote alluvial plot:")
message("  ", png_path)
message("  ", pdf_path)
message("Wrote flow counts:")
message("  ", tsv_path)
