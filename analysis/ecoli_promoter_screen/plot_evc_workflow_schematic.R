source(file.path("analysis", "_helpers.R"))

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package `ggplot2` is required for the EVC workflow schematic.", call. = FALSE)
}

ggplot2 <- asNamespace("ggplot2")

out_dir <- analysis_output_dir("binsfeld_modeling_steps")
map_file <- file.path(
  "assets",
  "brochadolab-Binsfeld2025-e5ada46",
  "Code",
  "ReporterScreen_Binsfeld25",
  "R-scripts",
  "Map.txt"
)

plate_map <- utils::read.table(map_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
plate_map$row <- substr(plate_map$Well, 1, 1)
plate_map$col <- as.integer(sub("^[A-P]", "", plate_map$Well))
plate_map$row_index <- match(plate_map$row, LETTERS[1:16])
plate_map$x <- plate_map$col
plate_map$y <- 17 - plate_map$row_index
plate_map$well_class <- "Other compounds"
plate_map$well_class[grepl("Water", plate_map$Drug)] <- "Water controls"
plate_map$well_class[plate_map$Drug == "Clarithromycin"] <- "Example compound"

selected_well <- "F13"
selected <- plate_map[plate_map$Well == selected_well, , drop = FALSE]
selected_text <- paste0(selected$Well, ": ", selected$Drug, ", dose ", selected$ConcMock)

plate_map$xmin <- 1.0 + plate_map$x - 0.46
plate_map$xmax <- 1.0 + plate_map$x + 0.46
plate_map$ymin <- 27.5 + plate_map$y - 0.46
plate_map$ymax <- 27.5 + plate_map$y + 0.46
plate_map$selected <- plate_map$Well == selected_well

matrix_x <- c(39, 49, 59, 69, 79)
matrix_labels <- c("well h", "reporter q", "EVC 1", "EVC 2", "row mean")
matrix_header <- data.frame(
  xmin = matrix_x - 4.4,
  xmax = matrix_x + 4.4,
  ymin = 42.0,
  ymax = 45.2,
  label = matrix_labels,
  stringsAsFactors = FALSE
)
matrix_rows <- data.frame(
  row_name = c("A1", "A2", selected_well, "I13", "P12"),
  y = c(38.4, 35.0, 31.6, 28.2, 24.8),
  compound = c("Water_1", "Water_1", "Clarithromycin", "Water_2", "Tetracycline"),
  reporter = c("r[A1,q]", "r[A2,q]", "r[F13,q]", "r[I13,q]", "r[P12,q]"),
  evc1 = c("r[A1,E1]", "r[A2,E1]", "r[F13,E1]", "r[I13,E1]", "r[P12,E1]"),
  evc2 = c("r[A1,E2]", "r[A2,E2]", "r[F13,E2]", "r[I13,E2]", "r[P12,E2]"),
  mean = c("e[A1]", "e[A2]", "e[F13]", "e[I13]", "e[P12]"),
  stringsAsFactors = FALSE
)
matrix_cells <- do.call(rbind, lapply(seq_len(nrow(matrix_rows)), function(i) {
  data.frame(
    xmin = matrix_x - 4.4,
    xmax = matrix_x + 4.4,
    ymin = matrix_rows$y[i] - 1.4,
    ymax = matrix_rows$y[i] + 1.4,
    label = c(
      paste0(matrix_rows$row_name[i], "\n", matrix_rows$compound[i]),
      matrix_rows$reporter[i],
      matrix_rows$evc1[i],
      matrix_rows$evc2[i],
      matrix_rows$mean[i]
    ),
    selected = matrix_rows$row_name[i] == selected_well,
    stringsAsFactors = FALSE
  )
}))

binsfeld_boxes <- data.frame(
  xmin = c(91, 91),
  xmax = c(126, 126),
  ymin = c(34.0, 24.5),
  ymax = c(41.0, 31.5),
  label = c(
    "Binsfeld et al.\nFor every well position h:\ne[h] = mean(EVC 1, EVC 2)",
    "For each reporter-replicate q:\nrobust regression r[h,q] on e[h]\nscore s[i] = residual"
  ),
  stringsAsFactors = FALSE
)

dstress_boxes <- data.frame(
  xmin = c(39, 70, 101),
  xmax = c(62, 93, 126),
  ymin = c(5.2, 5.2, 5.2),
  ymax = c(14.8, 14.8, 14.8),
  label = c(
    "Reporter and growth summaries\n(l[i], g[i])",
    "DStressR response scale\ny[i] = log2(l[i]) - alpha[a] log2(g[i])",
    "Optional background calibration\nthen linear-model inference"
  ),
  stringsAsFactors = FALSE
)

arrows <- data.frame(
  x = c(25.8, 83.6, 83.6, 62.5, 93.5),
  xend = c(34.0, 90.4, 90.4, 69.5, 100.5),
  y = c(35.0, 37.5, 28.0, 10.0, 10.0),
  yend = c(35.0, 37.5, 28.0, 10.0, 10.0),
  stringsAsFactors = FALSE
)

legend_data <- data.frame(
  x = c(2.5, 13.5, 24.5),
  y = c(23.2, 23.2, 23.2),
  well_class = c("Example compound", "Water controls", "Other compounds"),
  label = c("Clarithromycin wells", "Water controls", "Other wells"),
  stringsAsFactors = FALSE
)

plot <- ggplot2$ggplot() +
  ggplot2$geom_rect(
    data = plate_map,
    ggplot2$aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = well_class),
    color = "white",
    linewidth = 0.08
  ) +
  ggplot2$geom_rect(
    data = plate_map[plate_map$selected, ],
    ggplot2$aes(xmin = xmin - 0.12, xmax = xmax + 0.12, ymin = ymin - 0.12, ymax = ymax + 0.12),
    fill = NA,
    color = "#111827",
    linewidth = 0.7
  ) +
  ggplot2$annotate("text", x = 1.0, y = 46.0, label = "Actual 384-well compound layout",
                   hjust = 0, fontface = "bold", size = 4.6, color = "#111827") +
  ggplot2$annotate("text", x = 1.0, y = 24.8,
                   label = paste0("Highlighted well position: ", selected_text),
                   hjust = 0, size = 3.35, color = "#374151") +
  ggplot2$geom_point(
    data = legend_data,
    ggplot2$aes(x = x, y = y, fill = well_class),
    shape = 22,
    size = 4.8,
    color = "white"
  ) +
  ggplot2$geom_text(
    data = legend_data,
    ggplot2$aes(x = x + 1.7, y = y, label = label),
    hjust = 0,
    size = 3.1,
    color = "#374151"
  ) +
  ggplot2$geom_rect(
    data = matrix_header,
    ggplot2$aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "#e5e7eb",
    color = "white",
    linewidth = 0.35
  ) +
  ggplot2$geom_rect(
    data = matrix_cells,
    ggplot2$aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = ifelse(matrix_cells$selected, "#f5d0fe", "#f8fafc"),
    color = "white",
    linewidth = 0.35
  ) +
  ggplot2$geom_text(
    data = matrix_header,
    ggplot2$aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = label),
    fontface = "bold",
    size = 3.25,
    color = "#111827"
  ) +
  ggplot2$geom_text(
    data = matrix_cells,
    ggplot2$aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = label),
    size = 2.85,
    lineheight = 0.92,
    color = "#111827"
  ) +
  ggplot2$annotate("text", x = 35, y = 46.0,
                   label = "Same well positions become rows in a response table",
                   hjust = 0, fontface = "bold", size = 4.6, color = "#111827") +
  ggplot2$geom_rect(
    data = binsfeld_boxes,
    ggplot2$aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "#ffffff",
    color = "#374151",
    linewidth = 0.55
  ) +
  ggplot2$geom_text(
    data = binsfeld_boxes,
    ggplot2$aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = label),
    size = 3.15,
    lineheight = 0.95,
    color = "#111827"
  ) +
  ggplot2$annotate("text", x = 91, y = 46.0,
                   label = "Published Binsfeld score construction",
                   hjust = 0, fontface = "bold", size = 4.6, color = "#0f766e") +
  ggplot2$annotate("segment", x = 1, xend = 126, y = 19.7, yend = 19.7,
                   linewidth = 0.45, color = "#d1d5db") +
  ggplot2$annotate("text", x = 1.0, y = 17.5,
                   label = "DStressR generalization",
                   hjust = 0, fontface = "bold", size = 4.8, color = "#15803d") +
  ggplot2$geom_rect(
    data = dstress_boxes,
    ggplot2$aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "#ffffff",
    color = "#374151",
    linewidth = 0.55
  ) +
  ggplot2$geom_text(
    data = dstress_boxes,
    ggplot2$aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = label),
    size = 3.1,
    lineheight = 0.95,
    color = "#111827"
  ) +
  ggplot2$geom_segment(
    data = arrows,
    ggplot2$aes(x = x, xend = xend, y = y, yend = yend),
    arrow = grid::arrow(length = grid::unit(0.018, "npc"), type = "closed"),
    linewidth = 0.55,
    color = "#111827"
  ) +
  ggplot2$scale_fill_manual(
    values = c(
      "Example compound" = "#d946ef",
      "Water controls" = "#f59e0b",
      "Other compounds" = "#cbd5e1"
    )
  ) +
  ggplot2$coord_cartesian(xlim = c(0, 128), ylim = c(3.5, 48), expand = FALSE) +
  ggplot2$theme_void(base_size = 12) +
  ggplot2$theme(
    legend.position = "none",
    plot.background = ggplot2$element_rect(fill = "white", color = NA),
    panel.background = ggplot2$element_rect(fill = "white", color = NA),
    plot.margin = ggplot2$margin(8, 8, 8, 8)
  )

ggplot2$ggsave(
  file.path(out_dir, "binsfeld_vs_dstressr_evc_workflow_schematic.png"),
  plot,
  width = 14.0,
  height = 6.4,
  dpi = 320,
  bg = "white"
)
ggplot2$ggsave(
  file.path(out_dir, "binsfeld_vs_dstressr_evc_workflow_schematic.pdf"),
  plot,
  width = 14.0,
  height = 6.4,
  bg = "white"
)

message("Wrote EVC workflow schematic to: ", out_dir)
