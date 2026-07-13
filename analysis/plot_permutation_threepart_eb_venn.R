#!/usr/bin/env Rscript

source(file.path("analysis", "_helpers.R"))

suppressPackageStartupMessages({
  library(ggplot2)
})

out_dir <- file.path(getwd(), "analysis", "outputs", "empirical_replicate_pvalues", "venn")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

perm_file <- file.path(
  getwd(),
  "analysis",
  "outputs",
  "empirical_replicate_pvalues",
  "promoter_compound_empirical_replicate_pvalues.tsv"
)
three_file <- file.path(
  getwd(),
  "analysis",
  "outputs",
  "three_part_mixture",
  "promoter_compound_three_part_mixture_results.tsv"
)

perm <- read.delim(perm_file, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
three <- read.delim(three_file, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)

perm$pair_id <- paste(perm$promoter, perm$srn_code, sep = "__")
three$pair_id <- paste(three$promoter, three$srn_code, sep = "__")

base <- merge(
  perm[, c(
    "promoter",
    "srn_code",
    "compound_label",
    "pair_id",
    "permutation_pvalue",
    "permutation_padj_by_promoter",
    "estimated_alpha_eb_padj_by_promoter"
  )],
  three[, c(
    "pair_id",
    "posterior_class",
    "local_fdr",
    "local_fdr_qvalue_by_promoter"
  )],
  by = "pair_id",
  all.x = TRUE,
  sort = FALSE
)
base$eb_bh <- is.finite(base$estimated_alpha_eb_padj_by_promoter) &
  base$estimated_alpha_eb_padj_by_promoter < 0.05
base$permutation_bh <- is.finite(base$permutation_padj_by_promoter) &
  base$permutation_padj_by_promoter < 0.05
base$three_q05 <- is.finite(base$local_fdr_qvalue_by_promoter) &
  base$local_fdr_qvalue_by_promoter <= 0.05 &
  base$posterior_class %in% c("repressed", "activated")
base$three_q20 <- is.finite(base$local_fdr_qvalue_by_promoter) &
  base$local_fdr_qvalue_by_promoter <= 0.20 &
  base$posterior_class %in% c("repressed", "activated")

circle_df <- function(x0, y0, r, group, n = 240) {
  theta <- seq(0, 2 * pi, length.out = n)
  data.frame(
    x = x0 + r * cos(theta),
    y = y0 + r * sin(theta),
    group = group
  )
}

region_label_positions <- function(spec, resolution = 550) {
  grid <- expand.grid(
    x = seq(0.04, 0.96, length.out = resolution),
    y = seq(-0.08, 0.96, length.out = resolution)
  )
  inside_a <- (grid$x - spec$x0[1])^2 + (grid$y - spec$y0[1])^2 <= spec$r[1]^2
  inside_b <- (grid$x - spec$x0[2])^2 + (grid$y - spec$y0[2])^2 <= spec$r[2]^2
  inside_c <- (grid$x - spec$x0[3])^2 + (grid$y - spec$y0[3])^2 <= spec$r[3]^2
  masks <- list(
    A = inside_a & !inside_b & !inside_c,
    B = !inside_a & inside_b & !inside_c,
    C = !inside_a & !inside_b & inside_c,
    AB = inside_a & inside_b & !inside_c,
    AC = inside_a & !inside_b & inside_c,
    BC = !inside_a & inside_b & inside_c,
    ABC = inside_a & inside_b & inside_c
  )
  do.call(
    rbind,
    lapply(names(masks), function(code) {
      region_grid <- grid[masks[[code]], , drop = FALSE]
      if (nrow(region_grid) == 0) {
        return(data.frame(code = code, x = NA_real_, y = NA_real_))
      }
      data.frame(
        code = code,
        x = stats::median(region_grid$x),
        y = stats::median(region_grid$y)
      )
    })
  )
}

make_region <- function(d, a, b, c) {
  region <- paste0(
    ifelse(d[[a]], "A", ""),
    ifelse(d[[b]], "B", ""),
    ifelse(d[[c]], "C", "")
  )
  region[region == ""] <- "None"
  region
}

plot_venn <- function(membership, c_col, c_label, prefix) {
  membership$region_code <- make_region(membership, "eb_bh", "permutation_bh", c_col)
  membership$region <- factor(
    membership$region_code,
    levels = c("None", "A", "B", "C", "AB", "AC", "BC", "ABC"),
    labels = c(
      "Not rejected",
      "EB only",
      "Permutation only",
      paste0(c_label, " only"),
      "EB + permutation",
      paste0("EB + ", c_label),
      paste0("Permutation + ", c_label),
      "All three"
    )
  )
  write.table(
    membership,
    file.path(out_dir, paste0(prefix, "_membership.tsv")),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  write.table(
    membership[membership$region_code != "None", , drop = FALSE],
    file.path(out_dir, paste0(prefix, "_rejected_only.tsv")),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  counts <- as.data.frame(table(membership$region), stringsAsFactors = FALSE)
  names(counts) <- c("region", "n")
  counts <- counts[order(-counts$n, counts$region), ]
  write.table(
    counts,
    file.path(out_dir, paste0(prefix, "_region_counts.tsv")),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  code_counts <- as.integer(table(factor(
    membership$region_code,
    levels = c("A", "B", "C", "AB", "AC", "BC", "ABC")
  )))
  names(code_counts) <- c("A", "B", "C", "AB", "AC", "BC", "ABC")
  circle_spec <- data.frame(
    group = c("DStressR EB BH", "Permutation BH", c_label),
    x0 = c(0.40, 0.60, 0.50),
    y0 = c(0.57, 0.57, 0.34),
    r = c(0.30, 0.30, 0.30)
  )
  label_positions <- region_label_positions(circle_spec)
  labels <- data.frame(
    code = names(code_counts),
    n = as.integer(code_counts),
    stringsAsFactors = FALSE
  )
  labels <- merge(
    labels,
    label_positions,
    by = "code",
    all.x = TRUE,
    sort = FALSE
  )
  circles <- rbind(
    circle_df(circle_spec$x0[1], circle_spec$y0[1], circle_spec$r[1], circle_spec$group[1]),
    circle_df(circle_spec$x0[2], circle_spec$y0[2], circle_spec$r[2], circle_spec$group[2]),
    circle_df(circle_spec$x0[3], circle_spec$y0[3], circle_spec$r[3], circle_spec$group[3])
  )
  totals <- c(
    `DStressR EB BH` = sum(membership$eb_bh, na.rm = TRUE),
    `Permutation BH` = sum(membership$permutation_bh, na.rm = TRUE),
    stats::setNames(sum(membership[[c_col]], na.rm = TRUE), c_label)
  )
  set_labels <- data.frame(
    label = paste0(names(totals), "\n", totals),
    x = c(0.22, 0.78, 0.50),
    y = c(0.93, 0.93, -0.04)
  )

  p <- ggplot() +
    geom_polygon(
      data = circles,
      aes(x, y, group = group, fill = group),
      color = "#334155",
      alpha = 0.30,
      linewidth = 0.55
    ) +
    geom_text(data = labels, aes(x, y, label = n), size = 4.2, fontface = "bold") +
    geom_text(data = set_labels, aes(x, y, label = label), size = 3.7, lineheight = 0.95) +
    scale_fill_manual(values = c(
      `DStressR EB BH` = "#16a34a",
      `Permutation BH` = "#009E73",
      stats::setNames("#dc2626", c_label)
    )) +
    coord_equal(xlim = c(0.02, 0.98), ylim = c(-0.10, 1.00), expand = FALSE) +
    theme_void() +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    labs(
      title = "BH/significance overlap of promoter-compound pairs",
      subtitle = paste0("EB BH vs B=1000 permutation BH vs ", c_label)
    )
  ggsave(file.path(out_dir, paste0(prefix, "_venn.png")), p, width = 7.5, height = 6.2, dpi = 300, bg = "white")
  ggsave(file.path(out_dir, paste0(prefix, "_venn.pdf")), p, width = 7.5, height = 6.2, bg = "white")
  counts
}

counts_q05 <- plot_venn(base, "three_q05", "Three-part q<=0.05", "eb_permutation_threepart_q05")
counts_q20 <- plot_venn(base, "three_q20", "Three-part q<=0.20", "eb_permutation_threepart_q20")

message("Wrote EB/permutation/three-part Venn outputs to: ", out_dir)
print(counts_q05, row.names = FALSE)
print(counts_q20, row.names = FALSE)
