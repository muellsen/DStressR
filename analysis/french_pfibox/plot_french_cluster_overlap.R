source(file.path("analysis", "_helpers.R"))
source(file.path("analysis", "_paper_artifacts.R"))
ggplot2 <- asNamespace("ggplot2")
out_dir <- analysis_output_dir("french_pfibox", "promoter_activity_analysis")
run_tag <- Sys.getenv("PFIBOX_RUN_TAG", "reporters2000_perturbations100_T2-18h_treatment_replicate_scalable")
suffix <- Sys.getenv("PFIBOX_CLUSTER_OUTPUT_SUFFIX", "total_effect_ward")
suffix <- if(nzchar(suffix)) paste0("_",suffix) else ""
artifact <- read_paper_artifact(file.path(out_dir,paste0("pfibox_paper_mapping_",run_tag,suffix,".rds")))
long_counts <- artifact$diagnostics$long_counts
p <- ggplot2$ggplot(
  long_counts,
  ggplot2$aes(french_major_cluster, dstressr_reporter_cluster, fill = proportion)
) +
  ggplot2$geom_tile(color = "white", linewidth = 0.35) +
  ggplot2$geom_text(
    ggplot2$aes(label = ifelse(n > 0, n, "")),
    size = 2.6,
    color = "#111827"
  ) +
  ggplot2$scale_fill_gradient(
    low = "white",
    high = "#3f6fb5",
    labels = scales$percent_format(accuracy = 1),
    name = "Within-cluster\nfraction"
  ) +
  ggplot2$theme_light(base_size = 10) +
  ggplot2$theme(
    panel.grid = ggplot2$element_blank(),
    axis.text = ggplot2$element_text(color = "#4B5563"),
    axis.title = ggplot2$element_text(color = "#111827"),
    legend.position = "right",
    legend.title = ggplot2$element_text(size = 8),
    legend.text = ggplot2$element_text(size = 7)
  ) +
  ggplot2$labs(
    x = "French Mfuzz time-course cluster",
    y = "DStressR total-effect reporter cluster"
  )

base <- file.path(out_dir, "pfibox_destress_vs_french_mfuzz_cluster_overlap")
ggplot2$ggsave(paste0(base, ".png"), p, width = 7.2, height = 4.8, dpi = 300, bg = "white")
ggplot2$ggsave(paste0(base, ".pdf"), p, width = 7.2, height = 4.8, bg = "white")

