#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

out_dir <- file.path(getwd(), "analysis", "outputs", "eb_moderated_variance")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

variance_file <- file.path(out_dir, "promoter_empirical_bayes_variances.tsv")

if (!file.exists(variance_file)) {
  stop(
    "Missing promoter variance table: ", variance_file,
    "\nRun analysis/apply_eb_moderated_variances.R first.",
    call. = FALSE
  )
}

variance_df <- read.delim(variance_file, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
variance_df <- variance_df[order(-variance_df$promoter_s2), ]
variance_df$promoter <- factor(variance_df$promoter, levels = variance_df$promoter)

plot_df <- rbind(
  data.frame(
    promoter = variance_df$promoter,
    variance = variance_df$promoter_s2,
    sd = sqrt(variance_df$promoter_s2),
    estimate = "Raw variance",
    stringsAsFactors = FALSE
  ),
  data.frame(
    promoter = variance_df$promoter,
    variance = variance_df$moderated_s2,
    sd = sqrt(variance_df$moderated_s2),
    estimate = "EB-moderated variance",
    stringsAsFactors = FALSE
  )
)
plot_df$estimate <- factor(plot_df$estimate, levels = c("Raw variance", "EB-moderated variance"))

prior_s2 <- unique(variance_df$eb_prior_s2)
prior_df <- unique(variance_df$eb_prior_df)

p_var <- ggplot(plot_df, aes(promoter, variance, color = estimate, group = estimate)) +
  geom_hline(
    yintercept = prior_s2[1],
    linetype = "dashed",
    color = "#525252",
    linewidth = 0.45
  ) +
  geom_line(linewidth = 0.45, alpha = 0.75) +
  geom_point(size = 2.1, alpha = 0.95) +
  scale_color_manual(values = c("Raw variance" = "#64748b", "EB-moderated variance" = "#16a34a")) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
    legend.position = "bottom"
  ) +
  labs(
    title = "Empirical-Bayes moderation of promoter-specific variances",
    subtitle = paste0(
      "Dashed line: prior variance = ", signif(prior_s2[1], 3),
      " (prior df = ", round(prior_df[1], 2), ")"
    ),
    x = "Promoter, ordered high-to-low by raw DMSO residual variance",
    y = "DMSO residual variance",
    color = NULL
  )

ggsave(
  file.path(out_dir, "eb_variance_moderation_by_promoter.png"),
  p_var,
  width = 10,
  height = 5.5,
  dpi = 300,
  bg = "white"
)
ggsave(
  file.path(out_dir, "eb_variance_moderation_by_promoter.pdf"),
  p_var,
  width = 10,
  height = 5.5,
  bg = "white"
)

p_sd <- ggplot(plot_df, aes(promoter, sd, color = estimate, group = estimate)) +
  geom_hline(
    yintercept = sqrt(prior_s2[1]),
    linetype = "dashed",
    color = "#525252",
    linewidth = 0.45
  ) +
  geom_line(linewidth = 0.45, alpha = 0.75) +
  geom_point(size = 2.1, alpha = 0.95) +
  scale_color_manual(values = c("Raw variance" = "#64748b", "EB-moderated variance" = "#16a34a")) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
    legend.position = "bottom"
  ) +
  labs(
    title = "Empirical-Bayes moderation of promoter-specific residual SDs",
    subtitle = paste0(
      "Dashed line: prior SD = ", signif(sqrt(prior_s2[1]), 3),
      " (prior df = ", round(prior_df[1], 2), ")"
    ),
    x = "Promoter, ordered high-to-low by raw DMSO residual variance",
    y = "DMSO residual standard deviation",
    color = NULL
  )

ggsave(
  file.path(out_dir, "eb_sd_moderation_by_promoter.png"),
  p_sd,
  width = 10,
  height = 5.5,
  dpi = 300,
  bg = "white"
)
ggsave(
  file.path(out_dir, "eb_sd_moderation_by_promoter.pdf"),
  p_sd,
  width = 10,
  height = 5.5,
  bg = "white"
)

write.table(
  plot_df,
  file.path(out_dir, "eb_variance_moderation_by_promoter_plot_data.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

message("Wrote EB variance diagnostic plots to: ", out_dir)

replicate_file <- file.path(out_dir, "workflow_vs_destress_eb_replicate_pvalues.tsv")
if (file.exists(replicate_file)) {
  replicate_df <- read.delim(replicate_file, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
  dmso_codes <- unique(replicate_df$srn_code[grepl("_A1$|_A7$|_A24$|_B12$|_C10$|_D21$|_E3$|_F18$|_G6$|_H14$|_J17$|_K20$|_L11$|_M19$|_O9$|_O24$|_P1$|_P24$", replicate_df$srn_code)])

  # More robustly identify DMSO rows by using the fitted group DMSO mean:
  # rows from DMSO wells have finite group means and are retained in the full
  # replicate p-value table. The old LibMap is not needed here because each
  # promoter-libplate-replicate DMSO residual cloud was already centered.
  dmso_replicate <- replicate_df[replicate_df$srn_code %in% dmso_codes, ]
  if (nrow(dmso_replicate) > 0) {
    dmso_replicate$promoter_replicate <- paste(dmso_replicate$promoter, dmso_replicate$replicate, sep = "_")
    dmso_replicate$centered_specific <- dmso_replicate$destress_specific_effect -
      dmso_replicate$eb_dmso_group_mean

    replicate_var <- do.call(
      rbind,
      lapply(split(dmso_replicate, dmso_replicate$promoter_replicate), function(d) {
        data.frame(
          promoter = d$promoter[1],
          replicate = d$replicate[1],
          promoter_replicate = d$promoter_replicate[1],
          raw_replicate_s2 = stats::var(d$centered_specific, na.rm = TRUE),
          n_dmso = sum(is.finite(d$centered_specific)),
          stringsAsFactors = FALSE
        )
      })
    )

    replicate_var <- merge(
      replicate_var,
      variance_df[, c("promoter", "promoter_s2", "moderated_s2")],
      by = "promoter",
      all.x = TRUE,
      sort = FALSE
    )

    promoter_order <- as.character(variance_df$promoter)
    replicate_var$promoter <- factor(replicate_var$promoter, levels = promoter_order)
    replicate_var <- replicate_var[order(replicate_var$promoter, replicate_var$replicate), ]
    replicate_var$promoter_replicate <- factor(
      replicate_var$promoter_replicate,
      levels = replicate_var$promoter_replicate
    )

    replicate_plot_df <- rbind(
      data.frame(
        promoter_replicate = replicate_var$promoter_replicate,
        variance = replicate_var$raw_replicate_s2,
        estimate = "Raw replicate variance",
        stringsAsFactors = FALSE
      ),
      data.frame(
        promoter_replicate = replicate_var$promoter_replicate,
        variance = replicate_var$moderated_s2,
        estimate = "EB-moderated promoter variance",
        stringsAsFactors = FALSE
      )
    )
    replicate_plot_df$estimate <- factor(
      replicate_plot_df$estimate,
      levels = c("Raw replicate variance", "EB-moderated promoter variance")
    )

    p_rep <- ggplot(replicate_plot_df, aes(promoter_replicate, variance, color = estimate, group = estimate)) +
      geom_hline(
        yintercept = prior_s2[1],
        linetype = "dashed",
        color = "#525252",
        linewidth = 0.45
      ) +
      geom_point(size = 1.9, alpha = 0.95) +
      scale_color_manual(values = c("Raw replicate variance" = "#64748b", "EB-moderated promoter variance" = "#16a34a")) +
      theme_bw(base_size = 10) +
      theme(
        axis.text.x = element_text(angle = 75, hjust = 1, vjust = 1, size = 7),
        legend.position = "bottom"
      ) +
      labs(
        title = "Replicate-aware empirical-Bayes variance diagnostic",
        subtitle = "Promoters ordered high-to-low by pooled raw variance; replicates kept adjacent",
        x = "Promoter replicate",
        y = "DMSO residual variance",
        color = NULL
      )

    ggsave(
      file.path(out_dir, "eb_variance_moderation_by_promoter_replicate.png"),
      p_rep,
      width = 13,
      height = 5.8,
      dpi = 300,
      bg = "white"
    )
    ggsave(
      file.path(out_dir, "eb_variance_moderation_by_promoter_replicate.pdf"),
      p_rep,
      width = 13,
      height = 5.8,
      bg = "white"
    )

    write.table(
      replicate_var,
      file.path(out_dir, "eb_variance_moderation_by_promoter_replicate_plot_data.tsv"),
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
  }
}
