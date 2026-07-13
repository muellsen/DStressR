test_that("prepare_assay reproduces log2 lux over growth", {
  dat <- data.frame(
    promoter = c("P1", "P1"),
    compound = c("DMSO", "C1"),
    lux = c(16, 32),
    growth = c(2, 2)
  )
  assay <- prepare_assay(dat, promoter = "promoter", compound = "compound",
                         lux = "lux", growth = "growth", growth_exponent = 1)
  expect_equal(assay$.response, c(3, 4), tolerance = 1e-6)
})

test_that("estimate_growth_exponents recovers promoter-specific scaling", {
  dat <- expand.grid(
    promoter = c("P1", "P2"),
    growth = c(1, 2, 4, 8, 16),
    replicate = seq_len(3),
    stringsAsFactors = FALSE
  )
  dat$compound <- "DMSO"
  dat$lux <- ifelse(dat$promoter == "P1", 8 * dat$growth^1, 4 * dat$growth^0.5)
  est <- estimate_growth_exponents(dat, promoter = "promoter", compound = "compound",
                                   lux = "lux", growth = "growth", min_control_n = 5,
                                   shrink = FALSE)
  expect_equal(est$alpha_raw[match("P1", est$promoter)], 1, tolerance = 1e-6)
  expect_equal(est$alpha_raw[match("P2", est$promoter)], 0.5, tolerance = 1e-6)
})

test_that("growth exponent estimation adjusts technical covariates", {
  dat <- expand.grid(
    promoter = "P1",
    plate = c("A", "B"),
    replicate = seq_len(12),
    stringsAsFactors = FALSE
  )
  dat$compound <- "DMSO"
  dat$growth <- ifelse(dat$plate == "A", 1, 8) * rep(c(1, 1.2, 1.4, 1.6), length.out = nrow(dat))
  plate_effect <- ifelse(dat$plate == "A", 0.25, 8)
  dat$lux <- 4 * dat$growth^0.5 * plate_effect

  unadjusted <- estimate_growth_exponents(
    dat,
    promoter = "promoter",
    compound = "compound",
    lux = "lux",
    growth = "growth",
    min_control_n = 8,
    shrink = FALSE
  )
  adjusted <- estimate_growth_exponents(
    dat,
    promoter = "promoter",
    compound = "compound",
    lux = "lux",
    growth = "growth",
    covariates = "plate",
    min_control_n = 8,
    shrink = FALSE
  )

  expect_gt(abs(unadjusted$alpha_raw - 0.5), 0.5)
  expect_equal(adjusted$alpha_raw, 0.5, tolerance = 1e-6)
  expect_equal(adjusted$alpha_covariates, "plate")
})

test_that("fit_destress detects simulated specific effects", {
  dat <- simulate_screen(seed = 1, n_promoters = 8, n_compounds = 12, n_replicates = 3)
  assay <- prepare_assay(dat, promoter = "promoter", compound = "compound",
                         lux = "LUX.AUC_16", growth = "od_16h.measured",
                         batch = "batch", replicate = "replicate")
  fit <- fit_destress(assay, technical = c("batch", "replicate"))
  res <- results(fit)
  expect_true(all(c("specific_effect", "global_effect", "total_effect") %in% names(res)))
  truth <- unique(dat[dat$compound != "DMSO", c("promoter", "compound", "truth_specific")])
  joined <- merge(res, truth, by = c("promoter", "compound"))
  active <- abs(joined$truth_specific) > 0
  expect_gt(stats::cor(joined$specific_effect[active], joined$truth_specific[active]), 0.7)
})

test_that("call_hits adds interpretable classes", {
  tab <- data.frame(
    specific_effect = c(1, -1, 0.1),
    specific_padj = c(0.01, 0.01, 0.8)
  )
  out <- call_hits(tab)
  expect_equal(out$hit, c("Upregulated", "Downregulated", "Not DE"))
})

test_that("fit_effect_mixture separates three effect classes", {
  set.seed(2)
  tab <- data.frame(
    promoter = "P1",
    compound = paste0("C", seq_len(180)),
    truth = rep(c("repressed", "null", "activated"), c(35, 110, 35)),
    stringsAsFactors = FALSE
  )
  tab$specific_effect <- c(
    stats::rt(35, df = 5) * 0.08 - 0.75,
    stats::rt(110, df = 5) * 0.08,
    stats::rt(35, df = 5) * 0.08 + 0.75
  )

  out <- fit_effect_mixture(tab, df = 5)
  summary <- attr(out, "mixture_summary")

  expect_true(all(c(
    "prob_repressed",
    "prob_null",
    "prob_activated",
    "local_fdr",
    "posterior_nonnull",
    "local_fdr_qvalue_by_promoter"
  ) %in% names(out)))
  expect_equal(nrow(summary), 1)
  expect_lt(summary$location_repressed, summary$location_null)
  expect_lt(summary$location_null, summary$location_activated)
  expect_gt(mean(out$posterior_class == out$truth), 0.85)
  expect_true(all(out$empirical_null_padj_by_promoter >= 0 & out$empirical_null_padj_by_promoter <= 1))
  expect_true(all(out$local_fdr_qvalue_by_promoter >= 0 & out$local_fdr_qvalue_by_promoter <= 1))
})

test_that("empirical_replicate_pvalues compares replicate averages to matched controls", {
  dat <- expand.grid(
    promoter = "P1",
    libplate = "lp1",
    replicate = c("r1", "r2"),
    compound = c(paste0("DMSO", seq_len(6)), "C_high", "C_mid"),
    stringsAsFactors = FALSE
  )
  control_means <- c(-0.2, -0.1, -0.05, 0.05, 0.1, 0.2)
  dat$value <- 0
  for (i in seq_along(control_means)) {
    dat$value[dat$compound == paste0("DMSO", i)] <- control_means[i] +
      ifelse(dat$replicate[dat$compound == paste0("DMSO", i)] == "r1", -0.01, 0.01)
  }
  dat$value[dat$compound == "C_high"] <- c(0.95, 1.05)
  dat$value[dat$compound == "C_mid"] <- c(0.04, 0.06)

  out <- empirical_replicate_pvalues(
    dat,
    value = "value",
    promoter = "promoter",
    compound = "compound",
    control = paste0("DMSO", seq_len(6)),
    replicate = "replicate",
    strata = "libplate",
    min_replicates = 2,
    min_null = 5,
    permutation = TRUE,
    B = 200,
    seed = 1
  )

  expect_equal(nrow(out), 2)
  expect_lt(out$empirical_pvalue[out$compound == "C_high"], out$empirical_pvalue[out$compound == "C_mid"])
  expect_lt(out$permutation_pvalue[out$compound == "C_high"], out$permutation_pvalue[out$compound == "C_mid"])
  expect_equal(out$n_replicates[out$compound == "C_high"], 2)
  expect_equal(out$null_n[out$compound == "C_high"], 6)
  expect_true(all(out$permutation_pvalue >= 1 / 201 & out$permutation_pvalue <= 1))
})

test_that("add_dgrowthr_growth joins DGrowthR growth parameters", {
  if (!methods::isClass("DGrowthR")) {
    methods::setClass(
      "DGrowthR",
      slots = list(metadata = "data.frame", growth_parameters = "data.frame")
    )
  }
  object <- methods::new(
    "DGrowthR",
    metadata = data.frame(
      curve_id = c("c1", "c2"),
      strain_plate = c("p1", "p2"),
      stringsAsFactors = FALSE
    ),
    growth_parameters = data.frame(
      gpfit_id = c("p1", "p2"),
      OD_16 = c(0.31, 0.52),
      AUC = c(4.1, 5.2),
      stringsAsFactors = FALSE
    )
  )
  assay <- data.frame(
    curve_id = c("c2", "c1"),
    lux = c(20, 10),
    stringsAsFactors = FALSE
  )

  joined <- add_dgrowthr_growth(
    assay,
    object,
    by = "curve_id",
    model_covariate = "strain_plate",
    growth_metric = "OD_16",
    output = "dgrowthr_od16"
  )

  expect_equal(joined$curve_id, assay$curve_id)
  expect_equal(joined$dgrowthr_od16, c(0.52, 0.31))
})

test_that("plot_volcano returns a ggplot object", {
  skip_if_not_installed("ggplot2")
  tab <- data.frame(
    promoter = c("P1", "P1", "P2", "P3"),
    compound = c("C1", "C2", "C1", "C3"),
    compound_name = c("Drug A", "Drug B", "Drug A", "Drug C"),
    specific_effect = c(2.1, -0.3, -1.8, 0.5),
    specific_padj = c(0.001, 0.7, 0.02, 0.4),
    stringsAsFactors = FALSE
  )

  p <- plot_volcano(tab, compound_label = "compound_name", top_n = 2, top_promoters = 2)

  expect_s3_class(p, "ggplot")
})

test_that("plot_response_heatmap returns a ggplot with matrix attribute", {
  skip_if_not_installed("ggplot2")
  tab <- expand.grid(
    promoter = c("P1", "P2"),
    compound = c("C1", "C2", "C3"),
    stringsAsFactors = FALSE
  )
  tab$compound_name <- c("Drug A", "Drug A", "Drug B", "Drug B", "Drug C", "Drug C")
  tab$specific_effect <- c(1, -1, 0.2, -0.2, 2, -2)

  p <- plot_response_heatmap(
    tab,
    compound_label = "compound_name",
    top_n_compounds = Inf,
    cluster_rows = FALSE,
    cluster_cols = FALSE
  )
  mat <- attr(p, "response_matrix")

  expect_s3_class(p, "ggplot")
  expect_equal(dim(mat), c(2, 3))
  expect_equal(rownames(mat), c("P1", "P2"))
})

test_that("plot_effect_histogram returns pooled and promoter-faceted plots", {
  skip_if_not_installed("ggplot2")
  tab <- data.frame(
    promoter = rep(c("P1", "P2"), each = 6),
    specific_effect = c(-1, -0.5, -0.1, 0, 0.3, 1.2, -0.2, 0, 0.1, 0.4, 0.9, 1.4)
  )

  pooled <- plot_effect_histogram(tab, bins = 10)
  per_promoter <- plot_effect_histogram(tab, by = "promoter", bins = 10)

  expect_s3_class(pooled, "ggplot")
  expect_s3_class(per_promoter, "ggplot")
})

test_that("plot_response_cluster_blocks returns cluster summaries", {
  skip_if_not_installed("ggplot2")
  tab <- expand.grid(
    promoter = c("P1", "P2", "P3", "P4"),
    compound = c("C1", "C2", "C3", "C4", "C5"),
    stringsAsFactors = FALSE
  )
  tab$compound_name <- paste("Drug", tab$compound)
  tab$specific_effect <- c(
    1.2, 1.1, -0.1, -0.2,
    1.0, 0.9, -0.2, -0.1,
    -0.2, -0.1, 1.3, 1.1,
    -0.1, -0.2, 1.1, 1.2,
    0.1, 0.2, 0.0, -0.1
  )

  p <- plot_response_cluster_blocks(
    tab,
    compound_label = "compound_name",
    n_promoter_clusters = 2,
    n_compound_clusters = 2
  )

  expect_s3_class(p, "ggplot")
  expect_equal(dim(attr(p, "response_matrix")), c(4, 5))
  expect_equal(nrow(attr(p, "promoter_clusters")), 4)
  expect_equal(nrow(attr(p, "compound_clusters")), 5)
  expect_equal(nrow(attr(p, "block_summary")), 4)
})

test_that("plot_response_clustered_heatmap writes clustered heatmap output", {
  tab <- expand.grid(
    promoter = c("P1", "P2", "P3", "P4"),
    compound = c("C1", "C2", "C3", "C4", "C5"),
    stringsAsFactors = FALSE
  )
  tab$compound_name <- paste("Drug", tab$compound)
  tab$specific_effect <- c(
    1.2, 1.1, -0.1, -0.2,
    1.0, 0.9, -0.2, -0.1,
    -0.2, -0.1, 1.3, 1.1,
    -0.1, -0.2, 1.1, 1.2,
    0.1, 0.2, 0.0, -0.1
  )
  out_file <- tempfile(fileext = ".png")

  out <- plot_response_clustered_heatmap(
    tab,
    compound_label = "compound_name",
    n_promoter_clusters = 2,
    n_compound_clusters = 2,
    file = out_file,
    width = 8,
    height = 6,
    show_colnames = FALSE
  )

  expect_true(file.exists(out_file))
  expect_equal(dim(out$response_matrix), c(4, 5))
  expect_equal(nrow(out$promoter_clusters), 4)
  expect_equal(nrow(out$compound_clusters), 5)
})
