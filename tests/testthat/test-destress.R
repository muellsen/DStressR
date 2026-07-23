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
  expect_true(all(c("a_raw", "a_raw_se", "a_raw_df") %in% names(est)))
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

test_that("Binsfeld reporter data support DStressR model analysis", {
  data("binsfeld_reporter_auc", package = "DStressR")
  data("binsfeld_reporter_scores", package = "DStressR")

  expect_equal(nrow(binsfeld_reporter_auc), 24576)
  expect_true(all(c(
    "strain", "promoter", "compound", "od_auc", "lux_auc", "removed"
  ) %in% names(binsfeld_reporter_auc)))
  expect_equal(
    sort(unique(binsfeld_reporter_auc$compound[grepl("^Water_", binsfeld_reporter_auc$drug)])),
    "Water"
  )
  expect_true(all(c("Scores", "Z_scores") %in% unique(binsfeld_reporter_scores$statistic)))

  wt_auc <- binsfeld_reporter_auc[
    binsfeld_reporter_auc$strain == "WT" &
      binsfeld_reporter_auc$removed == "No" &
      binsfeld_reporter_auc$compound %in% c("Water", "Azithromycin", "Clarithromycin"),
  ]
  assay <- prepare_assay(
    wt_auc,
    promoter = "promoter",
    compound = "compound",
    control = "Water",
    lux = "lux_auc",
    growth = "od_auc",
    growth_exponent = "estimate",
    batch = "concentration_index",
    replicate = "replicate"
  )
  fit <- fit_destress(
    assay,
    technical = c("replicate", "concentration_index"),
    empirical_bayes = TRUE,
    adjustment = "by_promoter",
    interaction = FALSE,
    empty_vector_promoter = "EVC"
  )
  res <- results(fit)

  expect_s3_class(fit, "destress_fit")
  expect_true(all(c("empty_vector_effect", "specific_effect", "specific_padj") %in% names(res)))
  expect_false(is.null(fit$growth_exponents))
  expect_true(nrow(res) > 0)
})

test_that("fit_workflow dispatches to the model workflow", {
  dat <- simulate_screen(seed = 3, n_promoters = 5, n_compounds = 6, n_replicates = 2)
  assay <- prepare_assay(dat, promoter = "promoter", compound = "compound",
                         lux = "LUX.AUC_16", growth = "od_16h.measured",
                         batch = "batch", replicate = "replicate")

  direct <- fit_destress(assay, technical = c("batch", "replicate"))
  via_workflow <- fit_workflow(assay, workflow = "model", technical = c("batch", "replicate"))

  expect_s3_class(via_workflow, "destress_fit")
  expect_equal(attr(via_workflow, "destress_workflow"), "model")
  expect_equal(results(via_workflow), results(direct), tolerance = 1e-10)
})

test_that("fit_destress exposes staged model options", {
  dat <- simulate_screen(seed = 4, n_promoters = 5, n_compounds = 6, n_replicates = 2)
  assay <- prepare_assay(dat, promoter = "promoter", compound = "compound",
                         lux = "LUX.AUC_16", growth = "od_16h.measured",
                         batch = "batch", replicate = "replicate")

  fit <- fit_destress(
    assay,
    technical = c("batch", "replicate"),
    normalization = "model",
    testing = "student_t",
    aggregation = "none",
    adjustment = "by_promoter"
  )
  res <- results(fit)
  expected <- adjust_pvalues(res, pvalue = "specific_pvalue", output = "expected_specific_padj")

  expect_s3_class(fit, "destress_fit")
  expect_false(fit$empirical_bayes)
  expect_equal(fit$stages$normalization, "linear_model")
  expect_equal(fit$stages$testing, "student_t")
  expect_equal(fit$stages$adjustment, "by_promoter")
  expect_equal(res$specific_padj, expected$expected_specific_padj, tolerance = 1e-12)
})

test_that("fit_destress can fit scalable promoter-specific models", {
  dat <- simulate_screen(seed = 7, n_promoters = 5, n_compounds = 6, n_replicates = 3)
  assay <- prepare_assay(dat, promoter = "promoter", compound = "compound",
                         lux = "LUX.AUC_16", growth = "od_16h.measured",
                         batch = "batch", replicate = "replicate")

  fit <- fit_destress(
    assay,
    technical = c("batch", "replicate"),
    interaction = FALSE,
    empirical_bayes = FALSE,
    adjustment = "by_promoter"
  )
  res <- results(fit)

  expect_s3_class(fit, "destress_fit")
  expect_false(fit$interaction)
  expect_null(fit$full_fit)
  expect_equal(length(unique(fit$promoter_effects$promoter)), length(unique(dat$promoter)))
  expect_equal(nrow(res), length(unique(dat$promoter)) * length(setdiff(unique(dat$compound), "DMSO")))
  expect_true(all(c("specific_effect", "specific_pvalue", "specific_padj") %in% names(res)))
  expect_true(all(c("total_effect", "total_pvalue", "total_padj") %in% names(res)))
  expect_true(all(is.finite(res$total_pvalue)))
  expect_gt(sum(abs(res$specific_effect) > 1e-8), 0)
  expect_true(all(is.finite(res$specific_pvalue)))
  expect_true(all(res$specific_padj >= 0 & res$specific_padj <= 1))

  truth <- unique(dat[dat$compound != "DMSO", c("promoter", "compound", "truth_specific")])
  joined <- merge(res, truth, by = c("promoter", "compound"))
  active <- abs(joined$truth_specific) > 0
  expect_gt(stats::cor(joined$specific_effect[active], joined$truth_specific[active]), 0.7)
})

test_that("fit_destress separates global compound effects from promoter-specific effects", {
  dat <- expand.grid(
    promoter = paste0("P", seq_len(6)),
    compound = c("DMSO", "C_global", "C_specific"),
    replicate = paste0("r", seq_len(6)),
    stringsAsFactors = FALSE
  )
  baseline <- stats::setNames(seq(9.5, 10.5, length.out = 6), paste0("P", seq_len(6)))
  dat$value <- baseline[dat$promoter] +
    ifelse(dat$compound == "C_global", 1.5, 0) +
    ifelse(dat$compound == "C_specific" & dat$promoter == "P3", 1.5, 0)

  assay <- prepare_assay(
    dat,
    promoter = "promoter",
    compound = "compound",
    control = "DMSO",
    response = "value",
    replicate = "replicate"
  )
  fit <- fit_destress(assay, technical = "replicate", empirical_bayes = FALSE)
  res <- results(fit)

  global_rows <- res[res$compound == "C_global", ]
  expect_equal(global_rows$total_effect, rep(1.5, nrow(global_rows)), tolerance = 1e-8)
  expect_equal(unique(global_rows$global_effect), 1.5, tolerance = 1e-8)
  expect_equal(global_rows$specific_effect, rep(0, nrow(global_rows)), tolerance = 1e-8)
  expect_true(all(global_rows$specific_pvalue > 0.9))
  expect_true(all(global_rows$global_pvalue < 1e-8))

  specific_row <- res[res$compound == "C_specific" & res$promoter == "P3", ]
  expect_gt(specific_row$specific_effect, 1)
  expect_lt(specific_row$specific_pvalue, 1e-8)
})

test_that("fit_destress can remove a low-rank compound background", {
  promoters <- paste0("P", seq_len(6))
  compounds <- c("DMSO", "C_factor1", "C_factor2", "C_specific")
  dat <- expand.grid(
    promoter = promoters,
    compound = compounds,
    replicate = paste0("r", seq_len(5)),
    stringsAsFactors = FALSE
  )

  baseline <- stats::setNames(seq(9.5, 10.5, length.out = length(promoters)), promoters)
  loading <- stats::setNames(c(-2, -1, 0, 0, 1, 2), promoters)
  score <- c(DMSO = 0, C_factor1 = 1.2, C_factor2 = -0.8, C_specific = 0)
  sparse_specific <- stats::setNames(c(1, -2, 1, 1, -2, 1), promoters)
  dat$value <- baseline[dat$promoter] +
    loading[dat$promoter] * score[dat$compound] +
    ifelse(dat$compound == "C_specific", sparse_specific[dat$promoter], 0)

  assay <- prepare_assay(
    dat,
    promoter = "promoter",
    compound = "compound",
    control = "DMSO",
    response = "value",
    replicate = "replicate"
  )
  rank0 <- results(fit_destress(assay, technical = "replicate", empirical_bayes = FALSE))
  rank1 <- results(fit_destress(
    assay,
    technical = "replicate",
    empirical_bayes = FALSE,
    background_rank = 1
  ))

  factor0 <- rank0[rank0$compound %in% c("C_factor1", "C_factor2"), ]
  factor1 <- rank1[rank1$compound %in% c("C_factor1", "C_factor2"), ]
  expect_gt(stats::sd(factor0$specific_effect), 0.5)
  expect_lt(max(abs(factor1$specific_effect)), 1e-8)
  expect_gt(max(abs(factor1$low_rank_effect)), 0.5)

  sparse1 <- rank1[rank1$compound == "C_specific", ]
  expected_sparse <- sparse_specific[sparse1$promoter]
  expect_equal(sparse1$specific_effect, unname(expected_sparse), tolerance = 1e-8)
})

test_that("background_rank_diagnostics detects broad low-rank structure", {
  promoters <- paste0("P", seq_len(8))
  compounds <- paste0("C", seq_len(10))
  loading <- stats::setNames(seq(-1, 1, length.out = length(promoters)), promoters)
  score <- stats::setNames(c(seq(-2, 2, length.out = 6), rep(0, 4)), compounds)
  tab <- expand.grid(
    promoter = promoters,
    compound = compounds,
    stringsAsFactors = FALSE
  )
  tab$specific_effect <- loading[tab$promoter] * score[tab$compound]

  diag <- background_rank_diagnostics(
    tab,
    rank_max = 3,
    permutations = 30,
    seed = 1
  )

  expect_equal(nrow(diag), 3)
  expect_gt(diag$observed[1], diag$null_q95[1])
  expect_gt(diag$prop_variance[1], 0.9)
})

test_that("fit_destress subtracts model-based empty-vector background before centering", {
  dat <- expand.grid(
    promoter = c("EVC", "P1", "P2"),
    compound = c("DMSO", "C_background", "C_specific"),
    replicate = paste0("r", seq_len(5)),
    stringsAsFactors = FALSE
  )
  baseline <- c(EVC = 8, P1 = 10, P2 = 12)
  background <- c(DMSO = 0, C_background = 1, C_specific = 0)
  dat$value <- baseline[dat$promoter] + background[dat$compound] +
    ifelse(dat$promoter == "P2" & dat$compound == "C_specific", 2, 0)

  assay <- prepare_assay(
    dat,
    promoter = "promoter",
    compound = "compound",
    control = "DMSO",
    response = "value",
    replicate = "replicate"
  )
  fit <- fit_destress(
    assay,
    technical = "replicate",
    empirical_bayes = FALSE,
    empty_vector_promoter = "EVC"
  )
  res <- results(fit)
  fit_no_evc <- fit_destress(
    assay,
    technical = "replicate",
    empirical_bayes = FALSE
  )
  res_no_evc <- results(fit_no_evc, promoters = c("P1", "P2"))

  expect_false("EVC" %in% res$promoter)
  bg <- res[res$compound == "C_background", ]
  expect_equal(bg$empty_vector_effect, rep(1, nrow(bg)), tolerance = 1e-8)
  expect_equal(bg$background_adjusted_effect, rep(0, nrow(bg)), tolerance = 1e-8)
  expect_equal(bg$specific_effect, rep(0, nrow(bg)), tolerance = 1e-8)

  sig <- res[res$compound == "C_specific", ]
  expect_equal(sig$background_adjusted_effect[match(c("P1", "P2"), sig$promoter)], c(0, 2), tolerance = 1e-8)
  expect_equal(sig$specific_effect[match(c("P1", "P2"), sig$promoter)], c(-1, 1), tolerance = 1e-8)

  merged <- merge(
    res[, c("promoter", "compound", "specific_effect", "specific_se")],
    res_no_evc[, c("promoter", "compound", "specific_effect", "specific_se")],
    by = c("promoter", "compound"),
    suffixes = c("_evc", "_no_evc")
  )
  expect_equal(merged$specific_effect_evc, merged$specific_effect_no_evc, tolerance = 1e-8)
  expect_equal(merged$specific_se_evc, merged$specific_se_no_evc, tolerance = 1e-8)
})

test_that("prepare_assay can calibrate responses against a background promoter", {
  compounds <- c("DMSO", paste0("C_bg", seq_len(6)), "C_specific")
  dat <- expand.grid(
    promoter = c("EVC", "P1", "P2"),
    compound = compounds,
    replicate = paste0("r", seq_len(5)),
    stringsAsFactors = FALSE
  )
  background_score <- stats::setNames(c(0, seq(-1.5, 1.5, length.out = 6), 0.5), compounds)
  dat$value <- ifelse(
    dat$promoter == "EVC",
    5 + background_score[dat$compound],
    ifelse(dat$promoter == "P1", 10 + 2 * background_score[dat$compound], 12 - background_score[dat$compound])
  )
  dat$value <- dat$value + ifelse(dat$promoter == "P2" & dat$compound == "C_specific", 2, 0)

  plain <- prepare_assay(
    dat,
    promoter = "promoter",
    compound = "compound",
    control = "DMSO",
    response = "value",
    replicate = "replicate"
  )
  calibrated <- prepare_assay(
    dat,
    promoter = "promoter",
    compound = "compound",
    control = "DMSO",
    response = "value",
    replicate = "replicate",
    background_promoter = "EVC",
    background_method = "lm",
    background_by = c("compound", "replicate")
  )

  fit_plain <- fit_destress(plain, technical = "replicate", empirical_bayes = FALSE)
  fit_calibrated <- fit_destress(calibrated, technical = "replicate", empirical_bayes = FALSE)
  res_plain <- results(fit_plain, promoters = c("P1", "P2"))
  res_calibrated <- results(fit_calibrated)
  params <- model_parameters(fit_calibrated)

  expect_false("EVC" %in% res_calibrated$promoter)
  expect_true(all(c(".background_response", ".response_uncalibrated") %in% names(calibrated)))
  expect_true(all(c("P1", "P2") %in% params$background_calibration$promoter))
  expect_gt(max(abs(res_plain$specific_effect[res_plain$compound %in% paste0("C_bg", seq_len(6))])), 0.5)
  expect_lt(max(abs(res_calibrated$specific_effect[res_calibrated$compound %in% paste0("C_bg", seq_len(6))])), 0.25)

  specific <- res_calibrated[res_calibrated$compound == "C_specific", ]
  expect_gt(specific$specific_effect[match("P2", specific$promoter)], 0.8)
  expect_lt(specific$specific_effect[match("P1", specific$promoter)], -0.8)
})

test_that("prepare_assay defaults to Huber calibration when a background promoter is supplied", {
  skip_if_not_installed("MASS")
  dat <- expand.grid(
    promoter = c("EVC", "P1"),
    compound = c("DMSO", "C1", "C2"),
    replicate = paste0("r", seq_len(3)),
    stringsAsFactors = FALSE
  )
  score <- stats::setNames(c(0, 1, 2), unique(dat$compound))
  dat$value <- ifelse(dat$promoter == "EVC", 5 + score[dat$compound], 10 + 2 * score[dat$compound])

  assay <- prepare_assay(
    dat,
    promoter = "promoter",
    compound = "compound",
    control = "DMSO",
    response = "value",
    replicate = "replicate",
    background_promoter = "EVC",
    background_by = c("compound", "replicate")
  )

  expect_equal(attr(assay, "destress")$background_method, "huber")
  expect_equal(attr(assay, "destress")$background_fit$method, "huber")
})

test_that("prepare_assay supports Huber background calibration when MASS is available", {
  skip_if_not_installed("MASS")
  dat <- expand.grid(
    promoter = c("EVC", "P1"),
    compound = c("DMSO", paste0("C", seq_len(5))),
    replicate = paste0("r", seq_len(4)),
    stringsAsFactors = FALSE
  )
  score <- stats::setNames(c(0, -2, -1, 0.5, 1, 2), unique(dat$compound))
  dat$value <- ifelse(dat$promoter == "EVC", 5 + score[dat$compound], 10 + 2 * score[dat$compound])
  dat$value[dat$promoter == "P1" & dat$compound == "C5" & dat$replicate == "r4"] <- 30

  assay <- prepare_assay(
    dat,
    promoter = "promoter",
    compound = "compound",
    control = "DMSO",
    response = "value",
    replicate = "replicate",
    background_promoter = "EVC",
    background_method = "huber",
    background_by = c("compound", "replicate")
  )
  params <- attr(assay, "destress")$background_fit

  expect_equal(attr(assay, "destress")$background_method, "huber")
  expect_equal(params$method, "huber")
  expect_true(is.finite(params$slope))
})

test_that("fit_destress can prepare raw model data with growth-exponent options", {
  dat <- simulate_screen(seed = 6, n_promoters = 5, n_compounds = 6, n_replicates = 2)

  direct <- fit_destress(
    dat,
    preset = "model",
    promoter = "promoter",
    compound = "compound",
    control = "DMSO",
    lux = "LUX.AUC_16",
    growth = "od_16h.measured",
    growth_exponent = 1,
    batch = "batch",
    replicate = "replicate",
    technical = c("batch", "replicate")
  )
  assay <- prepare_assay(
    dat,
    promoter = "promoter",
    compound = "compound",
    control = "DMSO",
    lux = "LUX.AUC_16",
    growth = "od_16h.measured",
    growth_exponent = 1,
    batch = "batch",
    replicate = "replicate"
  )
  prepared <- fit_destress(assay, technical = c("batch", "replicate"))

  expect_s3_class(direct, "destress_fit")
  expect_equal(results(direct), results(prepared), tolerance = 1e-10)
  expect_equal(unique(direct$assay_info$growth_exponent), 1)

  estimated <- fit_destress(
    dat,
    preset = "model",
    promoter = "promoter",
    compound = "compound",
    control = "DMSO",
    lux = "LUX.AUC_16",
    growth = "od_16h.measured",
    growth_exponent = "estimate",
    batch = "batch",
    replicate = "replicate",
    technical = c("batch", "replicate")
  )
  params <- model_parameters(estimated)
  expect_true(all(c("growth_exponents", "promoter_effects") %in% names(params)))
  expect_true(all(c("a_raw", "alpha_raw", "alpha_shrunk") %in% names(params$growth_exponents)))
})

test_that("fit_destress rejects unimplemented stage combinations", {
  dat <- simulate_screen(seed = 5, n_promoters = 4, n_compounds = 5, n_replicates = 2)
  assay <- prepare_assay(dat, promoter = "promoter", compound = "compound",
                         lux = "LUX.AUC_16", growth = "od_16h.measured")

  expect_error(
    fit_destress(assay, normalization = "model", testing = "gaussian_z"),
    "currently supports"
  )
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

test_that("fit_median_polish reproduces legacy median-polish residuals", {
  dat <- expand.grid(
    promoter = c("P1", "P2"),
    libplate = "lp1",
    replicate = c("r1", "r2"),
    srn_code = c("DMSO1", "DMSO2", "C1", "C2"),
    stringsAsFactors = FALSE
  )
  dat$log2.auc.16hmeasured.normed <- c(
    10.0, 10.2, 11.4, 9.4,
    10.1, 10.3, 11.5, 9.5,
    12.0, 12.2, 13.6, 11.6,
    12.1, 12.3, 13.7, 11.7
  )

  out <- fit_median_polish(
    dat,
    control = c("DMSO1", "DMSO2"),
    maxiter = 1000,
    eps = 1e-8
  )

  group <- paste(dat$promoter, dat$libplate, dat$replicate, sep = "_")
  dmso_lookup <- tapply(
    dat$log2.auc.16hmeasured.normed[dat$srn_code %in% c("DMSO1", "DMSO2")],
    group[dat$srn_code %in% c("DMSO1", "DMSO2")],
    mean
  )
  expected_log2fc <- dat$log2.auc.16hmeasured.normed - dmso_lookup[group]
  expected_mat <- matrix(
    expected_log2fc,
    nrow = length(unique(group)),
    ncol = length(unique(dat$srn_code)),
    dimnames = list(sort(unique(group)), sort(unique(dat$srn_code)))
  )
  expected_mat[cbind(match(group, rownames(expected_mat)), match(dat$srn_code, colnames(expected_mat)))] <- expected_log2fc
  expected <- stats::medpolish(expected_mat, na.rm = TRUE, maxiter = 1000, eps = 1e-8, trace.iter = FALSE)

  expect_equal(out$polished_matrix, expected$residuals, tolerance = 1e-8)
  expect_true(all(c("log2FC.polished", "zscore", "pvalue") %in% names(out$replicate_results)))
  expect_true(all(c("pvalue.adj", "hit") %in% names(out$pair_results)))
  expect_false(any(out$pair_results$srn_code %in% c("DMSO1", "DMSO2")))
})

test_that("fit_median_polish can return DMSO normality tests", {
  dat <- expand.grid(
    promoter = "P1",
    libplate = "lp1",
    replicate = "r1",
    srn_code = c(paste0("DMSO", seq_len(5)), "C1"),
    stringsAsFactors = FALSE
  )
  dat$log2.auc.16hmeasured.normed <- c(10, 10.1, 9.9, 10.2, 9.8, 11.5)

  out <- fit_median_polish(
    dat,
    control = paste0("DMSO", seq_len(5)),
    normality = TRUE,
    normality_methods = "shapiro"
  )

  expect_true(all(c(
    "promoter_libplate_replicate",
    "promoter",
    "libplate",
    "replicate",
    "n",
    "shapiro.pval",
    "lillie.pval",
    "shapiro.pval.adj"
  ) %in% names(out$normality_results)))
  expect_equal(out$normality_results$n, 5)
  expect_true(is.finite(out$normality_results$shapiro.pval))
})

test_that("fit_workflow dispatches to the median-polish workflow", {
  dat <- expand.grid(
    promoter = c("P1", "P2"),
    libplate = "lp1",
    replicate = c("r1", "r2"),
    srn_code = c("DMSO1", "DMSO2", "C1", "C2"),
    stringsAsFactors = FALSE
  )
  dat$log2.auc.16hmeasured.normed <- c(
    10.0, 10.2, 11.4, 9.4,
    10.1, 10.3, 11.5, 9.5,
    12.0, 12.2, 13.6, 11.6,
    12.1, 12.3, 13.7, 11.7
  )

  direct <- fit_median_polish(dat, control = c("DMSO1", "DMSO2"))
  via_workflow <- fit_workflow(dat, workflow = "median-polish", control = c("DMSO1", "DMSO2"))

  expect_s3_class(via_workflow, "destress_median_polish")
  expect_equal(attr(via_workflow, "destress_workflow"), "median_polish")
  expect_equal(via_workflow$pair_results, direct$pair_results, tolerance = 1e-10)
})

test_that("fit_destress can run the median-polish preset", {
  dat <- expand.grid(
    promoter = c("P1", "P2"),
    libplate = "lp1",
    replicate = c("r1", "r2"),
    srn_code = c("DMSO1", "DMSO2", "C1", "C2"),
    stringsAsFactors = FALSE
  )
  dat$log2.auc.16hmeasured.normed <- c(
    10.0, 10.2, 11.4, 9.4,
    10.1, 10.3, 11.5, 9.5,
    12.0, 12.2, 13.6, 11.6,
    12.1, 12.3, 13.7, 11.7
  )

  out <- fit_destress(dat, preset = "median_polish", control = c("DMSO1", "DMSO2"))

  expect_s3_class(out, "destress_median_polish")
  expect_equal(attr(out, "destress_preset"), "median_polish_legacy")
  expect_equal(attr(out, "destress_stages")$normalization, "median_polish")
  expect_equal(attr(out, "destress_stages")$aggregation, "max_p")
})

test_that("fit_median_polish keeps the conservative replicate p-value", {
  dat <- expand.grid(
    promoter = "P1",
    libplate = "lp1",
    replicate = c("r1", "r2"),
    srn_code = c("DMSO1", "DMSO2", "C1"),
    stringsAsFactors = FALSE
  )
  dat$log2.auc.16hmeasured.normed <- c(10, 10.2, 13, 10.1, 10.3, 13.1)

  out <- fit_median_polish(dat, control = c("DMSO1", "DMSO2"))
  c1_replicates <- out$replicate_results[out$replicate_results$srn_code == "C1", ]
  c1_pair <- out$pair_results[out$pair_results$srn_code == "C1", ]

  expect_equal(c1_pair$pvalue, max(c1_replicates$pvalue), tolerance = 1e-12)
})

test_that("fit_empty_vector_control subtracts compound-specific EVC averages", {
  dat <- expand.grid(
    promoter = c("PEVC3", "P1", "P2"),
    replicate = c("r1", "r2"),
    srn_code = c("DMSO1", "DMSO2", "C1"),
    stringsAsFactors = FALSE
  )
  dat$value <- NA_real_
  dat$value[dat$promoter == "PEVC3" & dat$srn_code == "DMSO1"] <- c(1.0, 1.2)
  dat$value[dat$promoter == "PEVC3" & dat$srn_code == "DMSO2"] <- c(1.1, 1.3)
  dat$value[dat$promoter == "PEVC3" & dat$srn_code == "C1"] <- c(2.0, 2.2)
  dat$value[dat$promoter == "P1" & dat$srn_code == "DMSO1"] <- c(1.5, 1.7)
  dat$value[dat$promoter == "P1" & dat$srn_code == "DMSO2"] <- c(1.6, 1.8)
  dat$value[dat$promoter == "P1" & dat$srn_code == "C1"] <- c(4.5, 4.7)
  dat$value[dat$promoter == "P2" & dat$srn_code == "DMSO1"] <- c(0.8, 1.0)
  dat$value[dat$promoter == "P2" & dat$srn_code == "DMSO2"] <- c(0.9, 1.1)
  dat$value[dat$promoter == "P2" & dat$srn_code == "C1"] <- c(1.5, 1.7)

  out <- fit_empty_vector_control(
    dat,
    response = "value",
    control = c("DMSO1", "DMSO2")
  )

  p1_c1 <- out$replicate_results[
    out$replicate_results$promoter == "P1" &
      out$replicate_results$srn_code == "C1",
    ,
    drop = FALSE
  ]

  expect_equal(p1_c1$empty_vector_mean, c(2.1, 2.1), tolerance = 1e-12)
  expect_equal(p1_c1$log.evcfc, c(2.4, 2.6), tolerance = 1e-12)
  expect_false(any(out$replicate_results$promoter == "PEVC3"))
  expect_true(all(c("pvalue.adj", "hit") %in% names(out$pair_results)))
})

test_that("fit_workflow dispatches to the empty-vector workflow", {
  dat <- expand.grid(
    promoter = c("PEVC3", "P1", "P2"),
    replicate = c("r1", "r2"),
    srn_code = c("DMSO1", "DMSO2", "C1"),
    stringsAsFactors = FALSE
  )
  dat$value <- NA_real_
  dat$value[dat$promoter == "PEVC3" & dat$srn_code == "DMSO1"] <- c(1.0, 1.2)
  dat$value[dat$promoter == "PEVC3" & dat$srn_code == "DMSO2"] <- c(1.1, 1.3)
  dat$value[dat$promoter == "PEVC3" & dat$srn_code == "C1"] <- c(2.0, 2.2)
  dat$value[dat$promoter == "P1" & dat$srn_code == "DMSO1"] <- c(1.5, 1.7)
  dat$value[dat$promoter == "P1" & dat$srn_code == "DMSO2"] <- c(1.6, 1.8)
  dat$value[dat$promoter == "P1" & dat$srn_code == "C1"] <- c(4.5, 4.7)
  dat$value[dat$promoter == "P2" & dat$srn_code == "DMSO1"] <- c(0.8, 1.0)
  dat$value[dat$promoter == "P2" & dat$srn_code == "DMSO2"] <- c(0.9, 1.1)
  dat$value[dat$promoter == "P2" & dat$srn_code == "C1"] <- c(1.5, 1.7)

  direct <- fit_empty_vector_control(dat, response = "value", control = c("DMSO1", "DMSO2"))
  via_workflow <- fit_workflow(dat, workflow = "evc", response = "value", control = c("DMSO1", "DMSO2"))

  expect_s3_class(via_workflow, "destress_empty_vector")
  expect_equal(attr(via_workflow, "destress_workflow"), "empty_vector_control")
  expect_equal(via_workflow$pair_results, direct$pair_results, tolerance = 1e-10)
})

test_that("fit_destress can run the empty-vector preset", {
  dat <- expand.grid(
    promoter = c("PEVC3", "P1", "P2"),
    replicate = c("r1", "r2"),
    srn_code = c("DMSO1", "DMSO2", "C1"),
    stringsAsFactors = FALSE
  )
  dat$value <- NA_real_
  dat$value[dat$promoter == "PEVC3" & dat$srn_code == "DMSO1"] <- c(1.0, 1.2)
  dat$value[dat$promoter == "PEVC3" & dat$srn_code == "DMSO2"] <- c(1.1, 1.3)
  dat$value[dat$promoter == "PEVC3" & dat$srn_code == "C1"] <- c(2.0, 2.2)
  dat$value[dat$promoter == "P1" & dat$srn_code == "DMSO1"] <- c(1.5, 1.7)
  dat$value[dat$promoter == "P1" & dat$srn_code == "DMSO2"] <- c(1.6, 1.8)
  dat$value[dat$promoter == "P1" & dat$srn_code == "C1"] <- c(4.5, 4.7)
  dat$value[dat$promoter == "P2" & dat$srn_code == "DMSO1"] <- c(0.8, 1.0)
  dat$value[dat$promoter == "P2" & dat$srn_code == "DMSO2"] <- c(0.9, 1.1)
  dat$value[dat$promoter == "P2" & dat$srn_code == "C1"] <- c(1.5, 1.7)

  out <- fit_destress(dat, preset = "evc", response = "value", control = c("DMSO1", "DMSO2"))

  expect_s3_class(out, "destress_empty_vector")
  expect_equal(attr(out, "destress_preset"), "empty_vector_control")
  expect_equal(attr(out, "destress_stages")$normalization, "empty_vector")
  expect_equal(attr(out, "destress_stages")$testing, "gaussian_z")
})

test_that("fit_empty_vector_control keeps the conservative replicate p-value", {
  dat <- expand.grid(
    promoter = c("PEVC3", "P1"),
    replicate = c("r1", "r2"),
    srn_code = c("DMSO1", "DMSO2", "C1"),
    stringsAsFactors = FALSE
  )
  dat$value <- c(
    1.0, 1.2, 1.1, 1.3, 2.0, 2.2,
    1.4, 1.6, 1.5, 1.7, 3.0, 3.4
  )

  out <- fit_empty_vector_control(dat, response = "value", control = c("DMSO1", "DMSO2"))
  c1_replicates <- out$replicate_results[out$replicate_results$srn_code == "C1", ]
  c1_pair <- out$pair_results[out$pair_results$srn_code == "C1", ]

  expect_equal(c1_pair$pvalue, max(c1_replicates$pvalue), tolerance = 1e-12)
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
