test_that("prepare_assay reproduces log2 lux over growth", {
  dat <- data.frame(
    reporter = c("P1", "P1"),
    perturbation = c("DMSO", "C1"),
    lux = c(16, 32),
    growth = c(2, 2)
  )
  assay <- prepare_assay(dat, reporter = "reporter", perturbation = "perturbation",
                         lux = "lux", growth = "growth", growth_exponent = 1)
  expect_equal(assay$.response, c(3, 4), tolerance = 1e-6)
})

test_that("prepare_assay preserves requested numeric covariates", {
  dat <- data.frame(
    reporter = rep("P1", 4),
    perturbation = c("DMSO", "C1", "C1", "DMSO"),
    lux = c(16, 32, 64, 16),
    growth = c(2, 2, 2, 2),
    replicate = c(1, 1, 2, 2),
    dose_level = c(0, 1, 2, 0)
  )
  assay <- prepare_assay(
    dat,
    reporter = "reporter",
    perturbation = "perturbation",
    lux = "lux",
    growth = "growth",
    growth_exponent = 1,
    batch = "dose_level",
    replicate = "replicate",
    numeric_covariates = "dose_level"
  )

  expect_type(assay$dose_level, "double")
  expect_s3_class(assay$replicate, "factor")
})

test_that("growth-exponent estimation accepts numeric covariates", {
  dat <- expand.grid(
    reporter = c("P1", "P2"),
    dose_level = 0:3,
    replicate = 1:2,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  dat$perturbation <- "DMSO"
  dat$growth <- 2^(1 + 0.1 * dat$dose_level + 0.05 * dat$replicate +
    c(0.00, 0.04, -0.02, 0.03, -0.01, 0.05, -0.03, 0.02,
      0.01, -0.04, 0.02, -0.03, 0.04, -0.01, 0.03, -0.02))
  dat$lux <- 2^(0.5 + 0.8 * log2(dat$growth) + 0.2 * dat$dose_level)

  fit <- estimate_growth_exponents(
    dat,
    covariates = c("dose_level", "replicate"),
    numeric_covariates = "dose_level",
    min_control_n = 4,
    shrink = FALSE
  )

  expect_true(all(is.finite(fit$alpha_shrunk)))
  expect_true(all(grepl("dose_level", fit$alpha_covariates)))
})

test_that("estimate_growth_exponents recovers reporter-specific scaling", {
  dat <- expand.grid(
    reporter = c("P1", "P2"),
    growth = c(1, 2, 4, 8, 16),
    replicate = seq_len(3),
    stringsAsFactors = FALSE
  )
  dat$perturbation <- "DMSO"
  dat$lux <- ifelse(dat$reporter == "P1", 8 * dat$growth^1, 4 * dat$growth^0.5)
  est <- estimate_growth_exponents(dat, reporter = "reporter", perturbation = "perturbation",
                                   lux = "lux", growth = "growth", min_control_n = 5,
                                   shrink = FALSE)
  expect_true(all(c("a_raw", "a_raw_se", "a_raw_df") %in% names(est)))
  expect_equal(est$alpha_raw[match("P1", est$reporter)], 1, tolerance = 1e-6)
  expect_equal(est$alpha_raw[match("P2", est$reporter)], 0.5, tolerance = 1e-6)
})

test_that("growth exponent estimation adjusts technical covariates", {
  dat <- expand.grid(
    reporter = "P1",
    plate = c("A", "B"),
    replicate = seq_len(12),
    stringsAsFactors = FALSE
  )
  dat$perturbation <- "DMSO"
  dat$growth <- ifelse(dat$plate == "A", 1, 8) * rep(c(1, 1.2, 1.4, 1.6), length.out = nrow(dat))
  plate_effect <- ifelse(dat$plate == "A", 0.25, 8)
  dat$lux <- 4 * dat$growth^0.5 * plate_effect

  unadjusted <- estimate_growth_exponents(
    dat,
    reporter = "reporter",
    perturbation = "perturbation",
    lux = "lux",
    growth = "growth",
    min_control_n = 8,
    shrink = FALSE
  )
  adjusted <- estimate_growth_exponents(
    dat,
    reporter = "reporter",
    perturbation = "perturbation",
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
  dat <- simulate_screen(seed = 1, n_reporters = 8, n_perturbations = 12, n_replicates = 3)
  assay <- prepare_assay(dat, reporter = "reporter", perturbation = "perturbation",
                         lux = "LUX.AUC_16", growth = "od_16h.measured",
                         batch = "batch", replicate = "replicate")
  fit <- fit_destress(assay, technical = c("batch", "replicate"))
  res <- results(fit)
  expect_true(all(c("specific_effect", "global_effect", "total_effect") %in% names(res)))
  truth <- unique(dat[dat$perturbation != "DMSO", c("reporter", "perturbation", "truth_specific")])
  joined <- merge(res, truth, by = c("reporter", "perturbation"))
  active <- abs(joined$truth_specific) > 0
  expect_gt(stats::cor(joined$specific_effect[active], joined$truth_specific[active]), 0.7)
})

test_that("perturbation diagnostics rank mean effects and variance residuals", {
  tab <- data.frame(
    reporter = rep(paste0("P", 1:4), times = 3),
    perturbation = rep(c("quiet", "global", "mixed"), each = 4),
    total_effect = c(
      0.01, -0.01, 0.02, 0.00,
      1.00, 0.90, 1.10, 1.00,
      -1.00, 1.00, -0.80, 0.90
    )
  )

  diag <- perturbation_diagnostics(tab)

  expect_true(all(c(
    "mean_effect", "abs_mean_effect", "effect_variance",
    "rank_abs_mean_effect", "variance_residual"
  ) %in% names(diag)))
  expect_equal(diag$rank_abs_mean_effect[diag$perturbation == "quiet"], 1)
  expect_equal(diag$rank_abs_mean_effect[diag$perturbation == "global"], 3)
  expect_gt(
    diag$effect_variance[diag$perturbation == "mixed"],
    diag$effect_variance[diag$perturbation == "global"]
  )
})

test_that("mean-variance diagnostic plot stores its diagnostic table", {
  skip_if_not_installed("ggplot2")
  tab <- data.frame(
    reporter = rep(paste0("P", 1:4), times = 5),
    perturbation = rep(paste0("C", 1:5), each = 4),
    total_effect = c(
      0, 0.1, -0.1, 0,
      0.3, 0.4, 0.2, 0.3,
      0.8, 0.9, 0.7, 0.8,
      -1, 1, -1, 1,
      1.5, 1.6, 1.4, 1.5
    )
  )

  p <- plot_mean_variance_diagnostic(tab, label_by = "none")

  expect_s3_class(p, "ggplot")
  expect_true(is.data.frame(attr(p, "diagnostics")))
  expect_equal(nrow(attr(p, "diagnostics")), 5)
})

test_that("variance distribution diagnostics fit positive variance summaries", {
  set.seed(11)
  u <- stats::rbeta(80, shape1 = 1.8, shape2 = 4.5)
  diagnostics <- data.frame(
    perturbation = paste0("C", seq_along(u)),
    effect_variance = 0.25 * u / (1 - u)
  )

  fit <- fit_variance_distribution(
    diagnostics,
    seed = 2
  )

  expect_s3_class(fit, "destress_variance_distribution_fit")
  expect_true(all(c("beta-prime", "inverse-gamma", "log-normal", "log-t") %in% fit$distribution))
  expect_true(all(is.finite(fit$AIC)))
  expect_equal(attr(fit, "variance"), "effect_variance")

  mixture_fit <- fit_variance_distribution(
    diagnostics,
    distributions = "two_gaussian",
    seed = 2
  )
  expect_true("two-Gaussian mixture" %in% mixture_fit$distribution)
})

test_that("variance distribution plot stores fitted distributions", {
  skip_if_not_installed("ggplot2")
  set.seed(12)
  diagnostics <- data.frame(
    perturbation = paste0("C", seq_len(40)),
    effect_variance = exp(stats::rnorm(40, mean = -2, sd = 0.6))
  )

  p <- plot_variance_distribution(
    diagnostics,
    distributions = "beta_prime"
  )

  expect_s3_class(p, "ggplot")
  expect_s3_class(attr(p, "fit"), "destress_variance_distribution_fit")
})

test_that("Binsfeld reporter data support DStressR model analysis", {
  data("binsfeld_reporter_auc", package = "DStressR")
  data("binsfeld_reporter_scores", package = "DStressR")

  expect_equal(nrow(binsfeld_reporter_auc), 24576)
  expect_true(all(c(
    "strain", "reporter", "perturbation", "dose_level", "od_auc", "lux_auc", "removed"
  ) %in% names(binsfeld_reporter_auc)))
  expect_equal(
    sort(unique(binsfeld_reporter_auc$perturbation[grepl("^Water_", binsfeld_reporter_auc$drug)])),
    "Water"
  )
  expect_true(all(c("Scores", "Z_scores") %in% unique(binsfeld_reporter_scores$statistic)))

  wt_auc <- binsfeld_reporter_auc[
    binsfeld_reporter_auc$strain == "WT" &
      binsfeld_reporter_auc$removed == "No" &
      binsfeld_reporter_auc$perturbation %in% c("Water", "Azithromycin", "Clarithromycin"),
  ]
  assay <- prepare_assay(
    wt_auc,
    reporter = "reporter",
    perturbation = "perturbation",
    control = "Water",
    lux = "lux_auc",
    growth = "od_auc",
    growth_exponent = "estimate",
    batch = "dose_level",
    replicate = "replicate",
    growth_covariates = "replicate",
    numeric_covariates = "dose_level"
  )
  growth_fit <- attr(assay, "destress")$growth_exponent_fit
  expect_true(all(growth_fit$alpha_covariates == "replicate"))
  fit <- fit_destress(
    assay,
    technical = c("replicate", "dose_level"),
    empirical_bayes = TRUE,
    adjustment = "by_reporter",
    interaction = FALSE,
    empty_vector_reporter = "EVC"
  )
  res <- results(fit)

  expect_s3_class(fit, "destress_fit")
  expect_true(all(c("empty_vector_effect", "specific_effect", "specific_padj") %in% names(res)))
  expect_false(is.null(fit$growth_exponents))
  expect_true(nrow(res) > 0)
})

test_that("fit_workflow dispatches to the model workflow", {
  dat <- simulate_screen(seed = 3, n_reporters = 5, n_perturbations = 6, n_replicates = 2)
  assay <- prepare_assay(dat, reporter = "reporter", perturbation = "perturbation",
                         lux = "LUX.AUC_16", growth = "od_16h.measured",
                         batch = "batch", replicate = "replicate")

  direct <- fit_destress(assay, technical = c("batch", "replicate"))
  via_workflow <- fit_workflow(assay, workflow = "model", technical = c("batch", "replicate"))

  expect_s3_class(via_workflow, "destress_fit")
  expect_equal(attr(via_workflow, "destress_workflow"), "model")
  expect_equal(results(via_workflow), results(direct), tolerance = 1e-10)
})

test_that("fit_destress exposes staged model options", {
  dat <- simulate_screen(seed = 4, n_reporters = 5, n_perturbations = 6, n_replicates = 2)
  assay <- prepare_assay(dat, reporter = "reporter", perturbation = "perturbation",
                         lux = "LUX.AUC_16", growth = "od_16h.measured",
                         batch = "batch", replicate = "replicate")

  fit <- fit_destress(
    assay,
    technical = c("batch", "replicate"),
    normalization = "model",
    testing = "student_t",
    aggregation = "none",
    adjustment = "by_reporter"
  )
  res <- results(fit)
  expected <- adjust_pvalues(res, pvalue = "specific_pvalue", output = "expected_specific_padj")

  expect_s3_class(fit, "destress_fit")
  expect_false(fit$empirical_bayes)
  expect_equal(fit$stages$normalization, "linear_model")
  expect_equal(fit$stages$testing, "student_t")
  expect_equal(fit$stages$adjustment, "by_reporter")
  expect_equal(res$specific_padj, expected$expected_specific_padj, tolerance = 1e-12)
})

test_that("empirical-Bayes moderation estimates prior degrees of freedom", {
  raw_se <- sqrt(c(rep(0.02, 20), rep(0.08, 20), rep(0.20, 20)))
  moderated <- eb_moderate_se(raw_se, df = 8)

  expect_true(is.finite(moderated$prior_var))
  expect_gt(moderated$prior_df, 0)
  expect_gt(moderated$df, 8)
  expect_false(isTRUE(all.equal(moderated$se, raw_se)))
})

test_that("fit_destress can fit scalable reporter-specific models", {
  dat <- simulate_screen(seed = 7, n_reporters = 5, n_perturbations = 6, n_replicates = 3)
  assay <- prepare_assay(dat, reporter = "reporter", perturbation = "perturbation",
                         lux = "LUX.AUC_16", growth = "od_16h.measured",
                         batch = "batch", replicate = "replicate")

  fit <- fit_destress(
    assay,
    technical = c("batch", "replicate"),
    interaction = FALSE,
    empirical_bayes = FALSE,
    adjustment = "by_reporter"
  )
  res <- results(fit)

  expect_s3_class(fit, "destress_fit")
  expect_false(fit$interaction)
  expect_null(fit$full_fit)
  expect_equal(length(unique(fit$reporter_effects$reporter)), length(unique(dat$reporter)))
  expect_equal(nrow(res), length(unique(dat$reporter)) * length(setdiff(unique(dat$perturbation), "DMSO")))
  expect_true(all(c("specific_effect", "specific_pvalue", "specific_padj") %in% names(res)))
  expect_true(all(c("total_effect", "total_pvalue", "total_padj") %in% names(res)))
  expect_true(all(is.finite(res$total_pvalue)))
  expect_gt(sum(abs(res$specific_effect) > 1e-8), 0)
  expect_true(all(is.finite(res$specific_pvalue)))
  expect_true(all(res$specific_padj >= 0 & res$specific_padj <= 1))

  truth <- unique(dat[dat$perturbation != "DMSO", c("reporter", "perturbation", "truth_specific")])
  joined <- merge(res, truth, by = c("reporter", "perturbation"))
  active <- abs(joined$truth_specific) > 0
  expect_gt(stats::cor(joined$specific_effect[active], joined$truth_specific[active]), 0.7)
})

test_that("fit_destress interaction model supports numeric technical covariates", {
  dat <- simulate_screen(seed = 11, n_reporters = 4, n_perturbations = 5, n_replicates = 3)
  dat$dose_level <- rep(c(1, 2, 3), length.out = nrow(dat))
  dat$LUX.AUC_16 <- dat$LUX.AUC_16 + 0.03 * dat$dose_level

  assay <- prepare_assay(
    dat,
    reporter = "reporter",
    perturbation = "perturbation",
    lux = "LUX.AUC_16",
    growth = "od_16h.measured",
    batch = "dose_level",
    replicate = "replicate",
    numeric_covariates = "dose_level"
  )

  fit <- fit_destress(
    assay,
    technical = c("replicate", "dose_level"),
    interaction = TRUE,
    empirical_bayes = FALSE,
    adjustment = "by_reporter"
  )
  res <- results(fit)

  expect_true(fit$interaction)
  expect_equal(nrow(res), 4 * 5)
  expect_true(all(is.finite(res$specific_effect)))
  expect_true(all(is.finite(res$specific_pvalue)))
  expect_true(all(res$specific_padj_by_reporter >= 0 & res$specific_padj_by_reporter <= 1))
})

test_that("fit_destress separates global perturbation effects from reporter-specific effects", {
  dat <- expand.grid(
    reporter = paste0("P", seq_len(6)),
    perturbation = c("DMSO", "C_global", "C_specific"),
    replicate = paste0("r", seq_len(6)),
    stringsAsFactors = FALSE
  )
  baseline <- stats::setNames(seq(9.5, 10.5, length.out = 6), paste0("P", seq_len(6)))
  dat$value <- baseline[dat$reporter] +
    ifelse(dat$perturbation == "C_global", 1.5, 0) +
    ifelse(dat$perturbation == "C_specific" & dat$reporter == "P3", 1.5, 0)

  assay <- prepare_assay(
    dat,
    reporter = "reporter",
    perturbation = "perturbation",
    control = "DMSO",
    response = "value",
    replicate = "replicate"
  )
  fit <- fit_destress(assay, technical = "replicate", empirical_bayes = FALSE)
  res <- results(fit)

  global_rows <- res[res$perturbation == "C_global", ]
  expect_equal(global_rows$total_effect, rep(1.5, nrow(global_rows)), tolerance = 1e-8)
  expect_equal(unique(global_rows$global_effect), 1.5, tolerance = 1e-8)
  expect_equal(global_rows$specific_effect, rep(0, nrow(global_rows)), tolerance = 1e-8)
  expect_true(all(global_rows$specific_pvalue > 0.9))
  expect_true(all(global_rows$global_pvalue < 1e-8))

  specific_row <- res[res$perturbation == "C_specific" & res$reporter == "P3", ]
  expect_gt(specific_row$specific_effect, 1)
  expect_lt(specific_row$specific_pvalue, 1e-8)
})

test_that("fit_destress can remove a low-rank perturbation background", {
  reporters <- paste0("P", seq_len(6))
  perturbations <- c("DMSO", "C_factor1", "C_factor2", "C_specific")
  dat <- expand.grid(
    reporter = reporters,
    perturbation = perturbations,
    replicate = paste0("r", seq_len(5)),
    stringsAsFactors = FALSE
  )

  baseline <- stats::setNames(seq(9.5, 10.5, length.out = length(reporters)), reporters)
  loading <- stats::setNames(c(-2, -1, 0, 0, 1, 2), reporters)
  score <- c(DMSO = 0, C_factor1 = 1.2, C_factor2 = -0.8, C_specific = 0)
  sparse_specific <- stats::setNames(c(1, -2, 1, 1, -2, 1), reporters)
  dat$value <- baseline[dat$reporter] +
    loading[dat$reporter] * score[dat$perturbation] +
    ifelse(dat$perturbation == "C_specific", sparse_specific[dat$reporter], 0)

  assay <- prepare_assay(
    dat,
    reporter = "reporter",
    perturbation = "perturbation",
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

  factor0 <- rank0[rank0$perturbation %in% c("C_factor1", "C_factor2"), ]
  factor1 <- rank1[rank1$perturbation %in% c("C_factor1", "C_factor2"), ]
  expect_gt(stats::sd(factor0$specific_effect), 0.5)
  expect_lt(max(abs(factor1$rank_adjusted_total_effect)), 1e-8)
  expect_lt(max(abs(factor1$specific_effect)), 1e-8)
  expect_gt(max(abs(factor1$low_rank_effect)), 0.5)

  sparse1 <- rank1[rank1$perturbation == "C_specific", ]
  expected_sparse <- sparse_specific[sparse1$reporter]
  expect_equal(sparse1$specific_effect, unname(expected_sparse), tolerance = 1e-8)
})

test_that("background_rank_diagnostics detects broad low-rank structure", {
  reporters <- paste0("P", seq_len(8))
  perturbations <- paste0("C", seq_len(10))
  loading <- stats::setNames(seq(-1, 1, length.out = length(reporters)), reporters)
  score <- stats::setNames(c(seq(-2, 2, length.out = 6), rep(0, 4)), perturbations)
  tab <- expand.grid(
    reporter = reporters,
    perturbation = perturbations,
    stringsAsFactors = FALSE
  )
  tab$total_effect <- loading[tab$reporter] * score[tab$perturbation]

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
    reporter = c("EVC", "P1", "P2"),
    perturbation = c("DMSO", "C_background", "C_specific"),
    replicate = paste0("r", seq_len(5)),
    stringsAsFactors = FALSE
  )
  baseline <- c(EVC = 8, P1 = 10, P2 = 12)
  background <- c(DMSO = 0, C_background = 1, C_specific = 0)
  dat$value <- baseline[dat$reporter] + background[dat$perturbation] +
    ifelse(dat$reporter == "P2" & dat$perturbation == "C_specific", 2, 0)

  assay <- prepare_assay(
    dat,
    reporter = "reporter",
    perturbation = "perturbation",
    control = "DMSO",
    response = "value",
    replicate = "replicate"
  )
  fit <- fit_destress(
    assay,
    technical = "replicate",
    empirical_bayes = FALSE,
    empty_vector_reporter = "EVC"
  )
  res <- results(fit)
  fit_no_evc <- fit_destress(
    assay,
    technical = "replicate",
    empirical_bayes = FALSE
  )
  res_no_evc <- results(fit_no_evc, reporters = c("P1", "P2"))

  expect_false("EVC" %in% res$reporter)
  bg <- res[res$perturbation == "C_background", ]
  expect_equal(bg$empty_vector_effect, rep(1, nrow(bg)), tolerance = 1e-8)
  expect_equal(bg$background_adjusted_effect, rep(0, nrow(bg)), tolerance = 1e-8)
  expect_equal(bg$specific_effect, rep(0, nrow(bg)), tolerance = 1e-8)

  sig <- res[res$perturbation == "C_specific", ]
  expect_equal(sig$background_adjusted_effect[match(c("P1", "P2"), sig$reporter)], c(0, 2), tolerance = 1e-8)
  expect_equal(sig$specific_effect[match(c("P1", "P2"), sig$reporter)], c(-1, 1), tolerance = 1e-8)

  merged <- merge(
    res[, c("reporter", "perturbation", "specific_effect", "specific_se")],
    res_no_evc[, c("reporter", "perturbation", "specific_effect", "specific_se")],
    by = c("reporter", "perturbation"),
    suffixes = c("_evc", "_no_evc")
  )
  expect_equal(merged$specific_effect_evc, merged$specific_effect_no_evc, tolerance = 1e-8)
  expect_equal(merged$specific_se_evc, merged$specific_se_no_evc, tolerance = 1e-8)
})

test_that("prepare_assay can calibrate responses against a background reporter", {
  perturbations <- c("DMSO", paste0("C_bg", seq_len(6)), "C_specific")
  dat <- expand.grid(
    reporter = c("EVC", "P1", "P2"),
    perturbation = perturbations,
    replicate = paste0("r", seq_len(5)),
    stringsAsFactors = FALSE
  )
  background_score <- stats::setNames(c(0, seq(-1.5, 1.5, length.out = 6), 0.5), perturbations)
  dat$value <- ifelse(
    dat$reporter == "EVC",
    5 + background_score[dat$perturbation],
    ifelse(dat$reporter == "P1", 10 + 2 * background_score[dat$perturbation], 12 - background_score[dat$perturbation])
  )
  dat$value <- dat$value + ifelse(dat$reporter == "P2" & dat$perturbation == "C_specific", 2, 0)

  plain <- prepare_assay(
    dat,
    reporter = "reporter",
    perturbation = "perturbation",
    control = "DMSO",
    response = "value",
    replicate = "replicate"
  )
  calibrated <- prepare_assay(
    dat,
    reporter = "reporter",
    perturbation = "perturbation",
    control = "DMSO",
    response = "value",
    replicate = "replicate",
    background_reporter = "EVC",
    background_method = "lm",
    background_by = c("perturbation", "replicate")
  )

  fit_plain <- fit_destress(plain, technical = "replicate", empirical_bayes = FALSE)
  fit_calibrated <- fit_destress(calibrated, technical = "replicate", empirical_bayes = FALSE)
  res_plain <- results(fit_plain, reporters = c("P1", "P2"))
  res_calibrated <- results(fit_calibrated)
  params <- model_parameters(fit_calibrated)

  expect_false("EVC" %in% res_calibrated$reporter)
  expect_true(all(c(".background_response", ".response_uncalibrated") %in% names(calibrated)))
  expect_true(all(c("P1", "P2") %in% params$background_calibration$reporter))
  expect_gt(max(abs(res_plain$specific_effect[res_plain$perturbation %in% paste0("C_bg", seq_len(6))])), 0.5)
  expect_lt(max(abs(res_calibrated$specific_effect[res_calibrated$perturbation %in% paste0("C_bg", seq_len(6))])), 0.25)

  specific <- res_calibrated[res_calibrated$perturbation == "C_specific", ]
  expect_gt(specific$specific_effect[match("P2", specific$reporter)], 0.8)
  expect_lt(specific$specific_effect[match("P1", specific$reporter)], -0.8)
})

test_that("prepare_assay defaults to Huber calibration when a background reporter is supplied", {
  skip_if_not_installed("MASS")
  dat <- expand.grid(
    reporter = c("EVC", "P1"),
    perturbation = c("DMSO", "C1", "C2"),
    replicate = paste0("r", seq_len(3)),
    stringsAsFactors = FALSE
  )
  score <- stats::setNames(c(0, 1, 2), unique(dat$perturbation))
  dat$value <- ifelse(dat$reporter == "EVC", 5 + score[dat$perturbation], 10 + 2 * score[dat$perturbation])

  assay <- prepare_assay(
    dat,
    reporter = "reporter",
    perturbation = "perturbation",
    control = "DMSO",
    response = "value",
    replicate = "replicate",
    background_reporter = "EVC",
    background_by = c("perturbation", "replicate")
  )

  expect_equal(attr(assay, "destress")$background_method, "huber")
  expect_equal(attr(assay, "destress")$background_fit$method, "huber")
})

test_that("prepare_assay supports Huber background calibration when MASS is available", {
  skip_if_not_installed("MASS")
  dat <- expand.grid(
    reporter = c("EVC", "P1"),
    perturbation = c("DMSO", paste0("C", seq_len(5))),
    replicate = paste0("r", seq_len(4)),
    stringsAsFactors = FALSE
  )
  score <- stats::setNames(c(0, -2, -1, 0.5, 1, 2), unique(dat$perturbation))
  dat$value <- ifelse(dat$reporter == "EVC", 5 + score[dat$perturbation], 10 + 2 * score[dat$perturbation])
  dat$value[dat$reporter == "P1" & dat$perturbation == "C5" & dat$replicate == "r4"] <- 30

  assay <- prepare_assay(
    dat,
    reporter = "reporter",
    perturbation = "perturbation",
    control = "DMSO",
    response = "value",
    replicate = "replicate",
    background_reporter = "EVC",
    background_method = "huber",
    background_by = c("perturbation", "replicate")
  )
  params <- attr(assay, "destress")$background_fit

  expect_equal(attr(assay, "destress")$background_method, "huber")
  expect_equal(params$method, "huber")
  expect_true(is.finite(params$slope))
})

test_that("fit_destress can prepare raw model data with growth-exponent options", {
  dat <- simulate_screen(seed = 6, n_reporters = 5, n_perturbations = 6, n_replicates = 2)

  direct <- fit_destress(
    dat,
    preset = "model",
    reporter = "reporter",
    perturbation = "perturbation",
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
    reporter = "reporter",
    perturbation = "perturbation",
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
    reporter = "reporter",
    perturbation = "perturbation",
    control = "DMSO",
    lux = "LUX.AUC_16",
    growth = "od_16h.measured",
    growth_exponent = "estimate",
    batch = "batch",
    replicate = "replicate",
    technical = c("batch", "replicate")
  )
  params <- model_parameters(estimated)
  expect_true(all(c("growth_exponents", "reporter_effects") %in% names(params)))
  expect_true(all(c("a_raw", "alpha_raw", "alpha_shrunk") %in% names(params$growth_exponents)))
})

test_that("fit_destress rejects unimplemented stage combinations", {
  dat <- simulate_screen(seed = 5, n_reporters = 4, n_perturbations = 5, n_replicates = 2)
  assay <- prepare_assay(dat, reporter = "reporter", perturbation = "perturbation",
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
    reporter = "P1",
    perturbation = paste0("C", seq_len(180)),
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
    "local_fdr_qvalue_by_reporter"
  ) %in% names(out)))
  expect_equal(nrow(summary), 1)
  expect_lt(summary$location_repressed, summary$location_null)
  expect_lt(summary$location_null, summary$location_activated)
  expect_gt(mean(out$posterior_class == out$truth), 0.85)
  expect_true(all(out$empirical_null_padj_by_reporter >= 0 & out$empirical_null_padj_by_reporter <= 1))
  expect_true(all(out$local_fdr_qvalue_by_reporter >= 0 & out$local_fdr_qvalue_by_reporter <= 1))
})

test_that("fit_median_polish reproduces legacy median-polish residuals", {
  dat <- expand.grid(
    reporter = c("P1", "P2"),
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

  group <- paste(dat$reporter, dat$libplate, dat$replicate, sep = "_")
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
    reporter = "P1",
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
    "reporter",
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
    reporter = c("P1", "P2"),
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
    reporter = c("P1", "P2"),
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
    reporter = "P1",
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

test_that("fit_empty_vector_control subtracts perturbation-specific EVC averages", {
  dat <- expand.grid(
    reporter = c("PEVC3", "P1", "P2"),
    replicate = c("r1", "r2"),
    srn_code = c("DMSO1", "DMSO2", "C1"),
    stringsAsFactors = FALSE
  )
  dat$value <- NA_real_
  dat$value[dat$reporter == "PEVC3" & dat$srn_code == "DMSO1"] <- c(1.0, 1.2)
  dat$value[dat$reporter == "PEVC3" & dat$srn_code == "DMSO2"] <- c(1.1, 1.3)
  dat$value[dat$reporter == "PEVC3" & dat$srn_code == "C1"] <- c(2.0, 2.2)
  dat$value[dat$reporter == "P1" & dat$srn_code == "DMSO1"] <- c(1.5, 1.7)
  dat$value[dat$reporter == "P1" & dat$srn_code == "DMSO2"] <- c(1.6, 1.8)
  dat$value[dat$reporter == "P1" & dat$srn_code == "C1"] <- c(4.5, 4.7)
  dat$value[dat$reporter == "P2" & dat$srn_code == "DMSO1"] <- c(0.8, 1.0)
  dat$value[dat$reporter == "P2" & dat$srn_code == "DMSO2"] <- c(0.9, 1.1)
  dat$value[dat$reporter == "P2" & dat$srn_code == "C1"] <- c(1.5, 1.7)

  out <- fit_empty_vector_control(
    dat,
    response = "value",
    control = c("DMSO1", "DMSO2")
  )

  p1_c1 <- out$replicate_results[
    out$replicate_results$reporter == "P1" &
      out$replicate_results$srn_code == "C1",
    ,
    drop = FALSE
  ]

  expect_equal(p1_c1$empty_vector_mean, c(2.1, 2.1), tolerance = 1e-12)
  expect_equal(p1_c1$log.evcfc, c(2.4, 2.6), tolerance = 1e-12)
  expect_false(any(out$replicate_results$reporter == "PEVC3"))
  expect_true(all(c("pvalue.adj", "hit") %in% names(out$pair_results)))
})

test_that("fit_workflow dispatches to the empty-vector workflow", {
  dat <- expand.grid(
    reporter = c("PEVC3", "P1", "P2"),
    replicate = c("r1", "r2"),
    srn_code = c("DMSO1", "DMSO2", "C1"),
    stringsAsFactors = FALSE
  )
  dat$value <- NA_real_
  dat$value[dat$reporter == "PEVC3" & dat$srn_code == "DMSO1"] <- c(1.0, 1.2)
  dat$value[dat$reporter == "PEVC3" & dat$srn_code == "DMSO2"] <- c(1.1, 1.3)
  dat$value[dat$reporter == "PEVC3" & dat$srn_code == "C1"] <- c(2.0, 2.2)
  dat$value[dat$reporter == "P1" & dat$srn_code == "DMSO1"] <- c(1.5, 1.7)
  dat$value[dat$reporter == "P1" & dat$srn_code == "DMSO2"] <- c(1.6, 1.8)
  dat$value[dat$reporter == "P1" & dat$srn_code == "C1"] <- c(4.5, 4.7)
  dat$value[dat$reporter == "P2" & dat$srn_code == "DMSO1"] <- c(0.8, 1.0)
  dat$value[dat$reporter == "P2" & dat$srn_code == "DMSO2"] <- c(0.9, 1.1)
  dat$value[dat$reporter == "P2" & dat$srn_code == "C1"] <- c(1.5, 1.7)

  direct <- fit_empty_vector_control(dat, response = "value", control = c("DMSO1", "DMSO2"))
  via_workflow <- fit_workflow(dat, workflow = "evc", response = "value", control = c("DMSO1", "DMSO2"))

  expect_s3_class(via_workflow, "destress_empty_vector")
  expect_equal(attr(via_workflow, "destress_workflow"), "empty_vector_control")
  expect_equal(via_workflow$pair_results, direct$pair_results, tolerance = 1e-10)
})

test_that("fit_destress can run the empty-vector preset", {
  dat <- expand.grid(
    reporter = c("PEVC3", "P1", "P2"),
    replicate = c("r1", "r2"),
    srn_code = c("DMSO1", "DMSO2", "C1"),
    stringsAsFactors = FALSE
  )
  dat$value <- NA_real_
  dat$value[dat$reporter == "PEVC3" & dat$srn_code == "DMSO1"] <- c(1.0, 1.2)
  dat$value[dat$reporter == "PEVC3" & dat$srn_code == "DMSO2"] <- c(1.1, 1.3)
  dat$value[dat$reporter == "PEVC3" & dat$srn_code == "C1"] <- c(2.0, 2.2)
  dat$value[dat$reporter == "P1" & dat$srn_code == "DMSO1"] <- c(1.5, 1.7)
  dat$value[dat$reporter == "P1" & dat$srn_code == "DMSO2"] <- c(1.6, 1.8)
  dat$value[dat$reporter == "P1" & dat$srn_code == "C1"] <- c(4.5, 4.7)
  dat$value[dat$reporter == "P2" & dat$srn_code == "DMSO1"] <- c(0.8, 1.0)
  dat$value[dat$reporter == "P2" & dat$srn_code == "DMSO2"] <- c(0.9, 1.1)
  dat$value[dat$reporter == "P2" & dat$srn_code == "C1"] <- c(1.5, 1.7)

  out <- fit_destress(dat, preset = "evc", response = "value", control = c("DMSO1", "DMSO2"))

  expect_s3_class(out, "destress_empty_vector")
  expect_equal(attr(out, "destress_preset"), "empty_vector_control")
  expect_equal(attr(out, "destress_stages")$normalization, "empty_vector")
  expect_equal(attr(out, "destress_stages")$testing, "gaussian_z")
})

test_that("fit_empty_vector_control keeps the conservative replicate p-value", {
  dat <- expand.grid(
    reporter = c("PEVC3", "P1"),
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
    reporter = "P1",
    libplate = "lp1",
    replicate = c("r1", "r2"),
    perturbation = c(paste0("DMSO", seq_len(6)), "C_high", "C_mid"),
    stringsAsFactors = FALSE
  )
  control_means <- c(-0.2, -0.1, -0.05, 0.05, 0.1, 0.2)
  dat$value <- 0
  for (i in seq_along(control_means)) {
    dat$value[dat$perturbation == paste0("DMSO", i)] <- control_means[i] +
      ifelse(dat$replicate[dat$perturbation == paste0("DMSO", i)] == "r1", -0.01, 0.01)
  }
  dat$value[dat$perturbation == "C_high"] <- c(0.95, 1.05)
  dat$value[dat$perturbation == "C_mid"] <- c(0.04, 0.06)

  out <- empirical_replicate_pvalues(
    dat,
    value = "value",
    reporter = "reporter",
    perturbation = "perturbation",
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
  expect_lt(out$empirical_pvalue[out$perturbation == "C_high"], out$empirical_pvalue[out$perturbation == "C_mid"])
  expect_lt(out$permutation_pvalue[out$perturbation == "C_high"], out$permutation_pvalue[out$perturbation == "C_mid"])
  expect_equal(out$n_replicates[out$perturbation == "C_high"], 2)
  expect_equal(out$null_n[out$perturbation == "C_high"], 6)
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
    reporter = c("P1", "P1", "P2", "P3"),
    perturbation = c("C1", "C2", "C1", "C3"),
    compound_name = c("Drug A", "Drug B", "Drug A", "Drug C"),
    specific_effect = c(2.1, -0.3, -1.8, 0.5),
    specific_padj = c(0.001, 0.7, 0.02, 0.4),
    stringsAsFactors = FALSE
  )

  p <- plot_volcano(tab, perturbation_label = "compound_name", top_n = 2, top_reporters = 2)

  expect_s3_class(p, "ggplot")
})

test_that("plot_response_heatmap returns a ggplot with matrix attribute", {
  skip_if_not_installed("ggplot2")
  tab <- expand.grid(
    reporter = c("P1", "P2"),
    perturbation = c("C1", "C2", "C3"),
    stringsAsFactors = FALSE
  )
  tab$compound_name <- c("Drug A", "Drug A", "Drug B", "Drug B", "Drug C", "Drug C")
  tab$specific_effect <- c(1, -1, 0.2, -0.2, 2, -2)

  p <- plot_response_heatmap(
    tab,
    perturbation_label = "compound_name",
    top_n_perturbations = Inf,
    cluster_rows = FALSE,
    cluster_cols = FALSE
  )
  mat <- attr(p, "response_matrix")

  expect_s3_class(p, "ggplot")
  expect_equal(dim(mat), c(2, 3))
  expect_equal(rownames(mat), c("P1", "P2"))
})

test_that("plot_response_heatmap uses global reporter order unless clustered", {
  skip_if_not_installed("ggplot2")
  tab <- expand.grid(
    reporter = c("soxSp", "acrABp", "robp"),
    perturbation = c("C1", "C2"),
    stringsAsFactors = FALSE
  )
  tab$specific_effect <- seq_len(nrow(tab))

  p <- plot_response_heatmap(
    tab,
    top_n_perturbations = Inf,
    cluster_cols = FALSE
  )

  expect_equal(rownames(attr(p, "response_matrix")), c("acrABp", "robp", "soxSp"))
})

test_that("summarize_hits returns pair, reporter, and perturbation summaries", {
  tab <- expand.grid(
    reporter = c("P1", "P2"),
    perturbation = c("C1", "C2", "C3"),
    stringsAsFactors = FALSE
  )
  tab$compound_name <- c("Drug A", "Drug A", "Drug B", "Drug B", "Drug C", "Drug C")
  tab$specific_effect <- c(1.4, -0.4, 0.2, -1.3, 1.8, 0.1)
  tab$specific_padj_by_reporter <- c(0.01, 0.4, 0.7, 0.03, 0.001, 0.9)

  hit_summary <- summarize_hits(tab, perturbation_label = "compound_name", fdr = 0.05)

  expect_s3_class(hit_summary, "destress_hit_summary")
  expect_equal(sum(hit_summary$pairs$hit), 3)
  expect_equal(hit_summary$reporters$n_hits[hit_summary$reporters$reporter == "P1"], 2)
  expect_equal(hit_summary$perturbations$n_hits[hit_summary$perturbations$perturbation_label == "Drug B"], 1)
})

test_that("plot_hit_heatmap returns a ggplot with hit summaries", {
  skip_if_not_installed("ggplot2")
  tab <- expand.grid(
    reporter = c("P1", "P2"),
    perturbation = c("C1", "C2", "C3"),
    stringsAsFactors = FALSE
  )
  tab$compound_name <- c("Drug A", "Drug A", "Drug B", "Drug B", "Drug C", "Drug C")
  tab$specific_effect <- c(1.4, -0.4, 0.2, -1.3, 1.8, 0.1)
  tab$specific_padj_by_reporter <- c(0.01, 0.4, 0.7, 0.03, 0.001, 0.9)

  p <- plot_hit_heatmap(
    tab,
    perturbation_label = "compound_name",
    top_n_perturbations = Inf,
    order_rows = "input",
    order_cols = "input",
    show_perturbation_labels = TRUE
  )

  expect_s3_class(p, "ggplot")
  expect_equal(dim(attr(p, "hit_matrix")), c(2, 3))
  expect_s3_class(attr(p, "hit_summary"), "destress_hit_summary")
  expect_equal(sum(attr(p, "hit_matrix") != 0), 3)
})

test_that("plot_hit_heatmap uses global reporter order by default", {
  skip_if_not_installed("ggplot2")
  tab <- expand.grid(
    reporter = c("soxSp", "acrABp", "robp"),
    perturbation = c("C1", "C2"),
    stringsAsFactors = FALSE
  )
  tab$specific_effect <- c(1, -1, 0.5, -0.5, 1.2, -1.2)
  tab$specific_padj_by_reporter <- c(0.01, 0.02, 0.6, 0.7, 0.03, 0.04)

  p <- plot_hit_heatmap(
    tab,
    top_n_perturbations = Inf,
    order_cols = "input"
  )

  expect_equal(rownames(attr(p, "hit_matrix")), c("acrABp", "robp", "soxSp"))
})

test_that("plot_hit_heatmap drops perturbations without hits by default", {
  skip_if_not_installed("ggplot2")
  tab <- expand.grid(
    reporter = c("P1", "P2"),
    perturbation = c("C1", "C2", "C3"),
    stringsAsFactors = FALSE
  )
  tab$specific_effect <- c(1, -1, 0.2, -0.2, 1.2, -1.2)
  tab$specific_padj_by_reporter <- c(0.01, 0.02, 0.8, 0.9, 0.03, 0.04)

  p <- plot_hit_heatmap(tab, top_n_perturbations = Inf, order_cols = "input")

  expect_equal(colnames(attr(p, "hit_matrix")), c("C1 [C1]", "C3 [C3]"))
})

test_that("heatmaps label top perturbations by absolute column sum", {
  skip_if_not_installed("ggplot2")
  tab <- expand.grid(
    reporter = c("P1", "P2"),
    perturbation = paste0("C", 1:5),
    stringsAsFactors = FALSE
  )
  tab$specific_effect <- 0.1
  tab$specific_effect[tab$perturbation == "C5"] <- 10

  p <- plot_response_heatmap(
    tab,
    top_n_perturbations = Inf,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    top_perturbation_labels = 2
  )
  built <- ggplot2::ggplot_build(p)
  axis_labels <- built$layout$panel_params[[1]]$x$get_labels()

  expect_true("C5 [C5]" %in% axis_labels)
  expect_equal(sum(nzchar(axis_labels)), 2)
})

test_that("automatic perturbation labels enforce spacing", {
  skip_if_not_installed("ggplot2")
  tab <- expand.grid(
    reporter = c("P1", "P2"),
    perturbation = paste0("C", 1:8),
    stringsAsFactors = FALSE
  )
  tab$specific_effect <- 0.1
  tab$specific_effect[tab$perturbation == "C4"] <- 5
  tab$specific_effect[tab$perturbation == "C5"] <- 6
  tab$specific_effect[tab$perturbation == "C6"] <- 4

  p <- plot_response_heatmap(
    tab,
    top_n_perturbations = Inf,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    top_perturbation_labels = 3,
    perturbation_label_min_gap = 2
  )
  built <- ggplot2::ggplot_build(p)
  axis_labels <- built$layout$panel_params[[1]]$x$get_labels()

  expect_true("C5 [C5]" %in% axis_labels)
  expect_lt(sum(nzchar(axis_labels[4:6])), 3)
})

test_that("plot_effect_histogram returns pooled and reporter-faceted plots", {
  skip_if_not_installed("ggplot2")
  tab <- data.frame(
    reporter = rep(c("P1", "P2"), each = 6),
    specific_effect = c(-1, -0.5, -0.1, 0, 0.3, 1.2, -0.2, 0, 0.1, 0.4, 0.9, 1.4)
  )

  pooled <- plot_effect_histogram(tab, bins = 10)
  per_promoter <- plot_effect_histogram(tab, by = "reporter", bins = 10)

  expect_s3_class(pooled, "ggplot")
  expect_s3_class(per_promoter, "ggplot")
})

test_that("plot_response_cluster_blocks returns cluster summaries", {
  skip_if_not_installed("ggplot2")
  tab <- expand.grid(
    reporter = c("P1", "P2", "P3", "P4"),
    perturbation = c("C1", "C2", "C3", "C4", "C5"),
    stringsAsFactors = FALSE
  )
  tab$compound_name <- paste("Drug", tab$perturbation)
  tab$specific_effect <- c(
    1.2, 1.1, -0.1, -0.2,
    1.0, 0.9, -0.2, -0.1,
    -0.2, -0.1, 1.3, 1.1,
    -0.1, -0.2, 1.1, 1.2,
    0.1, 0.2, 0.0, -0.1
  )

  p <- plot_response_cluster_blocks(
    tab,
    perturbation_label = "compound_name",
    n_reporter_clusters = 2,
    n_perturbation_clusters = 2
  )

  expect_s3_class(p, "ggplot")
  expect_equal(dim(attr(p, "response_matrix")), c(4, 5))
  expect_equal(nrow(attr(p, "reporter_clusters")), 4)
  expect_equal(nrow(attr(p, "perturbation_clusters")), 5)
  expect_equal(nrow(attr(p, "block_summary")), 4)
})

test_that("plot_response_clustered_heatmap writes clustered heatmap output", {
  tab <- expand.grid(
    reporter = c("P1", "P2", "P3", "P4"),
    perturbation = c("C1", "C2", "C3", "C4", "C5"),
    stringsAsFactors = FALSE
  )
  tab$compound_name <- paste("Drug", tab$perturbation)
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
    perturbation_label = "compound_name",
    n_reporter_clusters = 2,
    n_perturbation_clusters = 2,
    file = out_file,
    width = 8,
    height = 6,
    show_colnames = FALSE
  )

  expect_true(file.exists(out_file))
  expect_equal(dim(out$response_matrix), c(4, 5))
  expect_equal(nrow(out$reporter_clusters), 4)
  expect_equal(nrow(out$perturbation_clusters), 5)
})
