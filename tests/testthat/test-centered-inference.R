test_that("total effects are the default target and reporter-relative effects remain selectable", {
  tab <- data.frame(
    reporter = c("P1", "P1"), perturbation = c("C1", "C2"),
    total_effect = c(2, -1), total_pvalue = c(0.001, 0.8),
    total_padj = c(0.002, 0.8),
    relative_effect = c(-2, 1), relative_pvalue = c(0.8, 0.001),
    relative_padj = c(0.8, 0.002)
  )
  expect_equal(call_hits(tab)$hit, c("Upregulated", "Not DE"))
  expect_equal(call_hits(tab, effect = "relative_effect")$hit,
    c("Not DE", "Upregulated"))
  expect_error(call_hits(tab[c("relative_effect", "relative_padj")]),
    "Missing required columns: total_effect, total_padj")
  expect_equal(adjust_pvalues(tab)$total_padj_by_reporter,
    p.adjust(tab$total_pvalue, "BH"))
  expect_equal(adjust_pvalues(tab, pvalue = "relative_pvalue")$relative_padj_by_reporter,
    p.adjust(tab$relative_pvalue, "BH"))
  expect_identical(formals(plot_response_heatmap)$value, "total_effect")
  expect_identical(formals(plot_effect_histogram)$value, "total_effect")
  expect_identical(formals(plot_volcano)$effect, "total_effect")
  expect_identical(formals(summarize_hits)$effect, "total_effect")
  expect_identical(formals(plot_hit_heatmap)$effect, "total_effect")
  expect_identical(formals(plot_response_cluster_blocks)$value, "total_effect")
  expect_identical(formals(plot_response_clustered_heatmap)$value, "total_effect")
  expect_identical(formals(fit_effect_mixture)$value, "total_effect")
})

test_that("centered contrast degrees of freedom recover Welch's two-sample formula", {
  v <- c(2, 5)
  df <- c(8, 19)
  expected <- sum(v)^2 / sum(v^2 / df)
  expect_equal(centered_contrast_df(v, df), rep(expected, 2))
  expect_equal(centered_contrast_df(rep(1, 7), 10), rep(10 * 294 / 217, 7))
  expect_equal(centered_contrast_df(v * 1e200, df), rep(expected, 2))
  expect_true(is.na(centered_contrast_df(1, 8)))
  expect_true(all(is.na(centered_contrast_df(c(0, 0), 8))))
})

test_that("equal-df empirical Bayes moments retain the established calculation", {
  v <- c(0.02, 0.04, 0.09, 0.2, 0.5, 0.8)
  df <- 12
  nu0 <- 2 * trigamma_inverse(var(log(v)) - psigamma(df / 2, 1))
  s0 <- exp(mean(log(v)) - digamma(df / 2) + log(df / 2) +
              digamma(nu0 / 2) - log(nu0 / 2))
  prior <- estimate_eb_prior_variance(v, df)
  expect_equal(prior$df, nu0)
  expect_equal(prior$var, s0)
  expect_equal(prior, estimate_eb_prior_variance(v, rep(df, length(v))))
})

test_that("unequal-df empirical Bayes moments correct each sampling distribution", {
  v <- c(0.02, 0.04, 0.09, 0.2, 0.5, 0.8)
  df <- c(4, 8, 10, 20, 40, 80)
  z <- log(v) - digamma(df / 2) + log(df / 2)
  nu0 <- 2 * trigamma_inverse(var(z) - mean(psigamma(df / 2, 1)))
  s0 <- exp(mean(z) + digamma(nu0 / 2) - log(nu0 / 2))
  prior <- estimate_eb_prior_variance(v, df)
  expect_equal(prior$df, nu0)
  expect_equal(prior$var, s0)
  expect_error(estimate_eb_prior_variance(v, c(4, 8)), "same length")
})

test_that("reporter-relative uncertainty is constructed from raw variances and moderated once", {
  dat <- simulate_screen(seed = 73, n_reporters = 5, n_perturbations = 8, n_replicates = 4)
  assay <- prepare_assay(dat, reporter = "reporter", perturbation = "perturbation",
                         lux = "LUX.AUC_16", growth = "od_16h.measured",
                         batch = "batch", replicate = "replicate")
  fit <- fit_destress(assay, technical = c("batch", "replicate"), empirical_bayes = FALSE)
  plain <- results(fit)
  fit$empirical_bayes <- TRUE
  moderated <- results(fit)
  key <- paste(plain$reporter, plain$perturbation)
  moderated <- moderated[match(key, paste(moderated$reporter, moderated$perturbation)), ]
  expect_equal(moderated$relative_effect, plain$relative_effect)
  expect_equal(moderated$relative_var, plain$relative_se^2)
  expect_equal(moderated$total_var_raw, plain$total_se^2)
  expected <- eb_moderate_se(plain$relative_se, plain$relative_df_raw)
  expect_equal(moderated$relative_se, expected$se)
  expect_equal(moderated$relative_df_for_test, expected$df)
  expected_p <- 2 * pt(abs(plain$relative_effect / expected$se), expected$df, lower.tail = FALSE)
  expect_equal(moderated$relative_pvalue, expected_p)
})
test_that("scientific preset names preserve the default model", {
  expect_identical(destress_presets()[1], "model")
  expect_identical(normalize_destress_preset("model"), "model")
  expect_identical(normalize_destress_preset("median_polish"), "median_polish")
  expect_identical(normalize_destress_preset("median_polish_legacy"), "median_polish")
  expect_identical(destress_preset_stages("model")$normalization, "linear_model")
})

test_that("interaction reporter-relative tests use the centered joint contrast covariance", {
  set.seed(101)
  dat <- expand.grid(
    reporter = paste0("R", 1:4), perturbation = c("Control", "P1", "P2"),
    replicate = 1:5, stringsAsFactors = FALSE
  )
  dat <- dat[-c(3, 17, 29), ]
  dat$dose <- runif(nrow(dat))
  dat$response <- rnorm(nrow(dat)) + dat$dose
  assay <- prepare_assay(dat, reporter = "reporter", perturbation = "perturbation",
    control = "Control", response = "response", replicate = "replicate",
    batch = "dose", numeric_covariates = "dose")
  fit <- fit_destress(assay, technical = c("replicate", "dose"),
    interaction = TRUE, empirical_bayes = FALSE)
  raw <- results(fit)
  reference <- data.frame(
    .reporter = factor(raw$reporter, levels = fit$levels$reporter),
    .perturbation = factor("Control", levels = fit$levels$perturbation)
  )
  mf <- fit_model_frame(fit$full_fit)
  for (column in fit$technical) {
    x <- mf[[column]]
    reference[[column]] <- if (is.factor(x)) {
      factor(names(sort(table(x), decreasing = TRUE))[1], levels = levels(x))
    } else stats::median(x)
  }
  treatment <- reference
  treatment$.perturbation <- factor(raw$perturbation, levels = fit$levels$perturbation)
  terms <- stats::delete.response(fit_terms(fit$full_fit))
  contrast <- stats::model.matrix(terms, treatment) - stats::model.matrix(terms, reference)
  for (idx in split(seq_len(nrow(raw)), raw$perturbation)) {
    contrast[idx, ] <- sweep(contrast[idx, , drop = FALSE], 2,
      colMeans(contrast[idx, , drop = FALSE]), "-")
  }
  vc <- safe_vcov(fit$full_fit)
  contrast <- contrast[, colnames(vc), drop = FALSE]
  expected_var <- unname(rowSums((contrast %*% vc) * contrast))
  expect_equal(raw$relative_var, expected_var)
  expect_equal(raw$relative_se, sqrt(expected_var))
  expect_equal(raw$relative_effect, as.numeric(contrast %*% fit_coef(fit$full_fit)))
  expect_false(isTRUE(all.equal(raw$relative_var, raw$total_var_raw)))
  expected_test <- wald_t_test(raw$relative_effect, sqrt(expected_var), raw$residual_df)
  expect_equal(raw$relative_pvalue, expected_test$pvalue)
  fit$empirical_bayes <- TRUE
  moderated <- results(fit)
  expected <- eb_moderate_se(sqrt(expected_var), raw$residual_df)
  expect_equal(moderated$relative_var, expected_var)
  expect_equal(moderated$relative_se, expected$se)
  expect_equal(moderated$relative_df_for_test, expected$df)
  expect_equal(moderated$relative_pvalue,
    wald_t_test(raw$relative_effect, expected$se, expected$df)$pvalue)
  fit$background_rank <- 1L
  rank_adjusted <- results(fit)
  expect_equal(rank_adjusted$relative_var, moderated$relative_var)
  expect_equal(rank_adjusted$relative_se, moderated$relative_se)
  expect_equal(rank_adjusted$relative_pvalue,
    wald_t_test(rank_adjusted$relative_effect, expected$se, expected$df)$pvalue)
  single <- results(fit, reporters = "R1")
  expect_equal(single$relative_effect, rep(0, nrow(single)), tolerance = 1e-12)
  expect_true(all(is.na(single$relative_pvalue)))
})

test_that("interaction total-effect tests honor moderation independently of relative tests", {
  dat <- simulate_screen(seed = 91, n_reporters = 4, n_perturbations = 6, n_replicates = 4)
  assay <- prepare_assay(dat, reporter = "reporter", perturbation = "perturbation",
                         lux = "LUX.AUC_16", growth = "od_16h.measured",
                         batch = "batch", replicate = "replicate")
  fit <- fit_destress(assay, technical = c("batch", "replicate"),
                       interaction = TRUE, empirical_bayes = FALSE)
  raw <- results(fit)
  fit$empirical_bayes <- TRUE
  mod <- results(fit)
  expect_equal(mod$total_effect, raw$total_effect)
  expect_equal(mod$relative_effect, raw$relative_effect)
  expect_equal(mod$total_var_raw, raw$total_se^2)
  expected <- eb_moderate_se(raw$total_se, raw$residual_df)
  expect_equal(mod$total_se, expected$se)
  expect_equal(mod$total_df_for_test, expected$df)
  expected_test <- wald_t_test(raw$total_effect, expected$se, expected$df)
  expect_equal(mod$total_pvalue, expected_test$pvalue)
  expect_equal(mod$total_padj_by_reporter,
    adjust_destress_pvalues(expected_test$pvalue, mod$reporter, "by_reporter"))
  expect_equal(mod$rank_adjusted_total_pvalue, mod$total_pvalue)
  expected_relative <- eb_moderate_se(raw$relative_se, raw$residual_df)
  expect_equal(mod$relative_se, expected_relative$se)
  expect_equal(raw$total_pvalue,
    wald_t_test(raw$total_effect, raw$total_se, raw$residual_df)$pvalue)
  fit$background_rank <- 1L
  rank_mod <- results(fit)
  expect_equal(rank_mod$total_effect, mod$total_effect)
  expect_equal(rank_mod$rank_adjusted_total_se, mod$total_se)
  expect_equal(rank_mod$rank_adjusted_total_pvalue,
    wald_t_test(rank_mod$rank_adjusted_total_effect,
      mod$total_se, mod$total_df_for_test)$pvalue)
})

