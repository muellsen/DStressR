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

test_that("specific uncertainty is constructed from raw variances and moderated once", {
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
  expect_equal(moderated$specific_effect, plain$specific_effect)
  expect_equal(moderated$specific_var, plain$specific_se^2)
  expect_equal(moderated$total_var_raw, plain$total_se^2)
  expected <- eb_moderate_se(plain$specific_se, plain$specific_df_raw)
  expect_equal(moderated$specific_se, expected$se)
  expect_equal(moderated$specific_df_for_test, expected$df)
  expected_p <- 2 * pt(abs(plain$specific_effect / expected$se), expected$df, lower.tail = FALSE)
  expect_equal(moderated$specific_pvalue, expected_p)
})
test_that("scientific preset names preserve the default model", {
  expect_identical(destress_presets()[1], "model")
  expect_identical(normalize_destress_preset("model"), "model")
  expect_identical(normalize_destress_preset("median_polish"), "median_polish")
  expect_identical(normalize_destress_preset("median_polish_legacy"), "median_polish")
  expect_identical(destress_preset_stages("model")$normalization, "linear_model")
})

