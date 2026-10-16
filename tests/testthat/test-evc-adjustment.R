evc_adjustment_assay <- function(method = "none") {
  set.seed(103)
  d <- expand.grid(reporter = c("EVC", "P1", "P2", "P3"),
                   perturbation = c("Water", paste0("C", 1:12)),
                   replicate = 1:6, stringsAsFactors = FALSE)
  d <- d[-c(1, 7, 18, 33, 49), ]
  d$value <- as.numeric(factor(d$reporter)) +
    as.numeric(factor(d$perturbation)) / 10 +
    rnorm(nrow(d), sd = as.numeric(factor(d$reporter)) / 3)
  prepare_assay(d, reporter = "reporter", perturbation = "perturbation",
    control = "Water", response = "value", replicate = "replicate",
    background_reporter = "EVC", background_method = method,
    background_by = c("perturbation", "replicate"))
}

test_that("response-level and effect-level EVC corrections fail early when combined", {
  for (method in c("subtract", "lm", "huber")) {
    assay <- evc_adjustment_assay(method)
    expect_error(fit_destress(assay, empty_vector_reporter = "EVC"),
                 "Conflicting EVC settings")
  }
  assay <- evc_adjustment_assay()
  expect_s3_class(fit_destress(assay, empty_vector_reporter = "EVC"), "destress_fit")
  expect_error(fit_destress(assay, empty_vector_reporter = "EVC", interaction = TRUE),
               "requires `interaction = FALSE`", fixed = TRUE)
})

test_that("fitted EVC subtraction combines raw variances and Welch degrees of freedom", {
  fit <- fit_destress(evc_adjustment_assay(), technical = "replicate",
                     empty_vector_reporter = "EVC", empirical_bayes = FALSE)
  res <- results(fit)
  raw <- fit$reporter_effects
  key <- function(d) paste(d$reporter, d$perturbation)
  reporter <- raw[match(key(res), key(raw)), ]
  evc <- raw[raw$reporter == "EVC", ]
  evc <- evc[match(res$perturbation, evc$perturbation), ]
  v <- reporter$total_se^2 + evc$total_se^2
  df <- v^2 / (reporter$total_se^4 / reporter$residual_df +
                 evc$total_se^4 / evc$residual_df)
  effect <- reporter$total_effect - evc$total_effect
  expect_equal(res$background_adjusted_effect, effect)
  expect_equal(res$background_adjusted_var_raw, v)
  expect_equal(res$background_adjusted_df_raw, df)
  expect_equal(res$rank_adjusted_total_se, sqrt(v))
  expect_equal(res$rank_adjusted_total_df_for_test, df)
  expect_equal(res$rank_adjusted_total_pvalue,
               2 * pt(-abs(effect / sqrt(v)), df))
  expect_false(any(res$reporter == "EVC"))
})

test_that("EVC-adjusted contrast variances are moderated once after combination", {
  fit <- fit_destress(evc_adjustment_assay(), technical = "replicate",
                     empty_vector_reporter = "EVC", empirical_bayes = FALSE)
  plain <- results(fit)
  fit$empirical_bayes <- TRUE
  moderated <- results(fit)
  expected <- eb_moderate_se(sqrt(plain$background_adjusted_var_raw),
                            plain$background_adjusted_df_raw)
  expect_equal(moderated$background_adjusted_effect, plain$background_adjusted_effect)
  expect_equal(moderated$background_adjusted_var_raw, plain$background_adjusted_var_raw)
  expect_equal(moderated$rank_adjusted_total_se, expected$se)
  expect_equal(moderated$rank_adjusted_total_df_for_test, expected$df)
  expect_equal(moderated$rank_adjusted_total_pvalue,
    2 * pt(-abs(plain$background_adjusted_effect / expected$se), expected$df))
  fit$background_rank <- 1
  rank_one <- results(fit)
  expect_equal(rank_one$rank_adjusted_total_se, moderated$rank_adjusted_total_se)
  expect_equal(rank_one$rank_adjusted_total_df_for_test,
               moderated$rank_adjusted_total_df_for_test)
})
