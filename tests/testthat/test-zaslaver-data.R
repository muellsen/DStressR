test_that("Zaslaver data preserve source dimensions and unique reporter IDs", {
  data_env <- new.env(parent = emptyenv())
  suppressWarnings(utils::data(
    list = c(
      "zaslaver_promoter_activity",
      "zaslaver_promoter_timecourse",
      "zaslaver_promoter_annotations"
    ),
    package = "DStressR",
    envir = data_env
  ))
  if (!exists("zaslaver_promoter_activity", envir = data_env, inherits = FALSE)) {
    load(test_path("../../data/zaslaver_promoter_activity.rda"), envir = data_env)
    load(test_path("../../data/zaslaver_promoter_timecourse.rda"), envir = data_env)
    load(test_path("../../data/zaslaver_promoter_annotations.rda"), envir = data_env)
  }

  zaslaver_promoter_activity <- data_env$zaslaver_promoter_activity
  zaslaver_promoter_timecourse <- data_env$zaslaver_promoter_timecourse
  zaslaver_promoter_annotations <- data_env$zaslaver_promoter_annotations

  expect_equal(nrow(zaslaver_promoter_activity), 23040)
  expect_equal(nrow(zaslaver_promoter_timecourse), 622080)
  expect_equal(nrow(zaslaver_promoter_annotations), 1920)

  expect_equal(length(unique(zaslaver_promoter_activity$reporter_id)), 1920)
  expect_equal(length(unique(zaslaver_promoter_timecourse$reporter_id)), 1920)
  expect_equal(length(unique(zaslaver_promoter_annotations$reporter_id)), 1920)

  expect_setequal(
    unique(zaslaver_promoter_activity$condition_label),
    c("Glucose", "no Glucose", "no AA", "Phosphate limited", "Nitrogen limited", "Ethanol")
  )
  expect_setequal(unique(zaslaver_promoter_activity$growth_rate), c(0.25, 0.8))
  expect_equal(length(unique(zaslaver_promoter_timecourse$time_index)), 54)
})
