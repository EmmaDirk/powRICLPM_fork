test_that("all columns from summary.powRICLMP(...) are named", {
  lagged_effects <- matrix(c(0.4, 0.15, 0.2, 0.3), ncol = 2, byrow = TRUE)
  out <- powRICLPM(
    target_power = 0.8,
    sample_size = c(500),
    time_points = c(3),
    ICC = c(0.4),
    reliability = c(1),
    RI_cor = 0.3,
    lagged_effects = lagged_effects,
    within_cor = -0.42,
    reps = 2,
    seed = 1234
  )

  table_parameter <- summary(out, parameter = "wB2~wA1")
  expect_equal(colnames(table_parameter), c("Sample size", "Time points", "ICC", "Reliability", "Population", "Avg", "Bias", "Min", "EmpSE", "SEAvg", "MSE", "Accuracy", "Cover", "Power", "Error", "Not converged", "Inadmissible"))

  table_condition <- summary(out, sample_size = 500, intraclass_correlation = 0.4, time_points = 3, reliability = 1)
  expect_equal(colnames(table_condition), c("Population", "Avg", "Bias", "Min", "EmpSE", "SEAvg", "MSE", "Accuracy", "Cover", "Power"))

  table_condition_inferred_reliability <- summary(out, sample_size = 500, intraclass_correlation = 0.4, time_points = 3)
  expect_equal(table_condition_inferred_reliability, table_condition)

  table_condition_icc <- summary(out, sample_size = 500, ICC = 0.4, time_points = 3, reliability = 1)
  expect_equal(table_condition_icc, table_condition)

  condition_summary_output <- capture.output(
    summary(out, sample_size = 500, ICC = 0.4, time_points = 3, reliability = 1)
  )
  expect_false(any(grepl("Errors:", condition_summary_output, fixed = TRUE)))

  out_with_error <- out
  out_with_error$conditions[[1]]$estimation_information$n_error <- 1
  condition_summary_error_output <- capture.output(
    summary(out_with_error, sample_size = 500, ICC = 0.4, time_points = 3, reliability = 1)
  )
  expect_true(any(grepl("Errors:", condition_summary_error_output, fixed = TRUE)))

  expect_error(
    summary(out, sample_size = 500, intraclass_correlation = 0.4, ICC = 0.4, time_points = 3, reliability = 1),
    "Both.*intraclass_correlation.*ICC"
  )
})

test_that("summary.powRICLPM labels intraclass correlation from the selected argument name", {
  lagged_effects <- matrix(c(0.4, 0.15, 0.2, 0.3), ncol = 2, byrow = TRUE)
  out <- powRICLPM(
    target_power = 0.8,
    sample_size = 500,
    time_points = 3,
    intraclass_correlation = 0.4,
    RI_cor = 0.3,
    lagged_effects = lagged_effects,
    within_cor = -0.42,
    reps = 2,
    seed = 1234
  )

  table_parameter <- summary(out, parameter = "wB2~wA1")
  expect_equal(colnames(table_parameter)[3], "Intraclass correlation")
})

test_that("summary.powRICLPM handles time-varying reliability condition labels", {
  lagged_effects <- matrix(c(0.4, 0.15, 0.2, 0.3), ncol = 2, byrow = TRUE)
  out <- powRICLPM(
    target_power = 0.8,
    sample_size = 500,
    time_points = 3,
    intraclass_correlation = 0.4,
    reliability = c(0.8, 0.9, 1),
    RI_cor = 0.3,
    lagged_effects = lagged_effects,
    within_cor = -0.42,
    reps = 2,
    seed = 1234
  )

  expect_equal(give(out, "conditions")$reliability, "time-varying")

  table_condition <- summary(
    out,
    sample_size = 500,
    intraclass_correlation = 0.4,
    time_points = 3
  )
  expect_equal(colnames(table_condition), c("Population", "Avg", "Bias", "Min", "EmpSE", "SEAvg", "MSE", "Accuracy", "Cover", "Power"))
})
