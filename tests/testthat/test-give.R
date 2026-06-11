test_that("give() works", {

  # Create valid powRICLPM() input
  lagged_effects <- matrix(c(0.4, 0.15, 0.2, 0.3), ncol = 2, byrow = TRUE)

  out1 <- powRICLPM(
    target_power = 0.8,
    sample_size = c(300, 400),
    time_points = 3,
    ICC = 0.5,
    RI_cor = 0.3,
    lagged_effects = lagged_effects,
    within_cor = 0.3,
    reps = 2,
    seed = 123456
  )

  # Execute give()
  df_conditions <- give(out1, "conditions")
  df_condition_alias <- give(out1, "intraclass_correlation")
  expect_message(
    df_condition_legacy <- give(out1, "ICC"),
    "Please consider using.*intraclass_correlation.*instead of.*ICC"
  )
  df_problems <- give(out1, "estimation_problems")
  df_results <- give(out1, what = "results", parameter = "wB2~wA1")
  df_uncertainty <- give(out1, what = "uncertainty", parameter = "wB2~wA1")
  df_names <- give(out1, "names")

  # Run tests
  expect_error(give(1, "conditions"))
  expect_error(give(out1, "results"))

  expect_s3_class(df_conditions, "data.frame")
  expect_equal(dim(df_conditions), c(2, 4))
  expect_equal(df_condition_alias, df_conditions)
  expect_equal(df_condition_legacy, df_conditions)

  expect_s3_class(df_problems, "data.frame")
  expect_equal(dim(df_problems), c(2, 7))

  expect_s3_class(df_results, "data.frame")
  expect_equal(dim(df_results), c(2, 14))

  expect_s3_class(df_uncertainty, "data.frame")
  expect_equal(nrow(df_uncertainty), 2)

  expect_type(df_names, "character")
  expect_equal(length(df_names), 20)
})
