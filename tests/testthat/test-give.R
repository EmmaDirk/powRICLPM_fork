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
  df_condition_icc <- give(out1, "ICC")
  df_problems <- give(out1, "estimation_problems")
  df_results <- give(out1, what = "results", parameter = "wB2~wA1")
  df_uncertainty <- give(out1, what = "uncertainty", parameter = "wB2~wA1")
  df_names <- give(out1, "names")

  # Run tests
  expect_error(give(1, "conditions"))
  expect_error(give(out1, "results"))
  expect_error(give(out1, "AF_proportion"), "only available for DPM")

  expect_s3_class(df_conditions, "data.frame")
  expect_equal(dim(df_conditions), c(2, 4))
  expect_equal(names(df_conditions)[3], "ICC")
  expect_equal(names(df_condition_alias)[3], "intraclass_correlation")
  expect_equal(names(df_condition_icc)[3], "ICC")
  names(df_condition_alias) <- names(df_conditions)
  expect_equal(df_condition_alias, df_conditions)
  expect_equal(df_condition_icc, df_conditions)

  expect_s3_class(df_problems, "data.frame")
  expect_equal(dim(df_problems), c(2, 7))

  expect_s3_class(df_results, "data.frame")
  expect_equal(dim(df_results), c(2, 14))

  expect_s3_class(df_uncertainty, "data.frame")
  expect_equal(nrow(df_uncertainty), 2)
  expect_true(all(c("MCSE_average", "MCSE_EmpSE", "MCSE_accuracy") %in% names(df_uncertainty)))
  expect_true("MCSE_SD" %in% names(df_uncertainty))

  expect_type(df_names, "character")
  expect_equal(length(df_names), 20)
})

test_that("give() labels DPM proportion conditions", {
  object <- list(
    conditions = list(list(
      sample_size = 1000,
      time_points = 3,
      ICC = 0.2,
      reliability = 1
    )),
    session = list(
      model = "DPM",
      version = "0.2.1",
      argument_names = list(intraclass_correlation = "AF_proportion")
    )
  )
  class(object) <- c("powRICLPM", "list")

  expect_equal(names(give(object, "conditions"))[3], "AF_proportion")
  expect_equal(names(give(object, "AF_proportion"))[3], "AF_proportion")
  expect_error(give(object, "intraclass_correlation"), "not available for DPM")
  expect_error(give(object, "ICC"), "not available for DPM")
  expect_false("reliability" %in% names(give(object, "conditions")))
  expect_error(give(object, "reliability"), "not available for DPM")

  print_output <- capture.output(print(object))
  expect_true(any(grepl("Dynamic Panel Model \\(DPM\\)", print_output)))
  expect_false(any(grepl("Reliability", print_output, fixed = TRUE)))
})

test_that("give() notes free-loading anchor for condition tables", {
  object_dpm <- list(
    conditions = list(list(
      sample_size = 1000,
      time_points = 3,
      AF_proportion = 0.2,
      reliability = 1,
      constraints = "AF_loadings_free"
    )),
    session = list(
      model = "DPM",
      argument_names = list(AF_proportion = "AF_proportion")
    )
  )
  class(object_dpm) <- c("powDPM", "powRICLPM", "list")

  object_riclpm <- list(
    conditions = list(list(
      sample_size = 1000,
      time_points = 3,
      ICC = 0.5,
      reliability = 1,
      constraints = "RI_loadings_free"
    )),
    session = list(
      model = "RICLPM",
      argument_names = list(intraclass_correlation = "ICC")
    )
  )
  class(object_riclpm) <- c("powRICLPM", "list")

  expect_silent(give(object_dpm, "conditions"))
  expect_silent(give(object_riclpm, "conditions"))
  expect_silent(inote_condition_loading_anchor(object_dpm, "AF_A=~A3"))
})
