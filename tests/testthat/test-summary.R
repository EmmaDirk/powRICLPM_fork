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
  expect_error(
    summary(out, sample_size = 500, AF_proportion = 0.4, time_points = 3),
    "only available for DPM"
  )
  expect_error(
    summary(out, parameter = c("wB2~wA1", "wA2~wB1")),
    "length 2"
  )
  expect_error(
    summary(out, parameter = 1),
    "type double"
  )
  expect_error(
    summary(out, parameter = "not_a_parameter"),
    "not found"
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

test_that("summary.powRICLPM handles vector reliability conditions", {
  lagged_effects <- matrix(c(0.4, 0.15, 0.2, 0.3), ncol = 2, byrow = TRUE)
  out <- suppressWarnings(
    powRICLPM(
      target_power = 0.8,
      sample_size = 500,
      time_points = 4,
      intraclass_correlation = 0.4,
      reliability = c(0.8, 0.9, 1),
      estimate_ME = TRUE,
      RI_cor = 0.3,
      lagged_effects = lagged_effects,
      within_cor = -0.42,
      reps = 2,
      seed = 1234
    )
  )

  expect_equal(give(out, "conditions")$reliability, c(0.8, 0.9, 1))

  table_condition <- summary(
    out,
    sample_size = 500,
    intraclass_correlation = 0.4,
    time_points = 4,
    reliability = 0.9
  )
  expect_equal(colnames(table_condition), c("Population", "Avg", "Bias", "Min", "EmpSE", "SEAvg", "MSE", "Accuracy", "Cover", "Power"))
})

test_that("matrix reliability labels print and filter consistently", {
  base_estimates <- data.frame(
    parameter = "wB2~wA1",
    population_value = 0.2,
    average = 0.2,
    bias = 0,
    minimum = 0.1,
    EmpSE = 0.05,
    SEAvg = 0.05,
    MSE = 0.0025,
    accuracy = 0.2,
    coverage = 0.95,
    power = 0.8
  )
  condition <- list(
    sample_size = 600,
    time_points = 4,
    ICC = 0.5,
    reliability = "A 0.8, B 0.7",
    reliability_matrix = matrix(c(0.8, 0.7), nrow = 2, ncol = 4),
    loadings = matrix(1, nrow = 2, ncol = 4),
    estimates = base_estimates,
    estimation_information = list(
      n_error = 0,
      n_nonconvergence = 0,
      n_inadmissible = 0,
      n_completed = 2
    ),
    skewness = 0,
    kurtosis = 0,
    constraints = "none",
    bounds = FALSE,
    estimate_ME = TRUE,
    significance_criterion = 0.05
  )
  object <- list(
    conditions = list(
      condition,
      modifyList(condition, list(
        reliability = "A 0.9, B 0.85",
        reliability_matrix = matrix(c(0.9, 0.85), nrow = 2, ncol = 4)
      ))
    ),
    session = list(
      model = "RICLPM",
      version = "0.2.1",
      reps = 2,
      constraints = "none",
      bounds = FALSE,
      estimate_ME = TRUE,
      argument_names = list(intraclass_correlation = "ICC")
    )
  )
  class(object) <- c("powRICLPM", "list")

  condition_table <- give(object, "conditions")
  expect_equal(condition_table$reliability_A, c("0.8", "0.9"))
  expect_equal(condition_table$reliability_B, c("0.7", "0.85"))
  expect_output(summary(
    object,
    sample_size = 600,
    time_points = 4,
    ICC = 0.5,
    reliability = "A = 0.8, B = 0.7"
  ), "SIMULATION RESULTS")
  expect_error(
    summary(object, sample_size = 600, time_points = 4, ICC = 0.5, reliability = "A=0.8, B=0.8"),
    "A 0.8, B 0.7"
  )
})

test_that("parameter-not-found and no-usable-estimates errors differ", {
  good_condition <- list(
    sample_size = 600,
    time_points = 4,
    ICC = 0.5,
    reliability = 1,
    estimates = data.frame(
      parameter = "wB2~wA1",
      population_value = 0.2,
      average = 0.2,
      bias = 0,
      minimum = 0.1,
      EmpSE = 0.05,
      SEAvg = 0.05,
      MSE = 0.0025,
      accuracy = 0.2,
      coverage = 0.95,
      power = 0.8
    ),
    MCSEs = data.frame(MCSE_average = 0),
    estimation_information = list(
      n_error = 0,
      n_nonconvergence = 0,
      n_inadmissible = 0,
      n_completed = 1
    )
  )
  failed_condition <- modifyList(good_condition, list(
    reliability = 0.8,
    estimates = NA,
    est_tab = data.frame(
      lhs = "wB2",
      op = "~",
      rhs = "wA1",
      free = TRUE
    ),
    estimation_information = list(
      n_error = 1,
      n_nonconvergence = 0,
      n_inadmissible = 0,
      n_completed = 0
    )
  ))
  object <- list(
    conditions = list(good_condition, failed_condition),
    session = list(
      model = "RICLPM",
      version = "0.2.1",
      reps = 1,
      estimate_ME = TRUE,
      argument_names = list(intraclass_correlation = "ICC")
    )
  )
  class(object) <- c("powRICLPM", "list")

  expect_error(give(object, "results", parameter = "not_a_parameter"), "not found")
  expect_error(give(object, "results", parameter = "wB2~wA1"), "No usable estimates")
})

test_that("summary.powRICLPM omits reliability column for DPM objects", {
  object <- list(
    conditions = list(list(
      sample_size = 800,
      time_points = 3,
      ICC = 0.2,
      reliability = 1,
      estimates = data.frame(
        parameter = "B2~A1",
        population_value = 0.2,
        average = 0.2,
        bias = 0,
        minimum = 0.1,
        EmpSE = 0.05,
        SEAvg = 0.05,
        MSE = 0.0025,
        accuracy = 0.2,
        coverage = 0.95,
        power = 0.8
      ),
      estimation_information = list(
        n_error = 0,
        n_nonconvergence = 0,
        n_inadmissible = 0,
        n_completed = 1
      )
    )),
    session = list(
      model = "DPM",
      version = "0.2.1",
      reps = 1,
      bounds = FALSE,
      constraints = "none",
      estimate_ME = FALSE,
      argument_names = list(intraclass_correlation = "AF_proportion")
    )
  )
  class(object) <- c("powRICLPM", "list")

  table_parameter <- summary(object, parameter = "B2~A1")
  expect_false("Reliability" %in% colnames(table_parameter))
  expect_false("Loadings" %in% colnames(table_parameter))
  expect_equal(colnames(table_parameter)[3], "AF proportion")

  overview_output <- capture.output(summary(object))
  expect_true(any(grepl("Dynamic Panel Model \\(DPM\\)", overview_output)))
  expect_false(any(grepl("Reliability", overview_output, fixed = TRUE)))
  expect_true(any(grepl("AF proportion", overview_output, fixed = TRUE)))
  expect_error(summary(object, reliability = 1), "not available for DPM")
  expect_error(
    summary(object, sample_size = 800, time_points = 3, ICC = 0.2),
    "not available for DPM"
  )
})
