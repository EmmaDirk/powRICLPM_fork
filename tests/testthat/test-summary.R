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
  expect_equal(colnames(table_parameter), c("Condition", "Population", "Avg", "Bias", "Min", "EmpSE", "SEAvg", "MSE", "Accuracy", "Cover", "Power"))

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

  overview <- summary(out)
  expect_equal(colnames(overview)[4], "Intraclass correlation")
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
  expect_output(summary(
    object,
    sample_size = 600,
    time_points = 4,
    ICC = 0.5,
    reliability = c(A = 0.8, B = 0.7)
  ), "SIMULATION RESULTS")
  expect_error(
    summary(object, sample_size = 600, time_points = 4, ICC = 0.5, reliability = "A=0.8, B=0.8"),
    "A = 0.8, B = 0.7"
  )
})

test_that("matrix reliability collapses when variables match in every condition", {
  base_estimates <- data.frame(
    parameter = "A1~~A1",
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
  base_condition <- list(
    sample_size = 500,
    time_points = 5,
    ICC = 0.4,
    reliability = "A = 0.7, B = 0.7",
    reliability_matrix = matrix(0.7, nrow = 2, ncol = 5),
    loadings = matrix(1, nrow = 2, ncol = 5),
    estimates = base_estimates,
    MCSEs = data.frame(
      MCSE_average = 0.01,
      MCSE_bias = 0.01,
      MCSE_MSE = 0.01,
      MCSE_coverage = 0.01,
      MCSE_SEAvg = 0.01,
      MCSE_EmpSE = 0.01,
      MCSE_SD = 0.01,
      MCSE_accuracy = 0.01,
      MCSE_power = 0.01
    ),
    estimation_information = list(
      n_error = 0,
      n_nonconvergence = 0,
      n_inadmissible = 0,
      n_completed = 2
    ),
    skewness = 0,
    kurtosis = 0,
    significance_criterion = 0.05
  )
  object <- list(
    conditions = list(
      base_condition,
      modifyList(base_condition, list(
        reliability = "A = 0.8, B = 0.8",
        reliability_matrix = matrix(0.8, nrow = 2, ncol = 5)
      )),
      modifyList(base_condition, list(
        reliability = "A = 0.9, B = 0.9",
        reliability_matrix = matrix(0.9, nrow = 2, ncol = 5)
      ))
    ),
    session = list(
      model = "RICLPM",
      target_power = 0.8,
      version = "0.2.1",
      reps = 2,
      constraints = "ME",
      bounds = FALSE,
      estimate_ME = TRUE,
      argument_names = list(intraclass_correlation = "ICC")
    )
  )
  class(object) <- c("powRICLPM", "list")

  conditions <- give(object, "conditions")
  expect_equal(names(conditions), c("sample_size", "time_points", "ICC", "reliability"))
  expect_equal(conditions$reliability, c("0.7", "0.8", "0.9"))

  overview <- summary(object)
  expect_equal(colnames(overview)[5], "Reliability")
  expect_false(any(grepl("Reliability A|Reliability B", colnames(overview))))

  parameter_summary <- summary(object, parameter = "A1~~A1")
  expect_equal(colnames(parameter_summary)[1], "Condition")
  expect_false("Reliability" %in% colnames(parameter_summary))

  selected <- summary(
    object,
    sample_size = 500,
    time_points = 5,
    ICC = 0.4,
    reliability = 0.8
  )
  expect_equal(rownames(selected), "A1~~A1")
})

test_that("mixed matrix reliability conditions keep stable condition columns", {
  base_estimates <- data.frame(
    parameter = "A1~~A1",
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
  base_condition <- list(
    sample_size = 500,
    time_points = 5,
    ICC = 0.4,
    reliability = "A = 0.8, B = 0.8",
    reliability_matrix = matrix(0.8, nrow = 2, ncol = 5),
    loadings = matrix(1, nrow = 2, ncol = 5),
    estimates = base_estimates,
    MCSEs = data.frame(
      MCSE_average = 0.01,
      MCSE_bias = 0.01,
      MCSE_MSE = 0.01,
      MCSE_coverage = 0.01,
      MCSE_SEAvg = 0.01,
      MCSE_EmpSE = 0.01,
      MCSE_SD = 0.01,
      MCSE_accuracy = 0.01,
      MCSE_power = 0.01
    ),
    estimation_information = list(
      n_error = 0,
      n_nonconvergence = 0,
      n_inadmissible = 0,
      n_completed = 2
    ),
    skewness = 0,
    kurtosis = 0,
    significance_criterion = 0.05
  )
  object <- list(
    conditions = list(
      base_condition,
      modifyList(base_condition, list(
        reliability = "A = 0.8, B = 0.7",
        reliability_matrix = matrix(c(0.8, 0.7), nrow = 2, ncol = 5)
      ))
    ),
    session = list(
      model = "RICLPM",
      target_power = 0.8,
      version = "0.2.1",
      reps = 2,
      constraints = "ME",
      bounds = FALSE,
      estimate_ME = TRUE,
      argument_names = list(intraclass_correlation = "ICC")
    )
  )
  class(object) <- c("powRICLPM", "list")

  conditions <- give(object, "conditions")
  expect_equal(names(conditions), c("sample_size", "time_points", "ICC", "reliability_A", "reliability_B"))
  expect_equal(conditions$reliability_A, c("0.8", "0.8"))
  expect_equal(conditions$reliability_B, c("0.8", "0.7"))

  parameter_summary <- summary(object, parameter = "A1~~A1")
  expect_equal(colnames(parameter_summary)[1], "Condition")
  expect_false(any(grepl("Reliability A|Reliability B", colnames(parameter_summary))))

  overview <- summary(object)
  expect_equal(colnames(overview)[5:6], c("Reliability A", "Reliability B"))
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
  expect_equal(colnames(table_parameter)[1], "Condition")

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

test_that("DPM custom loadings stay with condition columns in summaries", {
  object <- list(
    conditions = list(list(
      sample_size = 800,
      time_points = 4,
      ICC = 0.2,
      reliability = 1,
      loadings = matrix(
        c(1, 0.8, 1.1,
          1, 0.8, 1.1),
        nrow = 2,
        byrow = TRUE
      ),
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
        n_error = 1,
        n_nonconvergence = 0,
        n_inadmissible = 0,
        n_completed = 0
      )
    )),
    session = list(
      model = "DPM",
      version = "0.2.1",
      reps = 1,
      bounds = FALSE,
      constraints = "AF_loadings_free",
      estimate_ME = FALSE,
      argument_names = list(intraclass_correlation = "AF_proportion")
    )
  )
  class(object) <- c("powRICLPM", "list")

  estimation_problems <- give(object, "estimation_problems")
  expect_equal(
    names(estimation_problems),
    c("sample_size", "time_points", "AF_proportion", "loadings", "errors", "not_converged", "inadmissible")
  )
  expect_equal(estimation_problems$loadings, "c(1, 0.8, 1.1)")
  expect_equal(estimation_problems$errors, 1)
  expect_equal(estimation_problems$inadmissible, 0)

  overview <- NULL
  capture.output(overview <- summary(object))
  expect_equal(
    colnames(overview),
    c("Condition", "Sample size", "Time points", "AF proportion", "Loadings", "Error", "Not converged", "Inadmissible")
  )
  expect_equal(overview$Condition, 1)
  expect_equal(overview$Loadings, "c(1, 0.8, 1.1)")
  expect_equal(overview$Error, 1)
  expect_equal(overview$Inadmissible, 0)

  parameter_summary <- NULL
  capture.output(parameter_summary <- summary(object, parameter = "B2~A1"))
  expect_equal(
    colnames(parameter_summary),
    c("Condition", "Population", "Avg", "Bias", "Min", "EmpSE", "SEAvg", "MSE", "Accuracy", "Cover", "Power")
  )
  expect_equal(parameter_summary$Condition, 1)
})
