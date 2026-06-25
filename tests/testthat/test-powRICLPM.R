test_that("basic power analysis using lavaan runs", {
  out1 <- powRICLPM(
    target_power = 0.8,
    sample_size = 1000,
    time_points = 3,
    intraclass_correlation = 0.5,
    RI_cor = 0.3,
    lagged_effects = matrix(c(0.4, 0.15, 0.2, 0.3), ncol = 2, byrow = TRUE),
    within_cor = 0.3,
    reps = 2,
    seed = 123456
  )

  expect_equal(class(out1), c("powRICLPM", "list"))
  expect_equal(names(out1), c("conditions", "session"))
  expect_equal(length(out1$conditions), 1)
  expect_equal(
    c(
      "sample_size", "time_points", "ICC", "reliability", "RI_var", "RI_cov",
      "reliability_matrix", "pop_synt", "pop_tab", "est_synt", "est_tab", "estimate_ME", "skewness",
      "kurtosis", "significance_criterion", "estimates", "MCSEs", "reps",
      "condition_id", "estimation_information"
    ) %in% names(out1$conditions[[1]]),
    rep(TRUE, times = 20)
  )
  expect_type(out1$conditions[[1]]$estimates, "list")
  expect_type(out1$conditions[[1]]$MCSEs, "list")
  expect_type(out1$conditions[[1]]$estimation_information, "list")
  expect_true(all(c("MCSE_average", "MCSE_EmpSE", "MCSE_accuracy") %in% names(out1$conditions[[1]]$MCSEs)))
  expect_true("MCSE_SD" %in% names(out1$conditions[[1]]$MCSEs))
  expect_equal(out1$conditions[[1]]$MCSEs$MCSE_SD, out1$conditions[[1]]$MCSEs$MCSE_EmpSE)

  test_summary_condition <- summary(out1, sample_size = 1000, time_points = 3, intraclass_correlation = 0.5, reliability = 1)

  expect_equal(
    test_summary_condition$Population,
    c(1.000, 1.000, 0.300, 0.400, 0.150, 0.200, 0.300, 0.400, 0.150, 0.200, 0.300, 1.000, 1.000, 0.300, 0.781, 0.781, 0.834, 0.834, 0.130, 0.130)
  )
})

test_that("powRICLPM defaults to MLR for nonnormal lavaan data", {
  lagged_effects <- matrix(c(0.4, 0.15, 0.2, 0.3), ncol = 2, byrow = TRUE)
  out <- suppressWarnings(
    powRICLPM(
      target_power = 0.8,
      sample_size = 1000,
      time_points = 3,
      intraclass_correlation = 0.5,
      RI_cor = 0.3,
      lagged_effects = lagged_effects,
      within_cor = 0.3,
      skewness = 1,
      reps = 1,
      seed = 123456
    )
  )

  expect_equal(out$session$estimator, "MLR")

  out_explicit <- suppressWarnings(
    powRICLPM(
      target_power = 0.8,
      sample_size = 1000,
      time_points = 3,
      intraclass_correlation = 0.5,
      RI_cor = 0.3,
      lagged_effects = lagged_effects,
      within_cor = 0.3,
      skewness = 1,
      estimator = "ML",
      reps = 1,
      seed = 123456
    )
  )

  expect_equal(out_explicit$session$estimator, "ML")
})

test_that("lavaan supplied loadings require freed public estimation path", {
  lagged_effects <- matrix(c(0.4, 0.15, 0.2, 0.3), ncol = 2, byrow = TRUE)

  out_default <- suppressWarnings(
    powRICLPM(
      target_power = 0.8,
      sample_size = 1000,
      time_points = 3,
      ICC = 0.5,
      RI_cor = 0.3,
      lagged_effects = lagged_effects,
      within_cor = 0.3,
      reps = 1,
      seed = 123456
    )
  )

  expect_equal(class(out_default), c("powRICLPM", "list"))

  expect_error(
    powRICLPM(
      target_power = 0.8,
      sample_size = 1000,
      time_points = 3,
      ICC = 0.5,
      RI_cor = 0.3,
      lagged_effects = lagged_effects,
      within_cor = 0.3,
      loadings = c(1, 1, 1),
      reps = 1,
      seed = 123456
    ),
    "RI_loadings_free"
  )

  expect_error(
    powRICLPM(
      target_power = 0.8,
      sample_size = 1000,
      time_points = 3,
      ICC = 0.5,
      RI_cor = 0.3,
      lagged_effects = lagged_effects,
      within_cor = 0.3,
      loadings = matrix(1, nrow = 1, ncol = 3),
      constraints = "RI_loadings_free",
      reps = 1,
      seed = 123456
    ),
    "must have 2 rows.*has 1 row"
  )

  expect_error(
    powRICLPM(
      target_power = 0.8,
      sample_size = 1000,
      time_points = 3,
      ICC = 0.5,
      RI_cor = 0.3,
      lagged_effects = lagged_effects,
      within_cor = 0.3,
      loadings = matrix(1, nrow = 2, ncol = 4),
      constraints = "RI_loadings_free",
      reps = 1,
      seed = 123456
    ),
    "one column per time point.*has 4 columns.*time_points.*= 3"
  )

  expect_error(
    powRICLPM(
      target_power = 0.8,
      sample_size = 1000,
      time_points = Inf,
      ICC = 0.5,
      RI_cor = 0.3,
      lagged_effects = lagged_effects,
      within_cor = 0.3,
      loadings = c(1, 1, 1),
      constraints = "RI_loadings_free",
      reps = 1,
      seed = 123456
    ),
    "finite"
  )
})

test_that("lavaan custom loadings require freed public estimation path", {
  lagged_effects <- matrix(c(0.4, 0.15, 0.2, 0.3), ncol = 2, byrow = TRUE)
  loadings <- matrix(c(1, 0, -1.2, 1, 2.5, 0.25), nrow = 2, byrow = TRUE)

  expect_error(
    powRICLPM(
      target_power = 0.8,
      sample_size = 1000,
      time_points = 3,
      ICC = 0.5,
      RI_cor = 0.3,
      lagged_effects = lagged_effects,
      within_cor = 0.3,
      loadings = loadings,
      reps = 1,
      seed = 123456
    ),
    "RI_loadings_free"
  )

  out_free <- suppressWarnings(
    powRICLPM(
      target_power = 0.8,
      sample_size = 1000,
      time_points = 3,
      ICC = 0.5,
      RI_cor = 0.3,
      lagged_effects = lagged_effects,
      within_cor = 0.3,
      loadings = loadings,
      constraints = "RI_loadings_free",
      reps = 1,
      seed = 123456
    )
  )

  RI_loading_parameters <- c("RI_A=~A2", "RI_A=~A3", "RI_B=~B2", "RI_B=~B3")
  index <- match(RI_loading_parameters, out_free$conditions[[1]]$estimates$parameter)

  expect_true(grepl("RI_A=~0*A2", out_free$conditions[[1]]$pop_synt, fixed = TRUE))
  expect_true(grepl("RI_A=~-1.2*A3", out_free$conditions[[1]]$pop_synt, fixed = TRUE))
  expect_true(grepl("RI_B=~2.5*B2", out_free$conditions[[1]]$pop_synt, fixed = TRUE))
  expect_true(grepl("RI_B=~0.25*B3", out_free$conditions[[1]]$pop_synt, fixed = TRUE))
  expect_false(anyNA(index))
  expect_equal(out_free$conditions[[1]]$estimates$population_value[index], c(0, -1.2, 2.5, 0.25))
})

test_that("ICOV no convergence warnings are recognized", {
  expect_true(is_icov_nonconvergence_warning("lavaan->getICOV():\n   no convergence"))
  expect_false(is_icov_nonconvergence_warning("lavaan->getICOV():"))
  expect_false(is_icov_nonconvergence_warning("lavaan->lav_object_post_check(): some estimated lv variances are negative"))
})

test_that("fatal estimation errors are counted separately from other estimation issues", {
  estimation_information <- count_estimation_information(
    reps_completed = 2,
    n_error = 1,
    converged = list(TRUE, TRUE),
    admissible = list(TRUE, TRUE),
    n_data_icov_nonconvergence = 0,
    post_check_icov_nonconverged = c(FALSE, FALSE),
    n_icov_nonconvergence = 0
  )

  expect_equal(estimation_information$n_completed, 2)
  expect_equal(estimation_information$n_error, 1)
  expect_equal(estimation_information$n_nonconvergence, 0)
  expect_equal(estimation_information$n_inadmissible, 0)
})

test_that("ICOV failures remain counted as nonconverged estimation issues", {
  estimation_information <- count_estimation_information(
    reps_completed = 0,
    n_error = 0,
    converged = list(),
    admissible = list(),
    n_data_icov_nonconvergence = 1,
    post_check_icov_nonconverged = logical(),
    n_icov_nonconvergence = 1
  )

  expect_equal(estimation_information$n_completed, 0)
  expect_equal(estimation_information$n_error, 0)
  expect_equal(estimation_information$n_nonconvergence, 2)
  expect_equal(estimation_information$n_inadmissible, 0)
})

test_that("ICOV failures with freed RI loadings are counted as nonconverged", {
  lagged_effects <- matrix(c(0.4, 0.15, 0.2, 0.3), ncol = 2, byrow = TRUE)
  loadings <- matrix(c(
    1, 0, -1.2, 2,
    1, 2.5, 0.25, -0.5
  ), nrow = 2, byrow = TRUE)

  icov_warnings <- character()
  out <- withCallingHandlers(
    powRICLPM(
      target_power = 0.8,
      sample_size = 1000,
      time_points = 4,
      ICC = 0.5,
      RI_cor = 0.3,
      lagged_effects = lagged_effects,
      within_cor = 0.3,
      loadings = loadings,
      constraints = c("within", "RI_loadings_free"),
      reps = 2,
      seed = 123456
    ),
    warning = function(w) {
      if (is_icov_nonconvergence_warning(conditionMessage(w))) {
        icov_warnings <<- c(icov_warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    }
  )

  estimation_information <- out$conditions[[1]]$estimation_information

  expect_length(icov_warnings, 2)
  expect_equal(estimation_information$n_completed, 0)
  expect_equal(estimation_information$n_error, 0)
  expect_equal(estimation_information$n_nonconvergence, 2)
  expect_equal(estimation_information$n_inadmissible, 0)
  expect_true(all(is.na(out$conditions[[1]]$estimates$average)))
  expect_true(all(is.na(out$conditions[[1]]$estimates$SEAvg)))
  expect_true(all(is.na(out$conditions[[1]]$estimates$power)))
})

test_that("ICC and Phi argument names remain supported", {
  out <- powRICLPM(
    target_power = 0.8,
    sample_size = 1000,
    time_points = 3,
    ICC = 0.5,
    RI_cor = 0.3,
    lagged_effects = matrix(c(0.4, 0.15, 0.2, 0.3), ncol = 2, byrow = TRUE),
    within_cor = 0.3,
    reps = 1,
    seed = 123456
  )

  out_phi <- powRICLPM(
    target_power = 0.8,
    sample_size = 1000,
    time_points = 3,
    intraclass_correlation = 0.5,
    RI_cor = 0.3,
    Phi = matrix(c(0.4, 0.15, 0.2, 0.3), ncol = 2, byrow = TRUE),
    within_cor = 0.3,
    reps = 1,
    seed = 123456
  )

  expect_equal(out$conditions[[1]]$ICC, 0.5)
  expect_equal(out$session$argument_names$intraclass_correlation, "ICC")
  expect_equal(out_phi$conditions[[1]]$ICC, 0.5)
  expect_equal(out_phi$session$argument_names$lagged_effects, "Phi")
  expect_equal(out_phi$conditions[[1]]$estimates$population_value[out_phi$conditions[[1]]$estimates$parameter == "wB2~wA1"], 0.2)
})

test_that("powRICLPM requires sample_size or a complete search range", {
  lagged_effects <- matrix(c(0.4, 0.15, 0.2, 0.3), ncol = 2, byrow = TRUE)

  expect_error(
    powRICLPM(
      target_power = 0.8,
      time_points = 3,
      ICC = 0.5,
      RI_cor = 0.3,
      lagged_effects = lagged_effects,
      within_cor = 0.3,
      reps = 1,
      seed = 123456
    ),
    "Either.*sample_size.*search"
  )

  out <- powRICLPM(
    target_power = 0.8,
    search_lower = 1000,
    search_upper = 1000,
    search_step = 20,
    time_points = 3,
    ICC = 0.5,
    RI_cor = 0.3,
    lagged_effects = lagged_effects,
    within_cor = 0.3,
    reps = 1,
    seed = 123456
  )

  expect_equal(out$conditions[[1]]$sample_size, 1000)
})

test_that("powRICLPM validation errors use user-facing alias names", {
  lagged_effects <- matrix(c(0.4, 0.15, 0.2, 0.3), ncol = 2, byrow = TRUE)

  expect_error(
    powRICLPM(
      target_power = 0.8,
      sample_size = 1000,
      time_points = 3,
      RI_cor = 0.3,
      lagged_effects = lagged_effects,
      within_cor = 0.3,
      reps = 1
    ),
    "intraclass_correlation"
  )
  expect_error(
    powRICLPM(
      target_power = 0.8,
      sample_size = 1000,
      time_points = 3,
      ICC = "0.5",
      RI_cor = 0.3,
      lagged_effects = lagged_effects,
      within_cor = 0.3,
      reps = 1
    ),
    "ICC"
  )
  expect_error(
    powRICLPM(
      target_power = 0.8,
      sample_size = 1000,
      time_points = 3,
      intraclass_correlation = 0.5,
      RI_cor = 0.3,
      Phi = "lagged_effects",
      within_cor = 0.3,
      reps = 1
    ),
    "Phi"
  )
})

test_that("powDPM validation errors use DPM-specific argument names", {
  lagged_effects <- matrix(c(0.4, 0.15, 0.2, 0.3), ncol = 2, byrow = TRUE)

  expect_error(
    powDPM(
      target_power = 0.8,
      sample_size = 1000,
      time_points = 3,
      AF_proportion = 0.2,
      AF_cor = 0.3,
      lagged_effects = lagged_effects,
      wave_cor = 0.3,
      reps = 1
    ),
    "at least 4 time points"
  )
  expect_error(
    powDPM(
      target_power = 0.8,
      sample_size = 1000,
      time_points = 4,
      AF_proportion = "0.2",
      AF_cor = 0.3,
      lagged_effects = lagged_effects,
      wave_cor = 0.3,
      reps = 1
    ),
    "AF_proportion"
  )
  expect_error(
    powDPM(
      target_power = 0.8,
      sample_size = 1000,
      time_points = 4,
      AF_cor = 0.3,
      lagged_effects = lagged_effects,
      wave_cor = 0.3,
      reps = 1
    ),
    "AF_proportion.*missing"
  )
  expect_error(
    powDPM(
      target_power = 0.8,
      sample_size = 1000,
      time_points = 4,
      AF_proportion = 0.2,
      lagged_effects = lagged_effects,
      wave_cor = 0.3,
      reps = 1
    ),
    "AF_cor.*missing"
  )
  expect_error(
    powDPM(
      target_power = 0.8,
      sample_size = 1000,
      time_points = 4,
      AF_proportion = 0.2,
      AF_cor = 0.3,
      lagged_effects = lagged_effects,
      reps = 1
    ),
    "wave_cor.*missing"
  )
  expect_error(
    powDPM(
      target_power = 0.8,
      sample_size = 1000,
      time_points = 4,
      AF_proportion = 0.2,
      AF_cor = 0.3,
      lagged_effects = lagged_effects,
      wave_cor = 0.3,
      reps = 1,
      constraints = "ME"
    ),
    "estimate_ME = TRUE"
  )
  expect_warning(
    expect_error(
      powDPM(
        target_power = 0.8,
        sample_size = 1000,
        time_points = 4,
        AF_proportion = 0.2,
        AF_cor = 0.3,
        lagged_effects = lagged_effects,
        wave_cor = 0.3,
        constraints = "ME",
        reps = 1
      ),
      "estimate_ME = TRUE"
    ),
    NA
  )
  expect_error(
    powDPM(
      target_power = 0.8,
      sample_size = 1000,
      time_points = 4,
      AF_proportion = 0.2,
      AF_cor = 0.3,
      lagged_effects = lagged_effects,
      wave_cor = 0.3,
      reliability = 0.8,
      estimate_ME = TRUE,
      reps = 1
    ),
    "not identified.*39 parameters.*36 distinct"
  )
  dpm_psi_error <- tryCatch(
    powDPM(
      target_power = 0.8,
      sample_size = 1000,
      time_points = 4,
      AF_proportion = 0.4,
      AF_cor = 0.2,
      lagged_effects = matrix(c(0.4, 0.1, 0.2, 0.3), ncol = 2, byrow = TRUE),
      wave_cor = 0.1,
      reps = 1
    ),
    error = function(e) e
  )
  expect_s3_class(dpm_psi_error, "error")
  expect_match(conditionMessage(dpm_psi_error), "wave_cor.*AF_proportion.*AF_cor")
  expect_false(grepl("within_cor", conditionMessage(dpm_psi_error), fixed = TRUE))

  out <- suppressWarnings(
    powDPM(
      target_power = 0.8,
      sample_size = 1000,
      time_points = 4,
      AF_proportion = c(0.1, 0.2),
      AF_cor = 0.3,
      lagged_effects = lagged_effects,
      wave_cor = 0.3,
      reps = 1,
      seed = 123456
    )
  )
  expect_equal(out$session$model, "DPM")
  expect_equal(out$conditions[[1]]$ICC, 0.1)
  expect_equal(out$conditions[[2]]$ICC, 0.2)
  expect_length(out$conditions, 2)
  expect_equal(out$session$argument_names$intraclass_correlation, "AF_proportion")
})

test_that("simple identified powDPM results feed give, summary, and plot", {
  lagged_effects <- matrix(c(0.4, 0.1, 0.2, 0.3), ncol = 2, byrow = TRUE)

  out <- suppressWarnings(
    powDPM(
      target_power = 0.8,
      sample_size = c(600, 800),
      time_points = 4,
      AF_proportion = 0.2,
      AF_cor = 0.2,
      lagged_effects = lagged_effects,
      wave_cor = 0.1,
      reps = 2,
      seed = 123456
    )
  )

  dpm_results <- give(out, "results", parameter = "B2~A1")
  dpm_uncertainty <- give(out, "uncertainty", parameter = "B2~A1")
  dpm_summary <- summary(out, parameter = "B2~A1")

  expect_false(anyNA(dpm_results[, c("SEAvg", "accuracy", "coverage", "power")]))
  expect_false(anyNA(dpm_uncertainty[, c("MCSE_SEAvg", "MCSE_accuracy", "MCSE_coverage", "MCSE_power")]))
  expect_false(anyNA(dpm_summary[, c("SEAvg", "Accuracy", "Cover", "Power")]))
  expect_warning(
    plot(out, parameter = "B2~A1"),
    NA
  )
})

test_that("powDPM supports measurement error on latent true scores", {
  lagged_effects <- matrix(c(0.4, 0.1, 0.2, 0.3), ncol = 2, byrow = TRUE)

  out <- suppressWarnings(
    powDPM(
      target_power = 0.8,
      sample_size = 900,
      time_points = 4,
      AF_proportion = 0.2,
      AF_cor = 0.2,
      lagged_effects = lagged_effects,
      wave_cor = 0.1,
      reliability = c(0.8, 0.9),
      estimate_ME = TRUE,
      constraints = c("stationarity", "ME"),
      reps = 2,
      seed = 123456
    )
  )

  expect_true(out$session$estimate_ME)
  expect_equal(give(out, "conditions")$reliability, c(0.8, 0.9))
  expect_equal(out$conditions[[1]]$ME_var, matrix(0.25, nrow = 2, ncol = 4))
  expect_true(grepl("tA1=~1*A1", out$conditions[[1]]$pop_synt, fixed = TRUE))
  expect_true(grepl("AF_A=~1*tA2", out$conditions[[1]]$pop_synt, fixed = TRUE))
  expect_true(grepl("tB2~0.2*tA1", out$conditions[[1]]$pop_synt, fixed = TRUE))
  expect_true(grepl("A1~~0.25*A1", out$conditions[[1]]$pop_synt, fixed = TRUE))
  expect_true(grepl("A1~~MEvarA*start(0.25)*A1", out$conditions[[1]]$est_synt, fixed = TRUE))

  expect_true("tB2~tA1" %in% give(out, "names"))
  dpm_me_summary <- summary(out, parameter = "tB2~tA1")
  expect_false(anyNA(dpm_me_summary[, c("SEAvg", "Accuracy", "Cover", "Power")]))
  expect_warning(
    plot(out, parameter = "tB2~tA1"),
    NA
  )
})

test_that("powDPM validation errors catch common DPM loading mistakes", {
  lagged_effects <- matrix(c(0.4, 0.15, 0.2, 0.3), ncol = 2, byrow = TRUE)
  base <- list(
    target_power = 0.8,
    sample_size = 1000,
    time_points = 4,
    AF_proportion = 0.2,
    AF_cor = 0.3,
    lagged_effects = lagged_effects,
    wave_cor = 0.3,
    reps = 1,
    seed = 123456,
    constraints = "loadings_free"
  )

  expect_error(
    do.call(powDPM, c(base, list(loadings = c(1, 0.8)))),
    "one value per time point.*first value must be `NA`"
  )
  expect_error(
    do.call(powDPM, c(base, list(loadings = c(NA, 0.8, 1, 1)))),
    "second DPM loading must be 1"
  )
  expect_error(
    do.call(powDPM, c(base, list(loadings = c(NaN, 1, 0.8, 0.7)))),
    "first DPM loading must be `NA`"
  )
  expect_error(
    do.call(powDPM, c(base, list(loadings = c(NA, 1, Inf, 0.7)))),
    "finite values, except for the required first-wave NA"
  )
  expect_error(
    do.call(powDPM, c(base, list(loadings = matrix(c(NA, 1, 0.8, 0.7, 1, 1, 1.2, 0.9), nrow = 2, byrow = TRUE)))),
    "first DPM loading column must be `NA`"
  )
})

test_that("conflicting argument aliases error", {
  expect_error(
    powRICLPM(
      target_power = 0.8,
      sample_size = 1000,
      time_points = 3,
      intraclass_correlation = 0.5,
      ICC = 0.5,
      RI_cor = 0.3,
      lagged_effects = matrix(c(0.4, 0.15, 0.2, 0.3), ncol = 2, byrow = TRUE),
      within_cor = 0.3,
      reps = 1,
      seed = 123456
    ),
    "Both.*intraclass_correlation.*ICC"
  )

  expect_error(
    powRICLPM(
      target_power = 0.8,
      sample_size = 1000,
      time_points = 3,
      intraclass_correlation = 0.5,
      RI_cor = 0.3,
      lagged_effects = matrix(c(0.4, 0.15, 0.2, 0.3), ncol = 2, byrow = TRUE),
      Phi = matrix(c(0.4, 0.15, 0.2, 0.3), ncol = 2, byrow = TRUE),
      within_cor = 0.3,
      reps = 1,
      seed = 123456
    ),
    "Both.*lagged_effects.*Phi"
  )
})

test_that("basic power analysis with multiple experimental conditions works", {
  out1 <- powRICLPM(
    target_power = 0.8,
    sample_size = c(400, 500),
    time_points = c(3, 4),
    ICC = c(0.4, 0.6),
    RI_cor = 0.3,
    lagged_effects = matrix(c(0.4, 0.15, 0.2, 0.3), ncol = 2, byrow = TRUE),
    within_cor = 0.3,
    reps = 2,
    seed = 123456
  )

  expect_equal(class(out1), c("powRICLPM", "list"))
  expect_equal(names(out1), c("conditions", "session"))
  expect_equal(length(out1$conditions), 8)
  expect_equal(
    c(
      "sample_size", "time_points", "ICC", "reliability", "RI_var", "RI_cov",
      "reliability_matrix", "pop_synt", "pop_tab", "est_synt", "est_tab", "estimate_ME", "skewness",
      "kurtosis", "significance_criterion", "estimates", "MCSEs", "reps",
      "condition_id", "estimation_information"
    ) %in% names(out1$conditions[[1]]),
    rep(TRUE, times = 20)
  )

})

test_that("power analysis with constraints works", {
  out_within <- powRICLPM(
    target_power = 0.8,
    sample_size = 1000,
    time_points = 3,
    ICC = 0.5,
    RI_cor = 0.3,
    lagged_effects = matrix(c(0.4, 0.15, 0.2, 0.3), ncol = 2, byrow = TRUE),
    within_cor = 0.3,
    reps = 2,
    seed = 123456,
    constraints = "within"
  )
  test_summary_condition <- summary(out_within, sample_size = 1000, time_points = 3, ICC = 0.5, reliability = 1)

  expect_equal(
    test_summary_condition$Population,
    c(1.000, 1.000, 0.300, 0.400, 0.150, 0.200, 0.300, 0.400, 0.150, 0.200, 0.300, 1.000, 1.000, 0.300, 0.781, 0.781, 0.834, 0.834, 0.130, 0.130)
  )

  out_stationarity <- powRICLPM(
    target_power = 0.8,
    sample_size = 1000,
    time_points = 3,
    ICC = 0.5,
    RI_cor = 0.3,
    lagged_effects = matrix(c(0.4, 0.15, 0.2, 0.3), ncol = 2, byrow = TRUE),
    within_cor = 0.3,
    reps = 2,
    seed = 123456,
    constraints = "stationarity"
  )

  test_summary_condition_stationarity <- summary(out_stationarity, sample_size = 1000, time_points = 3, ICC = 0.5, reliability = 1)

  # No starting values, and hence zero population values for residual covariances (to aid convergence)
  expect_equal(
    test_summary_condition_stationarity$Population,
    c(1.000, 1.000, 0.300, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 0.400, 0.150, 0.200, 0.300, 0.400, 0.150, 0.200, 0.300, 0.300, 0.781, 0.781, 0.834, 0.834, 0, 0)
  )
})


test_that("power analysis for the STARTS model works", {
  expect_warning({
    out <- powRICLPM(
      target_power = 0.8,
      sample_size = c(500),
      time_points = 4,
      ICC = .5,
      RI_cor = 0.3,
      lagged_effects = matrix(c(0.4, 0.15, 0.2, 0.3), ncol = 2, byrow = TRUE),
      within_cor = 0.3,
      reliability = .85,
      estimate_ME = TRUE,
      reps = 1,
      seed = 1234
    )
  })

  expect_equal(out$session$estimate_ME, TRUE)
  expect_equal(
    c("A1~~A1", "A2~~A2", "B1~~B1", "B2~~B2") %in% out$conditions[[1]]$estimates$parameter,
    c(T, T, T, T)
  )
})

test_that("reliability scenario paths run", {
  lagged_effects <- matrix(c(0.4, 0.15, 0.2, 0.3), ncol = 2, byrow = TRUE)

  out_vector <- suppressWarnings(
    powRICLPM(
      target_power = 0.8,
      sample_size = 1000,
      time_points = 3,
      ICC = 0.5,
      RI_cor = 0.3,
      lagged_effects = lagged_effects,
      within_cor = 0.3,
      reliability = c(0.8, 0.7, 1),
      reps = 1,
      seed = 123456
    )
  )

  expect_equal(length(out_vector$conditions), 3)
  expect_equal(unname(vapply(out_vector$conditions, function(x) x$reliability, numeric(1))), c(0.8, 0.7, 1))
  expect_equal(out_vector$conditions[[2]]$reliability_matrix, matrix(0.7, nrow = 2, ncol = 3))
  expect_true(grepl("A2~~0.857142857142857*A2", out_vector$conditions[[2]]$pop_synt, fixed = TRUE))
  expect_true(grepl("B3~~0.857142857142857*B3", out_vector$conditions[[2]]$pop_synt, fixed = TRUE))

  reliability_matrix <- matrix(c(0.8, 0.7, 0.9, 0.85), nrow = 2, byrow = TRUE)
  out_matrix <- suppressWarnings(
    powRICLPM(
      target_power = 0.8,
      sample_size = 1000,
      time_points = 4,
      ICC = 0.5,
      RI_cor = 0.3,
      lagged_effects = lagged_effects,
      within_cor = 0.3,
      reliability = reliability_matrix,
      estimate_ME = TRUE,
      reps = 1,
      seed = 123456
    )
  )

  expect_equal(length(out_matrix$conditions), 2)
  expect_equal(
    unname(vapply(out_matrix$conditions, function(x) x$reliability, character(1))),
    c("A=0.8, B=0.9", "A=0.7, B=0.85")
  )
  expect_equal(out_matrix$conditions[[2]]$reliability_matrix, matrix(c(0.7, 0.85), nrow = 2, ncol = 4))
  expect_true(grepl("A2~~0.857142857142857*A2", out_matrix$conditions[[2]]$pop_synt, fixed = TRUE))
  expect_true(grepl("B2~~0.352941176470588*B2", out_matrix$conditions[[2]]$pop_synt, fixed = TRUE))
  expect_true(grepl("A2~~start(0.857142857142857)*A2", out_matrix$conditions[[2]]$est_synt, fixed = TRUE))
  expect_true(grepl("B2~~start(0.352941176470588)*B2", out_matrix$conditions[[2]]$est_synt, fixed = TRUE))
})

test_that("bounded estimation for STARTS model in powRICLPM() works", {

  expect_warning({
    out1 <- powRICLPM(
      target_power = 0.8,
      sample_size = c(500),
      time_points = 4,
      ICC = .5,
      RI_cor = 0.3,
      lagged_effects = matrix(c(0.4, 0.15, 0.2, 0.3), ncol = 2, byrow = TRUE),
      within_cor = 0.3,
      reliability = .85,
      estimate_ME = TRUE,
      bounds = TRUE,
      reps = 1,
      seed = 1234
    )
  })

  expect_true(out1$session$bounds)
})

test_that("power analysis using Mplus works", {
  powRICLPM(
    sample_size = 1000,
    time_points = 4,
    ICC = 0.5,
    RI_cor = 0.3,
    lagged_effects = matrix(c(.5, .1, .4, .5), ncol = 2, byrow = TRUE),
    within_cor = 0.3,
    reps = 1000,
    seed = 123456,
    save_path = tempdir(),
    software = "Mplus"
  )

  if (.Platform$OS.type %in% c("windows", "mac")) {
    path <- normalizePath(
      file.path(tempdir(), "Condition1.inp"),
      mustWork = FALSE
    )

    expect_true(file.exists(path))
  }

})

test_that("power analysis for the STARTS model using Mplus works", {
  out_unconstrained <- powRICLPM(
    target_power = 0.8,
    sample_size = c(2000),
    time_points = 8,
    ICC = .5,
    RI_cor = 0.3,
    lagged_effects = matrix(c(0.4, 0.15, 0.2, 0.3), ncol = 2, byrow = TRUE),
    within_cor = 0.3,
    reliability = .85,
    estimate_ME = TRUE,
    reps = 2,
    seed = 1234,
    software = "Mplus",
    save_path = tempdir()
  )

  out_constrained <- powRICLPM(
    target_power = 0.8,
    sample_size = c(2000),
    time_points = 8,
    ICC = .5,
    RI_cor = 0.3,
    lagged_effects = matrix(c(0.4, 0.15, 0.2, 0.3), ncol = 2, byrow = TRUE),
    within_cor = 0.3,
    reliability = .85,
    estimate_ME = TRUE,
    reps = 2,
    seed = 1234,
    software = "Mplus",
    constraints = "ME",
    save_path = tempdir()
  )

  expect_error(
    powRICLPM(
      target_power = 0.8,
      sample_size = c(2000),
      time_points = 10,
      ICC = .5,
      RI_cor = 0.3,
      lagged_effects = matrix(c(0.4, 0.15, 0.2, 0.3), ncol = 2, byrow = TRUE),
      within_cor = 0.3,
      reliability = .85,
      estimate_ME = TRUE,
      reps = 2,
      seed = 1234,
      software = "Mplus",
      bounds = TRUE,
      save_path = tempdir()
    )
  )

  vector_reliability_save_path <- tempfile()
  dir.create(vector_reliability_save_path)
  powRICLPM(
    target_power = 0.8,
    sample_size = c(2000),
    time_points = 4,
    ICC = .5,
    RI_cor = 0.3,
    lagged_effects = matrix(c(0.4, 0.15, 0.2, 0.3), ncol = 2, byrow = TRUE),
    within_cor = 0.3,
    reliability = c(.85, .8),
    estimate_ME = TRUE,
    reps = 2,
    seed = 1234,
    software = "Mplus",
    save_path = vector_reliability_save_path
  )
  expect_true(file.exists(file.path(vector_reliability_save_path, "Condition1.inp")))
  expect_true(file.exists(file.path(vector_reliability_save_path, "Condition2.inp")))

})
