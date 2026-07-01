test_that("DPM time-point checks work", {
  expect_null(icheck_T_DPM(3))
  expect_null(icheck_T_DPM(4))
  expect_error(icheck_T_DPM("4"), "integers")
  expect_error(icheck_T_DPM(NA_real_), "finite")
  expect_error(icheck_T_DPM(2.5), "integers")
  expect_error(icheck_T_DPM(2), "at least 3 time points")
})

test_that("DPM constraint checks reject unsupported options", {
  expect_null(icheck_DPM_constraints("none"))
  expect_null(icheck_DPM_constraints(c("lagged", "residuals")))
  expect_null(icheck_DPM_constraints("AF_loadings_free"))
  expect_null(icheck_DPM_constraints("ME", estimate_ME = TRUE))
  expect_true(has_constraint(c("stationarity", "AF_loadings_free"), "loadings_free"))
  expect_error(icheck_DPM_constraints("ME"), "estimate_ME = TRUE")
  expect_error(icheck_DPM_constraints("within"), "invalid DPM constraints")
  expect_error(icheck_DPM_constraints("loadings_free"), "invalid DPM constraints")
  expect_error(icheck_DPM_constraints("RI_loadings_free"), "invalid DPM constraints")
  expect_error(icheck_DPM_constraints(c("stationarity", "lagged")), "cannot combine")
})

test_that("DPM lagged-effects stationarity error is concise", {
  lagged_effects <- matrix(c(0.9, 0.7, 0.7, 0.9), nrow = 2)
  stationarity_error <- tryCatch(
    icheck_DPM_lagged_effects(lagged_effects),
    error = conditionMessage
  )

  expect_match(stationarity_error, "must specify a stationary process", fixed = TRUE)
  expect_match(stationarity_error, "Use smaller autoregressive and/or cross-lagged effects", fixed = TRUE)
  expect_match(stationarity_error, "largest absolute eigenvalue", fixed = TRUE)
  expect_false(grepl("Try smaller lagged effects", stationarity_error, fixed = TRUE))
  expect_false(grepl("unit circle", stationarity_error, fixed = TRUE))
})

test_that("DPM identification table is enforced", {
  for (i in seq_len(nrow(DPM_identification_table))) {
    row <- DPM_identification_table[i, ]
    constraints <- strsplit(row$constraints_key, "+", fixed = TRUE)[[1]]
    constraints <- constraints[constraints != "none"]
    if (isTRUE(row$AF_loadings_free)) {
      constraints <- c(constraints, "AF_loadings_free")
    }
    if (length(constraints) == 0L) {
      constraints <- "none"
    }
    if (row$min_waves > 3L) {
      expect_error(
        icheck_DPM_identification(
          time_points = row$min_waves - 1L,
          reliability = if (row$estimate_ME) .8 else 1,
          estimate_ME = row$estimate_ME,
          constraints = constraints
        ),
        "not identified"
      )
    }
    expect_null(icheck_DPM_identification(
      time_points = row$min_waves,
      reliability = if (row$estimate_ME) .8 else 1,
      estimate_ME = row$estimate_ME,
      constraints = constraints
    ))
  }
})

test_that("DPM misspecification detection separates restrictive and general cases", {
  restrictive <- detect_DPM_misspecification(
    reliability = .8,
    loadings = NULL,
    time_points = 4,
    constraints = "none",
    estimate_ME = FALSE
  )
  expect_true(restrictive$restrictive)
  expect_match(restrictive$restrictive_reasons, "Measurement error was generated, but not estimated", fixed = TRUE)

  general <- detect_DPM_misspecification(
    reliability = 1,
    loadings = NULL,
    time_points = 4,
    constraints = "none",
    estimate_ME = TRUE
  )
  expect_false(general$restrictive)
  expect_true(general$general)
})

test_that("check_loadings() writes DPM loading interpretation", {
  loading_matrix <- matrix(
    c(1, 0.8, 1.1, 1, 1.2, 0.9),
    nrow = 2,
    byrow = TRUE
  )

  expect_output(
    check_loadings(c(1, 0.8), model = "DPM"),
    "AF_A loads on A3 and AF_B loads on B3 with 0.8"
  )
  expect_output(
    check_loadings(loading_matrix, model = "DPM"),
    "AF_B loads on B4 with 0.9"
  )
  expect_output(
    check_loadings(loading_matrix, time_points = 4, model = "DPM"),
    "AF_A loads on A4 with 1.1"
  )
  expect_error(
    check_loadings(c(0.8, 1), model = "DPM"),
    "first DPM loading must be 1"
  )
  expect_error(
    check_loadings(c(1, NA), model = "DPM"),
    "finite values"
  )
  expect_error(
    check_loadings(c(1, Inf), model = "DPM"),
    "finite values"
  )
  expect_error(
    check_loadings(c(1, 0.8), time_points = 4, model = "DPM"),
    "requires 3 loadings"
  )
  expect_error(
    check_loadings(c(1), time_points = 4, model = "DPM"),
    "requires 3 loadings"
  )
  expect_error(
    check_loadings(matrix(c(1, 0.8, 1, 0.8), nrow = 2, byrow = TRUE), time_points = 4, model = "DPM"),
    "requires 3 loading columns"
  )
  expect_error(
    check_loadings(loading_matrix, model = "dpm"),
    "RICLPM.*DPM.*dpm"
  )
})

test_that("icheck_loadings() supports DPM loading specifications", {
  expect_null(icheck_loadings(c(1, .8), 3, "lavaan", "AF_loadings_free", model = "DPM"))
  expect_null(icheck_loadings(
    matrix(c(1, .8, 1, 1.2), nrow = 2, byrow = TRUE),
    3,
    "lavaan",
    "AF_loadings_free",
    model = "DPM"
  ))

  expect_error(
    icheck_loadings(c(1, .8, .7), 3, "lavaan", "AF_loadings_free", model = "DPM"),
    "requires 2 loadings"
  )
  expect_error(
    icheck_loadings(c(.8, 1), 3, "lavaan", "AF_loadings_free", model = "DPM"),
    "first DPM loading must be 1"
  )
  expect_error(
    icheck_loadings(matrix(c(.8, 1, 1, 1.2), nrow = 2, byrow = TRUE), 3, "lavaan", "AF_loadings_free", model = "DPM"),
    "first DPM loading must be 1"
  )
  expect_error(
    icheck_loadings(matrix(c(1, .8, 1, 1.2), nrow = 1, byrow = TRUE), 3, "lavaan", "AF_loadings_free", model = "DPM"),
    "must have 2 rows"
  )
  expect_error(
    icheck_loadings(matrix(c(1, .8, .7, 1, 1.2, .9), nrow = 2, byrow = TRUE), 3, "lavaan", "AF_loadings_free", model = "DPM"),
    "requires 2 loading columns"
  )
  expect_error(
    icheck_loadings(c(1, NA), 3, "lavaan", "AF_loadings_free", model = "DPM"),
    "finite values"
  )
  expect_error(
    icheck_loadings(matrix(c(NaN, .8, 1, 1.2), nrow = 2, byrow = TRUE), 3, "lavaan", "AF_loadings_free", model = "DPM"),
    "finite values"
  )
  expect_error(
    icheck_loadings(c(1, .8), c(3, 4), "lavaan", "AF_loadings_free", model = "DPM"),
    "one value of.*time_points"
  )
})

test_that("DPM parameter-counting checks catch underidentified sample-size paths", {
  expect_error(
    icheck_N(1000, 4, constraints = "none", ME = TRUE, model = "DPM"),
    "not identified.*39 parameters.*36 distinct"
  )
  expect_error(
    icheck_N(1000, 4, constraints = "AF_loadings_free", ME = TRUE, model = "DPM"),
    "not identified.*43 parameters.*36 distinct"
  )
  expect_error(
    icheck_N(1000, 4, constraints = c("ME", "AF_loadings_free"), ME = TRUE, model = "DPM"),
    "not identified.*37 parameters.*36 distinct"
  )
  expect_null(icheck_N(1000, 4, constraints = "ME", ME = TRUE, model = "DPM"))
  expect_null(icheck_N(
    1000, 4, constraints = c("stationarity", "AF_loadings_free"),
    ME = TRUE, model = "DPM"
  ))
  expect_equal(
    count_parameters(2, 4, "stationarity", TRUE, model = "DPM"),
    count_parameters(2, 4, c("stationarity", "ME"), TRUE, model = "DPM")
  )
})

test_that("DPM lavaan syntax uses accumulating factors and observed lagged effects", {
  lagged_effects <- matrix(c(0.3, 0.1, 0.2, 0.25), ncol = 2, byrow = TRUE)

  conditions <- create_conditions_DPM(
    target_power = 0.8,
    sample_size = 1000,
    time_points = 2,
    AF_proportion = 0.2,
    AF_cor = 0.3,
    lagged_effects = lagged_effects,
    dynamics_cor = 0.2,
    reliability = 1,
    loadings = NULL,
    skewness = 0,
    kurtosis = 0,
    estimate_ME = FALSE,
    significance_criterion = 0.05,
    reps = 1,
    seed = 123456,
    constraints = "none",
    bounds = FALSE,
    estimator = "ML"
  )

  condition <- conditions[[1]]
  expect_equal(condition$model, "DPM")
  expect_equal(condition$AF_var, 0.2)
  expect_equal(condition$AF_cov, 0.06)
  expect_true(grepl("AF_A=~1*A2", condition$pop_synt, fixed = TRUE))
  expect_false(grepl("AF_A=~1*A1", condition$pop_synt, fixed = TRUE))
  expect_true(grepl("A2~0.3*A1", condition$pop_synt, fixed = TRUE))
  expect_false(grepl("wA", condition$pop_synt, fixed = TRUE))
  expect_true(grepl("AF_A~~gAA1*start", condition$est_synt, fixed = TRUE))
})

test_that("DPM measurement error syntax uses latent true-score process", {
  lagged_effects <- matrix(c(0.3, 0.1, 0.2, 0.25), ncol = 2, byrow = TRUE)

  condition <- create_conditions_DPM(
    target_power = 0.8,
    sample_size = 1000,
    time_points = 4,
    AF_proportion = 0.2,
    AF_cor = 0.3,
    lagged_effects = lagged_effects,
    dynamics_cor = 0.2,
    reliability = 0.8,
    loadings = NULL,
    skewness = 0,
    kurtosis = 0,
    estimate_ME = TRUE,
    significance_criterion = 0.05,
    reps = 1,
    seed = 123456,
    constraints = "ME",
    bounds = FALSE,
    estimator = "ML"
  )[[1]]

  expect_equal(condition$ME_var, matrix(0.25, nrow = 2, ncol = 4))
  expect_true(grepl("tA1=~1*A1", condition$pop_synt, fixed = TRUE))
  expect_true(grepl("tB4=~1*B4", condition$pop_synt, fixed = TRUE))
  expect_true(grepl("AF_A=~1*tA2", condition$pop_synt, fixed = TRUE))
  expect_false(grepl("AF_A=~1*A2", condition$pop_synt, fixed = TRUE))
  expect_true(grepl("tA2~0.3*tA1", condition$pop_synt, fixed = TRUE))
  expect_true(grepl("A1~~0.25*A1", condition$pop_synt, fixed = TRUE))
  expect_true(grepl("A1~~MEvarA*start(0.25)*A1", condition$est_synt, fixed = TRUE))
  expect_true(grepl("tB2~start(0.2)*tA1", condition$est_synt, fixed = TRUE))
})

test_that("DPM stationarity syntax follows fixed and free loading equations", {
  lagged_effects <- matrix(c(0.3, 0.1, 0.2, 0.25), ncol = 2, byrow = TRUE)

  fixed <- create_conditions_DPM(
    target_power = 0.8,
    sample_size = 1000,
    time_points = 3,
    AF_proportion = 0.2,
    AF_cor = 0.3,
    lagged_effects = lagged_effects,
    dynamics_cor = 0.2,
    reliability = 1,
    loadings = NULL,
    skewness = 0,
    kurtosis = 0,
    estimate_ME = FALSE,
    significance_criterion = 0.05,
    reps = 1,
    seed = 123456,
    constraints = "stationarity",
    bounds = FALSE,
    estimator = "ML"
  )[[1]]

  free <- create_conditions_DPM(
    target_power = 0.8,
    sample_size = 1000,
    time_points = 4,
    AF_proportion = 0.2,
    AF_cor = 0.3,
    lagged_effects = lagged_effects,
    dynamics_cor = 0.2,
    reliability = 1,
    loadings = c(1, 0.8, 1.1),
    skewness = 0,
    kurtosis = 0,
    estimate_ME = FALSE,
    significance_criterion = 0.05,
    reps = 1,
    seed = 123456,
    constraints = c("stationarity", "AF_loadings_free"),
    bounds = FALSE,
    estimator = "ML"
  )[[1]]

  expect_true(grepl("gAA1==vfa + gAA1*a + gAB1*b", fixed$est_synt, fixed = TRUE))
  expect_true(grepl("rvarA2==1 - vfa", fixed$est_synt, fixed = TRUE))
  expect_true(grepl("AF_A=~lx3*start(0.8)*A3", free$est_synt, fixed = TRUE))
  expect_true(grepl("gAA3:=vfa*lx3 + gAA2*a + gAB2*b", free$est_synt, fixed = TRUE))
  expect_true(grepl("rcov4==wc - lx4*ly4*cff", free$est_synt, fixed = TRUE))
})
