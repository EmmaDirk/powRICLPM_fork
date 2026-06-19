test_that("icheck_target() works", {
  expect_null(icheck_target(0.8))
  expect_error(icheck_target(1.4))
  expect_error(icheck_target("a"))
})

test_that("icheck_T() works", {
  expect_null(icheck_T(c(3, 4), ME = FALSE), c(3, 4))
  expect_error(icheck_T(3.5, ME = FALSE))
  expect_error(icheck_T(Inf, ME = FALSE), "finite")
  expect_error(icheck_T(NA_real_, ME = FALSE), "finite")
  expect_error(icheck_T(c(2, 3), ME = FALSE))
  expect_warning(icheck_T(c(3:30), ME = FALSE))
  expect_error(icheck_T(c(3, 4), ME = TRUE))
})

test_that("icheck_ICC() works", {
  expect_null(icheck_ICC(c(0.5, 0.8)))
  expect_error(icheck_ICC(-0.5))
  expect_error(icheck_ICC(0))
  expect_error(icheck_ICC(1))
  expect_error(icheck_ICC(2))
  expect_error(icheck_ICC("0.5"))
})

test_that("icheck_cor() works", {
  expect_null(icheck_cor(0.4))
  expect_error(icheck_cor(1.3))
  expect_error(icheck_cor("a"))
  expect_error(icheck_cor(c(0.3, 0.4)))
})

test_that("icheck_lagged_effects() works", {
  m1 <- matrix(c(.3, .2, .15, .2), ncol = 2, byrow = TRUE)
  m2 <- matrix(c(.8, .5, .4, .9), ncol = 2, byrow = TRUE)

  expect_null(icheck_lagged_effects(m1))
  expect_error(icheck_lagged_effects("m1"))
  expect_error(icheck_lagged_effects(m2))
})

test_that("check_lagged_effects() writes lagged effect interpretation", {
  m1 <- matrix(c(.3, .2, .15, .2), ncol = 2, byrow = TRUE)

  expect_output(check_lagged_effects(m1), "According to `lagged_effects`")
  expect_error(check_lagged_effects("m1"), "lagged_effects")
  expect_output(check_lagged_effects(Phi = m1), "According to `Phi`")
  expect_error(check_lagged_effects(Phi = 1), "Phi")
  expect_error(check_lagged_effects(m1, extra = TRUE), "Unexpected argument")
  expect_error(
    check_lagged_effects(lagged_effects = m1, Phi = m1),
    "Both.*lagged_effects.*Phi"
  )
})

test_that("check_Phi() writes Phi interpretation", {
  m1 <- matrix(c(.3, .2, .15, .2), ncol = 2, byrow = TRUE)

  expect_output(check_Phi(lagged_effects = m1), "According to `lagged_effects`")
  expect_output(check_Phi(Phi = m1), "According to `Phi`")
  expect_error(check_Phi(), "Phi")
  expect_error(check_Phi(lagged_effects = 1), "lagged_effects")
  expect_error(
    check_Phi(lagged_effects = m1, Phi = m1),
    "Both.*lagged_effects.*Phi"
  )
})

test_that("check_loadings() writes loading interpretation", {
  loading_vector <- c(1, 0.5, -1.2, 2)
  loading_matrix_same <- matrix(
    c(1, 0.5, -1.2, 2, 1, 0.5, -1.2, 2),
    nrow = 2,
    byrow = TRUE
  )
  loading_matrix <- matrix(
    c(1, 0, -1.2, 2, 1, 2.5, 0.25, -0.5),
    nrow = 2,
    byrow = TRUE
  )

  expect_output(
    check_loadings(loading_vector),
    "According to `loadings`"
  )
  expect_output(
    check_loadings(loading_vector),
    "RI_A loads on A3 and RI_B loads on B3 with -1.2"
  )
  expect_output(
    check_loadings(loading_vector),
    "RI_A loads on A4 and RI_B loads on B4 with 2"
  )
  expect_output(
    check_loadings(loading_matrix_same),
    "RI_A loads on A2 and RI_B loads on B2 with 0.5"
  )
  loading_matrix_same_output <- capture.output(check_loadings(loading_matrix_same))
  expect_equal(
    sum(grepl("^\\s*[*\u2022] RI_", loading_matrix_same_output)),
    4
  )

  expect_output(
    check_loadings(loading_matrix, time_points = 4),
    "RI_A loads on A2 with 0"
  )
  expect_output(
    check_loadings(loading_matrix, time_points = 4),
    "RI_B loads on B2 with 2.5"
  )
  expect_output(
    check_loadings(loading_matrix),
    "RI_B loads on B4 with -0.5"
  )
  loading_matrix_output <- capture.output(check_loadings(loading_matrix))
  expect_equal(
    sum(grepl("^\\s*[*\u2022] RI_", loading_matrix_output)),
    8
  )

  expect_error(check_loadings(), "must be supplied")
  expect_error(check_loadings(loading_vector, time_points = 3), "length 4.*time_points.*= 3")
  expect_error(
    check_loadings(c(0.8, 1, 1), time_points = 3),
    "first entry.*0.8.*vector beginning with 1.*matrix.*c\\(1, 1\\)"
  )
  expect_error(
    check_loadings(
      suppressWarnings(matrix(c(1, 0, 1, -0.5, 0.25), nrow = 2, byrow = TRUE)),
      time_points = 3
    ),
    "first column.*RI_A = 1.*RI_B = -0.5.*matrix.*c\\(1, 1\\)"
  )
  expect_error(check_loadings(list(c(1, 0, 1))), "numeric vector.*numeric matrix.*list")
  expect_error(check_loadings(data.frame(a = c(1, 0, 1))), "numeric vector.*numeric matrix.*data frame")
  expect_error(check_loadings(array(c(1, 0, -1, 2.5), dim = c(2, 2, 1))), "numeric vector.*numeric matrix.*array")
  expect_error(check_loadings(numeric(0)), "at least one value.*length 0")
  expect_error(check_loadings(loading_vector, time_points = c(3, 4)), "multiple values.*powRICLPM.*calls.*time_points.*= 2")
  expect_error(check_loadings(loading_vector, time_points = 3.5), "positive whole number.*time_points.*= 3.5")
  expect_error(check_loadings(loading_vector, time_points = "4"), 'positive whole number.*time_points.*= "4"')
  expect_error(check_loadings(loading_vector, time_points = NA), "positive whole number.*time_points.*= NA")
  expect_error(check_loadings(loading_vector, time_points = Inf), "positive whole number.*time_points.*= Inf")
  expect_error(check_loadings(loading_vector, time_points = 0), "positive whole number.*time_points.*= 0")
  expect_error(check_loadings(loading_vector, time_points = -4), "positive whole number.*time_points.*= -4")
  expect_error(check_loadings(loading_vector, time_points = 4, extra = TRUE), "Unexpected argument")
})

test_that("check_reliability() writes reliability interpretation", {
  reliability_vector <- c(0.8, 0.7, 1, 0.9)
  reliability_matrix_same <- matrix(
    c(0.8, 0.7, 1, 0.9, 0.8, 0.7, 1, 0.9),
    nrow = 2,
    byrow = TRUE
  )
  reliability_matrix <- matrix(
    c(0.8, 0.7, 1, 0.9, 0.9, 0.85, 0.75, 1),
    nrow = 2,
    byrow = TRUE
  )

  expect_output(check_reliability(0.8), "Variables A and B have reliability 0.8")
  expect_output(check_reliability(reliability_vector), "A2 and B2 have reliability 0.7")
  expect_output(check_reliability(reliability_matrix_same), "A3 and B3 have reliability 1")
  reliability_matrix_same_output <- capture.output(check_reliability(reliability_matrix_same))
  expect_equal(
    sum(grepl("^\\s*[*\u2022] A", reliability_matrix_same_output)),
    4
  )

  expect_output(check_reliability(reliability_matrix, time_points = 4), "A2 has reliability 0.7")
  expect_output(check_reliability(reliability_matrix, time_points = 4), "B2 has reliability 0.85")
  reliability_matrix_output <- capture.output(check_reliability(reliability_matrix))
  expect_equal(
    sum(grepl("^\\s*[*\u2022] [AB]", reliability_matrix_output)),
    8
  )

  expect_error(check_reliability(), "must be supplied")
  expect_error(check_reliability(reliability_vector, time_points = 3), "length 4.*time_points.*= 3")
  expect_error(check_reliability(list(c(0.8, 0.7, 1))), "numeric vector.*numeric matrix.*list")
  expect_error(check_reliability(data.frame(a = c(0.8, 0.7, 1))), "numeric vector.*numeric matrix.*data frame")
  expect_error(check_reliability(array(c(0.8, 0.7, 1, 0.9), dim = c(2, 2, 1))), "numeric vector.*numeric matrix.*array")
  expect_error(check_reliability(numeric(0)), "at least one value.*length 0")
  expect_error(check_reliability(c(0.8, 0.7, 1), time_points = c(3, 4)), "multiple values.*powRICLPM.*calls.*time_points.*= 2")
  expect_error(check_reliability(c(0.8, 0.7, 1), time_points = 3.5), "positive whole")
  expect_error(check_reliability(c(0.8, 0.7, 1), time_points = "3"), 'positive whole.*"3"')
  expect_error(check_reliability(c(0.8, 0.7, 1), time_points = NA), "positive whole.*NA")
  expect_error(check_reliability(c(0.8, 0.7, 1), time_points = Inf), "positive whole.*Inf")
  expect_error(check_reliability(c(0.8, 0.7, 1), time_points = 0), "positive whole.*0")
  expect_error(check_reliability(c(0.8, 0.7, 1), time_points = -3), "positive whole.*-3")
  expect_error(check_reliability(reliability_vector, time_points = 4, extra = TRUE), "Unexpected argument")
})

test_that("icheck_reliability() works", {
  expect_null(icheck_rel(.8))
  expect_null(icheck_rel(.8, c(3, 4), "Mplus"))
  expect_null(icheck_rel(c(.8, .9, 1), 3, "lavaan"))
  expect_null(icheck_rel(matrix(c(.8, .9, 1, .7, .8, .9), nrow = 2, byrow = TRUE), 3, "lavaan"))
  expect_error(icheck_rel(8), "larger than 0.1.*at most 1")
  expect_error(icheck_rel("a"), "numeric vector.*numeric matrix")
  expect_error(icheck_rel(data.frame(a = c(.8, .9, 1)), 3, "lavaan"), "data frame")
  expect_error(icheck_rel(array(c(.8, .9, 1, .7), dim = c(2, 2, 1)), 2, "lavaan"), "array")
  expect_error(icheck_rel(numeric(0)), "at least one value.*length 0")
  expect_error(icheck_rel(-.8), "larger than 0.1.*at most 1")
  expect_error(icheck_rel(.1), "larger than 0.1.*at most 1")
  expect_error(icheck_rel(c(.8, Inf, 1), 3, "lavaan"), "contains: c\\(0.8, Inf, 1\\)")
  expect_error(icheck_rel(c(.8, .9, 1), c(3, 4), "lavaan"), "compare power.*multiple values.*time_points")
  expect_error(icheck_rel(c(.8, .9, 1), 3, "Mplus"), "Time-varying reliability.*lavaan")
  expect_error(icheck_rel(matrix(.8, nrow = 1, ncol = 3), 3, "lavaan"), "must have 2 rows.*has 1 row")
  expect_error(icheck_rel(matrix(.8, nrow = 3, ncol = 3), 3, "lavaan"), "must have 2 rows.*has 3 rows")
  expect_error(icheck_rel(matrix(.8, nrow = 2, ncol = 4), 3, "lavaan"), "one column per time point.*has 4 columns.*time_points.*= 3")
  expect_error(icheck_rel(c(.8, .9), 3, "lavaan"), "length 2.*time_points.*= 3")
  expect_error(icheck_rel(c(.8, .9, 1), 3.5, "lavaan"), "positive whole")
  expect_error(icheck_rel(c(.8, .9, 1), "3", "lavaan"), "positive whole")
  expect_error(icheck_rel(c(.8, .9, 1), NA, "lavaan"), "positive whole")
  expect_error(icheck_rel(c(.8, .9, 1), Inf, "lavaan"), "positive whole")
  expect_error(icheck_rel(c(.8, .9, 1), 0, "lavaan"), "positive whole")
  expect_error(icheck_rel(c(.8, .9, 1), -3, "lavaan"), "positive whole")
})

test_that("icheck_loadings() works", {
  expect_null(icheck_loadings(NULL, 3, "lavaan"))
  expect_null(icheck_loadings(c(1, 1, 1), 3, "lavaan", "RI_loadings_free"))
  expect_null(icheck_loadings(c(1, 0, -2), 3, "lavaan", "RI_loadings_free"))
  expect_null(icheck_loadings(
    matrix(c(1, 0.5, -1, 1, 2, 0), nrow = 2, byrow = TRUE),
    3,
    "lavaan",
    "RI_loadings_free"
  ))

  expect_error(
    icheck_loadings(c(1, 2, 3), c(3, 4), "lavaan"),
    "compare power.*multiple values.*time_points.*multiple.*powRICLPM.*calls.*time_points.*= 2"
  )
  expect_error(
    icheck_loadings(c(1, 2, 3), 3, "Mplus"),
    "Time-varying random-intercept loadings"
  )
  expect_error(icheck_loadings("loadings", 3, "lavaan"), "numeric vector.*numeric matrix")
  expect_error(icheck_loadings(data.frame(a = c(1, 2, 3)), 3, "lavaan"), "data frame")
  expect_error(icheck_loadings(array(c(1, 2, 3, 4), dim = c(2, 2, 1)), 4, "lavaan"), "array")
  expect_error(icheck_loadings(c(1, 2, 3), 3.5, "lavaan"), "positive whole")
  expect_error(icheck_loadings(c(1, 2, 3), "3", "lavaan"), "positive whole")
  expect_error(icheck_loadings(c(1, 2, 3), NA, "lavaan"), "positive whole")
  expect_error(icheck_loadings(c(1, 2, 3), Inf, "lavaan"), "positive whole")
  expect_error(icheck_loadings(c(1, 2, 3), 0, "lavaan"), "positive whole")
  expect_error(icheck_loadings(c(1, 2, 3), -3, "lavaan"), "positive whole")
  expect_error(
    icheck_loadings(c(1, 2), 3, "lavaan"),
    "length 2.*time_points.*= 3"
  )
  expect_error(
    icheck_loadings(matrix(1, nrow = 3, ncol = 3), 3, "lavaan"),
    "must have 2 rows.*has 3 rows"
  )
  expect_error(
    icheck_loadings(matrix(1, nrow = 2, ncol = 4), 3, "lavaan"),
    "one column per time point.*has 4 columns.*time_points.*= 3"
  )
  expect_error(
    icheck_loadings(c(1, Inf, 2), 3, "lavaan"),
    "contains: c\\(1, Inf, 2\\)"
  )
  expect_error(
    icheck_loadings(c(0.8, 1, 1), 3, "lavaan"),
    "first entry.*0.8.*vector beginning with 1.*matrix.*c\\(1, 1\\)"
  )
  expect_error(
    icheck_loadings(matrix(c(1, 0.5, 1, 0.8, 1, 1), nrow = 2, byrow = TRUE), 3, "lavaan"),
    "first column.*RI_A = 1.*RI_B = 0.8.*matrix.*c\\(1, 1\\)"
  )
  expect_error(
    icheck_loadings(c(1, 0, -2), 3, "lavaan", "none"),
    "RI_loadings_free"
  )
  expect_error(
    icheck_loadings(c(1, 1, 1), 3, "lavaan", "lagged"),
    "in a constraint vector"
  )
  expect_error(
    icheck_loadings(c(1, 1, 1), 3, "lavaan", "lagged"),
    "supplied.*loadings.*constraints.*= lagged"
  )
})

test_that("icheck_moment() works", {
  expect_null(icheck_moment(0.3))
  expect_error(icheck_moment("a"))
  expect_error(icheck_moment(c(0.2, 0.5)))
})

test_that("icheck_significance_criterion() works", {
  expect_null(icheck_significance_criterion(0.05))
  expect_error(icheck_significance_criterion(c(0.05, 0.10)))
  expect_error(icheck_significance_criterion(-0.05))
  expect_error(icheck_significance_criterion("a"))
})

test_that("icheck_ME() works", {
  expect_null(icheck_ME(TRUE))
  expect_error(icheck_ME(c(T, F)))
  expect_error(icheck_ME("T"))
  expect_error(icheck_ME(1))
})

test_that("icheck_reps() works", {
  expect_null(icheck_reps(1000))
  expect_error(icheck_reps("1000"))
  expect_error(icheck_reps(1000.5))
  expect_error(icheck_reps(0), "positive")
  expect_error(icheck_reps(-1000))
})

test_that("icheck_seed() works", {
  expect_equal(icheck_seed(1234), 1234)
  expect_warning(icheck_seed(NA))
  expect_error(icheck_seed("1234"))
  expect_error(icheck_seed(1234.5))
})

test_that("icheck_constraints() works", {
  expect_null(icheck_constraints("lagged", ME = F))
  expect_null(icheck_constraints("ME", ME = T))
  expect_error(icheck_constraints(1, ME = F))
  expect_error(icheck_constraints("a", ME = F))
  expect_error(icheck_constraints(TRUE, ME = F))
  expect_error(icheck_constraints(c("none", "ME"), ME = F))
  expect_error(icheck_constraints("ME", ME = F))
})

test_that("vector constraints validate and preserve within compatibility", {
  expect_null(icheck_constraints(c("lagged", "residuals"), ME = FALSE))
  expect_null(icheck_constraints("RI_loadings_free", ME = FALSE))
  expect_null(icheck_constraints(c("lagged", "RI_loadings_free"), ME = FALSE))
  expect_null(icheck_constraints(c("residuals", "RI_loadings_free"), ME = FALSE))
  expect_null(icheck_constraints(c("stationarity", "RI_loadings_free"), ME = FALSE))
  expect_null(icheck_constraints(c("ME", "RI_loadings_free"), ME = TRUE))
  expect_silent(icheck_constraints("within", ME = FALSE))
  expect_silent(icheck_constraints(c("within", "RI_loadings_free"), ME = FALSE))

  expect_error(
    icheck_constraints(c("within", "lagged"), ME = FALSE),
    "cannot combine 'within'"
  )
  expect_error(
    icheck_constraints(c("within", "residuals"), ME = FALSE),
    "cannot combine 'within'"
  )
  expect_error(
    icheck_constraints(c("stationarity", "residuals"), ME = FALSE),
    "cannot combine 'stationarity'"
  )
  expect_error(
    icheck_constraints(c("ME", "RI_loadings_free"), ME = FALSE),
    "estimate_ME = TRUE"
  )

  expect_error(
    icheck_constraints_software("RI_loadings_free", "Mplus"),
    "not available for Mplus"
  )
  expect_null(icheck_constraints_software(c("lagged", "residuals"), "Mplus"))
  expect_null(icheck_constraints_software(c("lagged", "ME"), "Mplus"))
  expect_null(icheck_constraints_software(c("residuals", "ME"), "Mplus"))
  expect_null(icheck_constraints_software(c("stationarity", "ME"), "Mplus"))
  expect_error(
    icheck_constraints_software(c("lagged", "RI_loadings_free"), "Mplus"),
    "c\\('lagged', 'RI_loadings_free'\\)"
  )
})

test_that("constraint helpers interpret legacy within and explicit vectors", {
  expect_true(has_constraint("within", "lagged"))
  expect_true(has_constraint("within", "residuals"))
  expect_false(has_constraint("within", "stationarity"))
  expect_true(has_constraint(c("lagged", "residuals"), "lagged"))
  expect_true(has_constraint(c("lagged", "residuals"), "residuals"))
  expect_true(has_constraint(c("stationarity", "RI_loadings_free"), "lagged"))
  expect_true(has_constraint(c("stationarity", "RI_loadings_free"), "RI_loadings_free"))
  expect_false(has_constraint(c("lagged", "RI_loadings_free"), "residuals"))

  expect_equal(format_constraints("within"), "within")
  expect_equal(format_constraints(c("lagged", "residuals")), "c('lagged', 'residuals')")
  expect_equal(
    format_constraints(c("lagged", "RI_loadings_free")),
    "c('lagged', 'RI_loadings_free')"
  )
})

test_that("icheck_estimator() works", {
  expect_equal(icheck_estimator(NA, skewness = 0, kurtosis = 1), "MLR")
  expect_equal(icheck_estimator(NA, 0, 0), "ML")
  expect_equal(icheck_estimator(NA, skewness = 1, kurtosis = 0), "MLR")
  expect_equal(icheck_estimator("ML", skewness = 1, kurtosis = 1), "ML")
  expect_error(icheck_estimator("a", skewness = 0, kurtosis = 0))
  expect_error(icheck_estimator(1, skewness = 0, kurtosis = 0))
  expect_error(icheck_estimator(c("ML", "MLR"), skewness = 0, kurtosis = 0))
})

test_that("icheck_path() works", {
  expect_equal(icheck_path(tempdir(), "lavaan"), tempdir())
  expect_equal(icheck_path(tempdir(), "Mplus"), tempdir())
  expect_error(icheck_path(1, "lavaan"))
  expect_error(icheck_path("non_existing_path"))
  expect_equal(icheck_path(NULL, "Mplus"), getwd())
})

test_that("icheck_N() works", {
  expect_null(icheck_N(c(200, 300), 3, constraints = "none", ME = FALSE))
  expect_error(icheck_N(c(200.4, 300), 3, constraints = "none", ME = FALSE))
  expect_error(icheck_N(c(-200, 300), 3, constraints = "none", ME = FALSE))
  expect_error(icheck_N(10, 3, constraints = "none", ME = FALSE))
  expect_null(icheck_N(17, 3, constraints = "lagged", ME = FALSE))
  expect_error(icheck_N(20, 3, constraints = "none", ME = TRUE))
  expect_null(icheck_N(20, 3, constraints = "within", ME = TRUE))
})

test_that("icheck_sample_size_search() works", {
  expect_null(icheck_sample_size_search(100, 200, 20))
  expect_error(icheck_sample_size_search(NULL, 200, 20), "sample_size")
  expect_error(icheck_sample_size_search(100, NULL, 20), "search_upper")
  expect_error(icheck_sample_size_search(100, 200, NULL), "search_step")
  expect_error(icheck_sample_size_search("100", 200, 20), "single numeric")
  expect_error(icheck_sample_size_search(100.5, 200, 20), "integer")
  expect_error(icheck_sample_size_search(100, 200, 0), "positive")
  expect_error(icheck_sample_size_search(200, 100, 20), "search_upper")
})

test_that("icheck_bounds() works", {
  expect_null(icheck_bounds(TRUE, "none", "lavaan"))
  expect_null(icheck_bounds(TRUE, "RI_loadings_free", "lavaan"))
  expect_error(icheck_bounds("TRUE", "none", "Mplus"))
  expect_error(icheck_bounds(TRUE, "lagged", "lavaan"))
  expect_error(icheck_bounds(TRUE, c("lagged", "RI_loadings_free"), "lavaan"))
})

test_that("icheck_software() works", {
  expect_null(icheck_software("lavaan", 0, 0))
  expect_null(icheck_software("Mplus", 0, 0))
  expect_error(icheck_software(c("lavaan", "Mplus"), 0, 0))
  expect_error(icheck_software(TRUE))
  expect_error(icheck_software("Mplus", 1, 1))
})

test_that("vector constraints are one condition and match within shorthand", {
  lagged_effects <- matrix(c(0.4, 0.15, 0.2, 0.3), ncol = 2, byrow = TRUE)
  Psi <- compute_Psi(lagged_effects, within_cor = 0.3)

  conditions_within <- create_conditions(
    target_power = 0.8,
    sample_size = 1000,
    time_points = 3,
    intraclass_correlation = 0.5,
    RI_cor = 0.3,
    lagged_effects = lagged_effects,
    within_cor = 0.3,
    Psi = Psi,
    reliability = 1,
    skewness = 0,
    kurtosis = 0,
    estimate_ME = FALSE,
    significance_criterion = 0.05,
    reps = 1,
    bootstrap_reps = NULL,
    seed = 123456,
    constraints = "within",
    bounds = FALSE,
    estimator = "ML",
    save_path = NULL,
    software = "lavaan"
  )

  conditions_vector <- create_conditions(
    target_power = 0.8,
    sample_size = 1000,
    time_points = 3,
    intraclass_correlation = 0.5,
    RI_cor = 0.3,
    lagged_effects = lagged_effects,
    within_cor = 0.3,
    Psi = Psi,
    reliability = 1,
    skewness = 0,
    kurtosis = 0,
    estimate_ME = FALSE,
    significance_criterion = 0.05,
    reps = 1,
    bootstrap_reps = NULL,
    seed = 123456,
    constraints = c("lagged", "residuals"),
    bounds = FALSE,
    estimator = "ML",
    save_path = NULL,
    software = "lavaan"
  )

  expect_length(conditions_vector, 1)
  expect_equal(conditions_vector[[1]]$est_synt, conditions_within[[1]]$est_synt)

  conditions_mplus_within <- create_conditions(
    target_power = 0.8,
    sample_size = 1000,
    time_points = 3,
    intraclass_correlation = 0.5,
    RI_cor = 0.3,
    lagged_effects = lagged_effects,
    within_cor = 0.3,
    Psi = Psi,
    reliability = 1,
    skewness = 0,
    kurtosis = 0,
    estimate_ME = FALSE,
    significance_criterion = 0.05,
    reps = 1,
    bootstrap_reps = NULL,
    seed = 123456,
    constraints = "within",
    bounds = FALSE,
    estimator = "ML",
    save_path = tempdir(),
    software = "Mplus"
  )

  conditions_mplus_vector <- create_conditions(
    target_power = 0.8,
    sample_size = 1000,
    time_points = 3,
    intraclass_correlation = 0.5,
    RI_cor = 0.3,
    lagged_effects = lagged_effects,
    within_cor = 0.3,
    Psi = Psi,
    reliability = 1,
    skewness = 0,
    kurtosis = 0,
    estimate_ME = FALSE,
    significance_criterion = 0.05,
    reps = 1,
    bootstrap_reps = NULL,
    seed = 123456,
    constraints = c("lagged", "residuals"),
    bounds = FALSE,
    estimator = "ML",
    save_path = tempdir(),
    software = "Mplus"
  )

  expect_equal(
    conditions_mplus_vector[[1]]$Mplus_synt,
    conditions_mplus_within[[1]]$Mplus_synt
  )

  conditions_mplus_ME_lagged <- create_conditions(
    target_power = 0.8,
    sample_size = 1000,
    time_points = 4,
    intraclass_correlation = 0.5,
    RI_cor = 0.3,
    lagged_effects = lagged_effects,
    within_cor = 0.3,
    Psi = Psi,
    reliability = 0.8,
    skewness = 0,
    kurtosis = 0,
    estimate_ME = TRUE,
    significance_criterion = 0.05,
    reps = 1,
    bootstrap_reps = NULL,
    seed = 123456,
    constraints = c("ME", "lagged"),
    bounds = FALSE,
    estimator = "ML",
    save_path = tempdir(),
    software = "Mplus"
  )

  expect_true(grepl("(alpha)", conditions_mplus_ME_lagged[[1]]$Mplus_synt, fixed = TRUE))
  expect_true(grepl("(MEvarA)", conditions_mplus_ME_lagged[[1]]$Mplus_synt, fixed = TRUE))

  conditions_mplus_ME_stationarity <- create_conditions(
    target_power = 0.8,
    sample_size = 1000,
    time_points = 4,
    intraclass_correlation = 0.5,
    RI_cor = 0.3,
    lagged_effects = lagged_effects,
    within_cor = 0.3,
    Psi = Psi,
    reliability = 0.8,
    skewness = 0,
    kurtosis = 0,
    estimate_ME = TRUE,
    significance_criterion = 0.05,
    reps = 1,
    bootstrap_reps = NULL,
    seed = 123456,
    constraints = c("stationarity", "ME"),
    bounds = FALSE,
    estimator = "ML",
    save_path = tempdir(),
    software = "Mplus"
  )

  expect_true(grepl("MODEL CONSTRAINT", conditions_mplus_ME_stationarity[[1]]$Mplus_synt, fixed = TRUE))
  expect_true(grepl("(MEvarA)", conditions_mplus_ME_stationarity[[1]]$Mplus_synt, fixed = TRUE))
})

test_that("RI_loadings_free changes lavaan estimation syntax only", {
  lagged_effects <- matrix(c(0.4, 0.15, 0.2, 0.3), ncol = 2, byrow = TRUE)
  Psi <- compute_Psi(lagged_effects, within_cor = 0.3)

  conditions <- create_conditions(
    target_power = 0.8,
    sample_size = 1000,
    time_points = 3,
    intraclass_correlation = 0.5,
    RI_cor = 0.3,
    lagged_effects = lagged_effects,
    within_cor = 0.3,
    Psi = Psi,
    reliability = 1,
    skewness = 0,
    kurtosis = 0,
    estimate_ME = FALSE,
    significance_criterion = 0.05,
    reps = 1,
    bootstrap_reps = NULL,
    seed = 123456,
    constraints = "RI_loadings_free",
    bounds = FALSE,
    estimator = "ML",
    save_path = NULL,
    software = "lavaan"
  )

  condition <- conditions[[1]]

  expect_true(grepl("RI_A=~1*A2", condition$pop_synt, fixed = TRUE))
  expect_true(grepl("RI_B=~1*B2", condition$pop_synt, fixed = TRUE))
  expect_true(grepl("RI_A=~lx2*start(1)*A2", condition$est_synt, fixed = TRUE))
  expect_true(grepl("RI_A=~lx3*start(1)*A3", condition$est_synt, fixed = TRUE))
  expect_true(grepl("RI_B=~ly2*start(1)*B2", condition$est_synt, fixed = TRUE))
  expect_true(grepl("RI_B=~ly3*start(1)*B3", condition$est_synt, fixed = TRUE))
})

test_that("loadings update lavaan data generation syntax", {
  lagged_effects <- matrix(c(0.4, 0.15, 0.2, 0.3), ncol = 2, byrow = TRUE)
  Psi <- compute_Psi(lagged_effects, within_cor = 0.3)

  conditions_vector <- create_conditions(
    target_power = 0.8,
    sample_size = 1000,
    time_points = 3,
    intraclass_correlation = 0.5,
    RI_cor = 0.3,
    lagged_effects = lagged_effects,
    within_cor = 0.3,
    Psi = Psi,
    reliability = 1,
    loadings = c(1, 0.5, -1),
    skewness = 0,
    kurtosis = 0,
    estimate_ME = FALSE,
    significance_criterion = 0.05,
    reps = 1,
    bootstrap_reps = NULL,
    seed = 123456,
    constraints = "none",
    bounds = FALSE,
    estimator = "ML",
    save_path = NULL,
    software = "lavaan"
  )

  condition_vector <- conditions_vector[[1]]
  expect_true(grepl("RI_A=~0.5*A2", condition_vector$pop_synt, fixed = TRUE))
  expect_true(grepl("RI_A=~-1*A3", condition_vector$pop_synt, fixed = TRUE))
  expect_true(grepl("RI_B=~0.5*B2", condition_vector$pop_synt, fixed = TRUE))
  expect_true(grepl("RI_B=~-1*B3", condition_vector$pop_synt, fixed = TRUE))
  expect_true(grepl("RI_A=~1*A2", condition_vector$est_synt, fixed = TRUE))
  expect_true(grepl("RI_B=~1*B2", condition_vector$est_synt, fixed = TRUE))

  conditions_matrix <- create_conditions(
    target_power = 0.8,
    sample_size = 1000,
    time_points = 3,
    intraclass_correlation = 0.5,
    RI_cor = 0.3,
    lagged_effects = lagged_effects,
    within_cor = 0.3,
    Psi = Psi,
    reliability = 1,
    loadings = matrix(c(1, 0, -1.2, 1, 2, 0.25), nrow = 2, byrow = TRUE),
    skewness = 0,
    kurtosis = 0,
    estimate_ME = FALSE,
    significance_criterion = 0.05,
    reps = 1,
    bootstrap_reps = NULL,
    seed = 123456,
    constraints = "RI_loadings_free",
    bounds = FALSE,
    estimator = "ML",
    save_path = NULL,
    software = "lavaan"
  )

  condition_matrix <- conditions_matrix[[1]]
  expect_true(grepl("RI_A=~0*A2", condition_matrix$pop_synt, fixed = TRUE))
  expect_true(grepl("RI_A=~-1.2*A3", condition_matrix$pop_synt, fixed = TRUE))
  expect_true(grepl("RI_B=~2*B2", condition_matrix$pop_synt, fixed = TRUE))
  expect_true(grepl("RI_B=~0.25*B3", condition_matrix$pop_synt, fixed = TRUE))
  expect_true(grepl("RI_A=~lx2*start(0)*A2", condition_matrix$est_synt, fixed = TRUE))
  expect_true(grepl("RI_A=~lx3*start(-1.2)*A3", condition_matrix$est_synt, fixed = TRUE))
  expect_true(grepl("RI_B=~ly2*start(2)*B2", condition_matrix$est_synt, fixed = TRUE))
  expect_true(grepl("RI_B=~ly3*start(0.25)*B3", condition_matrix$est_synt, fixed = TRUE))
  expect_true(grepl("RI_A~~start(1)*RI_A", condition_matrix$est_synt, fixed = TRUE))
  expect_true(grepl("RI_A~~start(0.3)*RI_B", condition_matrix$est_synt, fixed = TRUE))
})

test_that("reliability updates lavaan measurement-error syntax", {
  lagged_effects <- matrix(c(0.4, 0.15, 0.2, 0.3), ncol = 2, byrow = TRUE)
  Psi <- compute_Psi(lagged_effects, within_cor = 0.3)

  conditions_scalar <- create_conditions(
    target_power = 0.8,
    sample_size = 1000,
    time_points = 3,
    intraclass_correlation = 0.5,
    RI_cor = 0.3,
    lagged_effects = lagged_effects,
    within_cor = 0.3,
    Psi = Psi,
    reliability = 0.8,
    loadings = NULL,
    skewness = 0,
    kurtosis = 0,
    estimate_ME = FALSE,
    significance_criterion = 0.05,
    reps = 1,
    bootstrap_reps = NULL,
    seed = 123456,
    constraints = "none",
    bounds = FALSE,
    estimator = "ML",
    save_path = NULL,
    software = "lavaan"
  )
  condition_scalar <- conditions_scalar[[1]]
  expect_equal(condition_scalar$reliability, 0.8)
  expect_equal(condition_scalar$reliability_matrix, matrix(0.8, nrow = 2, ncol = 3))
  expect_equal(condition_scalar$ME_var, matrix(0.5, nrow = 2, ncol = 3))
  expect_true(grepl("A1~~0.5*A1", condition_scalar$pop_synt, fixed = TRUE))
  expect_true(grepl("B3~~0.5*B3", condition_scalar$pop_synt, fixed = TRUE))

  conditions_vector <- create_conditions(
    target_power = 0.8,
    sample_size = 1000,
    time_points = 3,
    intraclass_correlation = 0.5,
    RI_cor = 0.3,
    lagged_effects = lagged_effects,
    within_cor = 0.3,
    Psi = Psi,
    reliability = c(0.8, 0.7, 1),
    loadings = NULL,
    skewness = 0,
    kurtosis = 0,
    estimate_ME = TRUE,
    significance_criterion = 0.05,
    reps = 1,
    bootstrap_reps = NULL,
    seed = 123456,
    constraints = "none",
    bounds = FALSE,
    estimator = "ML",
    save_path = NULL,
    software = "lavaan"
  )
  condition_vector <- conditions_vector[[1]]
  expected_vector_ME <- matrix(c(0.5, 6 / 7, 0, 0.5, 6 / 7, 0), nrow = 2, byrow = TRUE)
  expect_equal(condition_vector$reliability, "structured")
  expect_equal(condition_vector$reliability_matrix, matrix(c(0.8, 0.7, 1, 0.8, 0.7, 1), nrow = 2, byrow = TRUE))
  expect_equal(condition_vector$ME_var, expected_vector_ME)
  expect_true(grepl("A2~~0.857142857142857*A2", condition_vector$pop_synt, fixed = TRUE))
  expect_true(grepl("B3~~0*B3", condition_vector$pop_synt, fixed = TRUE))
  expect_true(grepl("A2~~start(0.857142857142857)*A2", condition_vector$est_synt, fixed = TRUE))
  expect_true(grepl("B3~~start(0)*B3", condition_vector$est_synt, fixed = TRUE))

  reliability_matrix <- matrix(c(0.8, 0.7, 1, 0.9, 0.85, 0.75), nrow = 2, byrow = TRUE)
  conditions_matrix <- create_conditions(
    target_power = 0.8,
    sample_size = 1000,
    time_points = 3,
    intraclass_correlation = 0.5,
    RI_cor = 0.3,
    lagged_effects = lagged_effects,
    within_cor = 0.3,
    Psi = Psi,
    reliability = reliability_matrix,
    loadings = NULL,
    skewness = 0,
    kurtosis = 0,
    estimate_ME = TRUE,
    significance_criterion = 0.05,
    reps = 1,
    bootstrap_reps = NULL,
    seed = 123456,
    constraints = "ME",
    bounds = FALSE,
    estimator = "ML",
    save_path = NULL,
    software = "lavaan"
  )
  condition_matrix <- conditions_matrix[[1]]
  expected_matrix_ME <- ((1 - reliability_matrix) * 2) / reliability_matrix
  expected_A_start <- mean(expected_matrix_ME[1, ])
  expected_B_start <- mean(expected_matrix_ME[2, ])
  expect_equal(condition_matrix$reliability, "structured")
  expect_equal(condition_matrix$reliability_matrix, reliability_matrix)
  expect_equal(condition_matrix$ME_var, expected_matrix_ME)
  expect_true(grepl("A2~~0.857142857142857*A2", condition_matrix$pop_synt, fixed = TRUE))
  expect_true(grepl("B2~~0.352941176470588*B2", condition_matrix$pop_synt, fixed = TRUE))
  expect_true(grepl(paste0("A1~~MEvarA*start(", expected_A_start, ")*A1"), condition_matrix$est_synt, fixed = TRUE))
  expect_true(grepl(paste0("B1~~MEvarB*start(", expected_B_start, ")*B1"), condition_matrix$est_synt, fixed = TRUE))
})

test_that("RI_loadings_free combines with compatible lavaan constraints", {
  lagged_effects <- matrix(c(0.4, 0.15, 0.2, 0.3), ncol = 2, byrow = TRUE)
  Psi <- compute_Psi(lagged_effects, within_cor = 0.3)

  compatible_constraints <- list(
    c("lagged", "RI_loadings_free"),
    c("residuals", "RI_loadings_free"),
    c("stationarity", "RI_loadings_free"),
    c("within", "RI_loadings_free"),
    c("ME", "RI_loadings_free")
  )

  for (constraints in compatible_constraints) {
    conditions <- create_conditions(
      target_power = 0.8,
      sample_size = 1000,
      time_points = 3,
      intraclass_correlation = 0.5,
      RI_cor = 0.3,
      lagged_effects = lagged_effects,
      within_cor = 0.3,
      Psi = Psi,
      reliability = 1,
      skewness = 0,
      kurtosis = 0,
      estimate_ME = TRUE,
      significance_criterion = 0.05,
      reps = 1,
      bootstrap_reps = NULL,
      seed = 123456,
      constraints = constraints,
      bounds = FALSE,
      estimator = "ML",
      save_path = NULL,
      software = "lavaan"
    )

    expect_true(grepl("RI_A=~lx2*start(1)*A2", conditions[[1]]$est_synt, fixed = TRUE))
    expect_true(grepl("RI_B=~ly2*start(1)*B2", conditions[[1]]$est_synt, fixed = TRUE))
  }
})

test_that("RI_loadings_free updates parameter counting and powRICLPM output", {
  lagged_effects <- matrix(c(0.4, 0.15, 0.2, 0.3), ncol = 2, byrow = TRUE)

  expect_equal(
    count_parameters(2, 3, "RI_loadings_free", FALSE),
    count_parameters(2, 3, "none", FALSE) + 4
  )

  out <- suppressWarnings(
    powRICLPM(
      target_power = 0.8,
      sample_size = 1000,
      time_points = 3,
      ICC = 0.5,
      RI_cor = 0.3,
      lagged_effects = lagged_effects,
      within_cor = 0.3,
      reps = 1,
      seed = 123456,
      constraints = "RI_loadings_free",
      bounds = TRUE
    )
  )

  expect_true(out$session$bounds)
  expect_true(all(
    c("RI_A=~A2", "RI_A=~A3", "RI_B=~B2", "RI_B=~B3") %in%
      out$conditions[[1]]$estimates$parameter
  ))
})

