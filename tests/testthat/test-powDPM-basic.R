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
      dynamics_cor = 0.3,
      reps = 1
    ),
    "not identified with 3 waves.*requires at least 4 waves"
  )
  expect_error(
    powDPM(
      target_power = 0.8,
      sample_size = 1000,
      time_points = 4,
      AF_proportion = "0.2",
      AF_cor = 0.3,
      lagged_effects = lagged_effects,
      dynamics_cor = 0.3,
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
      dynamics_cor = 0.3,
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
      dynamics_cor = 0.3,
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
    "dynamics_cor.*missing"
  )
  expect_error(
    powDPM(
      target_power = 0.8,
      sample_size = 1000,
      time_points = 4,
      AF_proportion = 0.2,
      AF_cor = 0.3,
      lagged_effects = lagged_effects,
      dynamics_cor = 0.3,
      wave_cor = 0.3,
      reps = 1
    ),
    "unused argument"
  )
  expect_error(
    powDPM(
      target_power = 0.8,
      sample_size = 1000,
      time_points = 4,
      AF_proportion = 0.2,
      AF_cor = 0.3,
      lagged_effects = lagged_effects,
      dynamics_cor = 0.3,
      Phi = lagged_effects,
      reps = 1
    ),
    "unused argument"
  )
  expect_error(
    powDPM(
      target_power = 0.8,
      sample_size = 1000,
      time_points = 4,
      AF_proportion = 0.2,
      AF_cor = 0.3,
      lagged_effects = lagged_effects,
      dynamics_cor = 0.3,
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
        dynamics_cor = 0.3,
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
      dynamics_cor = 0.3,
      reliability = 0.8,
      estimate_ME = TRUE,
      reps = 1
    ),
    "not identified with 4 waves.*requires at least 5 waves"
  )
  expect_error(
    powDPM(
      target_power = 0.8,
      sample_size = 1000,
      time_points = 4,
      AF_proportion = 0.2,
      AF_cor = 0.3,
      lagged_effects = lagged_effects,
      dynamics_cor = 0.3,
      reliability = 0.8,
      estimate_ME = FALSE,
      constraints = "none",
      reps = 1,
      seed = 123456
    ),
    "Measurement error was generated but not estimated"
  )
  dpm_psi_error <- tryCatch(
    powDPM(
      target_power = 0.8,
      sample_size = 1000,
      time_points = 4,
      AF_proportion = 0.4,
      AF_cor = 0.2,
      lagged_effects = matrix(c(0.4, 0.1, 0.2, 0.3), ncol = 2, byrow = TRUE),
      dynamics_cor = 0.1,
      reps = 1
    ),
    error = function(e) e
  )
  expect_s3_class(dpm_psi_error, "error")
  expect_match(conditionMessage(dpm_psi_error), "dynamics_cor.*AF_proportion.*AF_cor")
  expect_false(grepl("within_cor", conditionMessage(dpm_psi_error), fixed = TRUE))

  out <- suppressWarnings(
    powDPM(
      target_power = 0.8,
      sample_size = 1000,
      time_points = 4,
      AF_proportion = c(0.1, 0.2),
      AF_cor = 0.3,
      lagged_effects = lagged_effects,
      dynamics_cor = 0.3,
      reps = 1,
      seed = 123456
    )
  )
  expect_s3_class(out, "powDPM")
  expect_equal(out$session$model, "DPM")
  expect_equal(out$conditions[[1]]$AF_proportion, 0.1)
  expect_equal(out$conditions[[2]]$AF_proportion, 0.2)
  expect_false("ICC" %in% names(out$conditions[[1]]))
  expect_false("RI_cor" %in% names(out$conditions[[1]]))
  expect_false("within_cor" %in% names(out$conditions[[1]]))
  expect_length(out$conditions, 2)
  expect_equal(out$session$argument_names$AF_proportion, "AF_proportion")
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
      dynamics_cor = 0.1,
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
      dynamics_cor = 0.1,
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
    dynamics_cor = 0.3,
    reps = 1,
    seed = 123456,
    constraints = "AF_loadings_free"
  )

  expect_error(
    do.call(powDPM, c(base, list(loadings = c(1, 0.8)))),
    "requires 3 loadings"
  )
  expect_error(
    do.call(powDPM, c(base, list(loadings = c(0.8, 1, 1)))),
    "first DPM loading must be 1"
  )
  expect_error(
    do.call(powDPM, c(base, list(loadings = c(NaN, 0.8, 0.7)))),
    "finite values"
  )
  expect_error(
    do.call(powDPM, c(base, list(loadings = c(1, Inf, 0.7)))),
    "finite values"
  )
  expect_error(
    do.call(powDPM, c(base, list(loadings = matrix(c(0.8, 0.7, 1, 1.2, 0.9, 1), nrow = 2, byrow = TRUE)))),
    "first DPM loading must be 1"
  )
})

test_that("powDPM notes custom data-generating loadings at argument checking", {
  lagged_effects <- matrix(c(0.4, 0.15, 0.2, 0.3), ncol = 2, byrow = TRUE)
  base <- list(
    target_power = 0.8,
    sample_size = 1000,
    time_points = 4,
    AF_proportion = 0.2,
    AF_cor = 0.3,
    lagged_effects = lagged_effects,
    dynamics_cor = 0.3,
    reps = 1,
    seed = 123456,
    constraints = "AF_loadings_free"
  )

  custom_messages <- capture.output(
    suppressWarnings(do.call(powDPM, c(base, list(loadings = c(1, 0.8, 1.1))))),
    type = "message"
  )
  default_messages <- capture.output(
    suppressWarnings(do.call(powDPM, base)),
    type = "message"
  )

  expect_true(any(grepl(
    "specified AF_proportion can only be interpreted as the accumulating-factor proportion at wave 2",
    custom_messages,
    fixed = TRUE
  )))
  expect_false(any(grepl(
    "specified AF_proportion can only be interpreted as the accumulating-factor proportion at wave 2",
    default_messages,
    fixed = TRUE
  )))
})
