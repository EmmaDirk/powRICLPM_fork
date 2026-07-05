test_that("compute_MCSE_MSE() returns a standard error", {
  thetas_hat <- matrix(
    c(
      1, 3, 5,
      2, 4, 8
    ),
    nrow = 2,
    byrow = TRUE
  )
  population_values <- c(0, 1)
  squared_errors <- (thetas_hat - population_values)^2
  MSE <- rowMeans(squared_errors)

  expect_equal(
    compute_MCSE_MSE(thetas_hat, population_values, MSE),
    apply(squared_errors, 1, stats::sd) / sqrt(ncol(thetas_hat))
  )
})

test_that("run_condition_monteCarlo() computes power MCSE from completed replications", {
  expect_true(any(grepl(
    "compute_MCSE_power(power, reps_completed)",
    deparse(body(run_condition_monteCarlo)),
    fixed = TRUE
  )))
})

test_that("run_condition_monteCarlo() handles failed fits without handler argument mismatch", {
  body_text <- paste(deparse(body(run_condition_monteCarlo)), collapse = "\n")

  expect_true(grepl("error = function(e)", body_text, fixed = TRUE))
  expect_false(grepl("error = function(e, p)", body_text, fixed = TRUE))
})

test_that("run_condition_monteCarlo() does not use an initial lavaan slot cache fit", {
  body_text <- paste(deparse(body(run_condition_monteCarlo)), collapse = "\n")

  expect_false(grepl("use_lavaan_slot_cache", body_text, fixed = TRUE))
  expect_false(grepl("slotOptions = lav_options", body_text, fixed = TRUE))
  expect_true(grepl("estimator = estimator", body_text, fixed = TRUE))
})
