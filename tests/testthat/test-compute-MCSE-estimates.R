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
