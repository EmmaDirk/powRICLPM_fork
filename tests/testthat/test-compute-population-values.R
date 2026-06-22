test_that("compute_Psi() works", {
  # Set lagged effects
  lagged_effects <- matrix(c(.2, .15, .10, .3), ncol = 2, byrow = TRUE)

  # Compute residual (co)variances
  output <- compute_Psi(lagged_effects = lagged_effects, within_cor = 0.3)

  # Run tests
  expect_type(output, "double")
  expect_equal(dim(output), c(2, 2))
  expect_equal(eigen(output)$values > 0, c(TRUE, TRUE))
})

test_that("compute_RI_var() works", {
  expect_equal(compute_RI_var(0), 0)
  expect_equal(compute_RI_var(.99), 99)
})

test_that("compute_RI_cov() works", {
  expect_equal(compute_RI_cov(0.3, compute_RI_var(0.5)), 0.3)
  expect_equal(compute_RI_cov(0.5, compute_RI_var(0.20)), 0.125)
})

test_that("compute_DPM_values() returns stationary DPM residual values", {
  lagged_effects <- matrix(c(.3, .1, .2, .25), ncol = 2, byrow = TRUE)
  loadings <- matrix(1, nrow = 2, ncol = 2)

  output <- compute_DPM_values(
    lagged_effects = lagged_effects,
    wave_cor = .2,
    AF_var = compute_AF_var(.2),
    AF_cov = compute_AF_cov(.3, compute_AF_var(.2)),
    loadings = loadings
  )

  expect_equal(dim(output$Psi), c(2, 2, 2))
  expect_equal(output$sigma_af, matrix(c(.2, .06, .06, .2), nrow = 2, byrow = TRUE))
  expect_true(all(eigen(output$Psi[, , 1])$values > 0))
  expect_equal(output$Psi[, , 1], output$Psi[, , 2])
})
