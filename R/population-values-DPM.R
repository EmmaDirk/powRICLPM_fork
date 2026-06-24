#' Compute Accumulating Factor Variance
#'
#' \code{compute_AF_var()} computes the accumulating factor variance for the
#' dynamic panel model. The value is interpreted as the proportion of observed
#' variance attributed to the accumulating factors.
#'
#' @noRd
compute_AF_var <- function(AF_proportion) {
  AF_proportion
}

#' Compute Accumulating Factor Covariance
#'
#' @noRd
compute_AF_cov <- function(AF_cor, AF_var) {
  AF_cor * AF_var
}

#' Compute Dynamic Panel Model Population Values
#'
#' @noRd
compute_DPM_values <- function(lagged_effects, wave_cor, AF_var, AF_cov, loadings) {
  sigma_s <- matrix(c(1, wave_cor, wave_cor, 1), nrow = 2, byrow = TRUE)
  sigma_af <- matrix(c(AF_var, AF_cov, AF_cov, AF_var), nrow = 2, byrow = TRUE)
  gamma <- sigma_af %*% solve(diag(2) - t(lagged_effects))

  n_lagged_waves <- ncol(loadings)
  psi <- array(NA_real_, dim = c(2, 2, n_lagged_waves))
  gamma_start <- gamma
  gamma_values <- array(NA_real_, dim = c(2, 2, n_lagged_waves + 1))
  gamma_values[, , 1] <- gamma_start

  for (i in seq_len(n_lagged_waves)) {
    lambda <- diag(loadings[, i], nrow = 2)
    psi[, , i] <- sigma_s -
      lambda %*% sigma_af %*% t(lambda) -
      lagged_effects %*% sigma_s %*% t(lagged_effects) -
      lambda %*% gamma %*% t(lagged_effects) -
      lagged_effects %*% t(gamma) %*% t(lambda)
    gamma <- sigma_af %*% t(lambda) + gamma %*% t(lagged_effects)
    gamma_values[, , i + 1] <- gamma
  }

  list(
    sigma_s = sigma_s,
    sigma_af = sigma_af,
    gamma_start = gamma_start,
    gamma = gamma_values,
    Psi = psi
  )
}
