#' Compute Residual Variances of Lagged Within-Components
#'
#' @description
#' \code{compute_Psi()} computes the variance-covariance matrix of the within-unit residuals (within a random intercept cross-lagged panel model) from wave 2 and later from specified \code{lagged_effects} and specified correlations between the within-components \code{within_cor}.
#'
#' @inheritParams powRICLPM
#'
#' @return A variance-covariance matrix for within-unit residuals from wave 2 and later.
#'
#' @details
#' The function is based on Equation (3.26) in \href{https://doi.org/10.7551/mitpress/6444.001.0001}{Kim and Nelson (1999, p. 27)}.
#'
#' @noRd
compute_Psi <- function(lagged_effects, within_cor) {
  wSigma <- matrix(c(1, within_cor, within_cor, 1), ncol = 2, byrow = TRUE)
  Psi <- matrix(
    data = (diag(length(wSigma)) - lagged_effects %x% lagged_effects) %*% c(wSigma),
    nrow = nrow(lagged_effects)
  )
  return(Psi)
}

#' Compute Random Intercept Variance
#'
#' \code{compute_RI_var()} computes the variance of the random intercept based on the proportion of between-unit variance \code{ICC} and conditional on within-unit variances of 1.
#'
#' @inheritParams powRICLPM
#'
#' @return A scalar representing the random intercept variance.
#'
#' @noRd
compute_RI_var <- function(ICC) {
  1 / (1 - ICC) - 1
}

#' Compute Random Intercept Correlation
#'
#' @param RI_cov A numeric value denoting the covariance between the random intercepts.
#' @param RI_var A numeric value denoting the variance of the random intercepts.
#'
#' @return A scalar representing the correlation between the random intercepts.
#'
#' @noRd
compute_RI_cov <- function(RI_cor, RI_var) {
  RI_cor * RI_var
}

#' Compute Measurement Error Variance
#'
#' @inheritParams compute_RI_cov
#' @param reliability A numeric value or matrix between 0 and 1, denoting the reliability of the variables.
#'
#' @noRd
compute_ME_var <- function(RI_var, reliability) {
  ((1 - reliability) * (1 + RI_var)) / reliability
}

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
