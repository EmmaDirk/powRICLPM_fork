#' Count number of estimated parameters
#'
#' \code{count_parameters()} counts the number of parameters that is estimated in the condition with most time points, taking imposed constraints and the potential estimation of measurement errors into account.
#'
#' @inheritParams powRICLPM
#' @param k An integer denoting the number of variables in the RI-CLPM.
#'
#' @noRd
count_parameters <- function(k, time_points, constraints, est_ME, model = "RICLPM") {
  time_points_max <- max(time_points)

  if (model == "DPM") {
    n_parameters <- sum(
      factorial(1 + k) / (2 * (k - 1)), # Baseline observed (co)variances
      factorial(1 + k) / (2 * (k - 1)), # Accumulating factor (co)variances
      k^2, # Accumulating factor covariances with wave 1
      k^2 * (time_points_max - 1), # Lagged effects
      factorial(1 + k) / (2 * (k - 1)) * (time_points_max - 1) # Residual (co)variances
    )

    if (has_constraint(constraints, "lagged")) {
      n_parameters <- n_parameters - (k^2 * (time_points_max - 2))
    }
    if (has_constraint(constraints, "residuals")) {
      n_parameters <- n_parameters - ((k + sum((k - 1):1)) * (time_points_max - 2))
    }
    if (has_constraint(constraints, "stationarity")) {
      n_parameters <- n_parameters - (k + 1)
    }
    if (has_constraint(constraints, "loadings_free")) {
      n_parameters <- n_parameters + (k * max(time_points_max - 2, 0))
    }
    if (est_ME) {
      n_parameters <- n_parameters + k * time_points_max
    }
    if (est_ME && (has_constraint(constraints, "ME") ||
                   has_constraint(constraints, "stationarity"))) {
      n_parameters <- n_parameters - (k * (time_points_max - 1))
    }
    return(n_parameters)
  }

  # Default k-variate RI-CLPM
  n_parameters <- sum(
    factorial(1 + k) / (2 * (k - 1)) * (time_points_max + 1), # Number of (co)variances within- and between-level
    k^2 * (time_points_max - 1), # Number of lagged effects
    ifelse(est_ME, k * time_points_max, 0) # ME
  )

  # Influence of constraints
  if (has_constraint(constraints, "lagged")) {
    n_parameters <- n_parameters - (k^2 * (time_points_max - 2))
  }
  if (has_constraint(constraints, "residuals")) {
    n_parameters <- n_parameters - ((k + sum((k - 1):1)) * (time_points_max - 2))
  }
  if (has_constraint(constraints, "stationarity")) {
    n_parameters <- n_parameters - ((time_points_max - 1) * k)
  }
  if (est_ME && (has_constraint(constraints, "ME") ||
                 has_constraint(constraints, "stationarity"))) {
    n_parameters <- n_parameters - ((time_points_max - 1) * k)
  }
  if (has_constraint(constraints, "loadings_free")) {
    n_parameters <- n_parameters + (k * (time_points_max - 1))
  }
  return(n_parameters)
}


#' Count distinct variance-covariance elements
#'
#' @noRd
count_distinct_information <- function(k, time_points) {
  n_observed <- k * time_points
  n_observed * (n_observed + 1) / 2
}
