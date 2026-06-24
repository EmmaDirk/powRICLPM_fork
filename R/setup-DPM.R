#' Set Up \code{powDPM} Analysis
#'
#' \code{create_conditions_DPM()} restructures the arguments of \code{powDPM()}
#' in a list, such that it can be used by \code{run_condition_monteCarlo()}.
#'
#' @inheritParams powDPM
#'
#' @noRd
create_conditions_DPM <- function(
  model = "DPM",
  target_power,
  sample_size,
  time_points,
  intraclass_correlation,
  RI_cor,
  lagged_effects,
  within_cor,
  Psi,
  reliability,
  loadings = NULL,
  skewness,
  kurtosis,
  estimate_ME,
  significance_criterion,
  reps,
  bootstrap_reps,
  seed,
  constraints,
  bounds,
  estimator,
  save_path,
  software
) {
  AF_proportion <- intraclass_correlation
  constraints <- normalize_constraints_for_software(constraints, software)

  conditions <- expand.grid(
    sample_size = sample_size,
    time_points = time_points,
    ICC = AF_proportion,
    RI_cor = RI_cor,
    within_cor = within_cor,
    reliability = 1,
    skewness = skewness,
    kurtosis = kurtosis,
    significance_criterion = significance_criterion,
    estimate_ME = FALSE,
    model = model,
    stringsAsFactors = FALSE
  )

  conditions$constraints <- I(replicate(nrow(conditions), constraints, simplify = FALSE))
  conditions$lagged_effects <- replicate(nrow(conditions), lagged_effects, simplify = FALSE)
  conditions$loadings <- I(lapply(conditions$time_points, function(time_points) {
    inormalize_loadings(loadings, time_points, model)
  }))
  conditions$reliability_matrix <- I(lapply(conditions$time_points, function(time_points) {
    inormalize_reliability(1, time_points)
  }))

  conditions$condition_id <- 1:nrow(conditions)
  conditions$AF_var <- sapply(conditions$ICC, compute_AF_var)
  conditions$AF_cov <- mapply(compute_AF_cov, conditions$RI_cor, conditions$AF_var)
  conditions$DPM_values <- Map(
    compute_DPM_values,
    lagged_effects = conditions$lagged_effects,
    wave_cor = conditions$within_cor,
    AF_var = conditions$AF_var,
    AF_cov = conditions$AF_cov,
    loadings = conditions$loadings
  )

  conditions <- split(conditions, seq(nrow(conditions)))
  conditions <- lapply(conditions, function(condition) {
    condition <- as.list(condition)
    condition$constraints <- condition$constraints[[1]]
    condition$loadings <- condition$loadings[[1]]
    condition$reliability_matrix <- condition$reliability_matrix[[1]]
    condition$DPM_values <- condition$DPM_values[[1]]
    condition
  })

  lapply(conditions, create_lavaan_DPM)
}
