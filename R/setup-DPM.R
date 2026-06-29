#' Set Up \code{powDPM} Analysis
#'
#' \code{create_conditions_DPM()} restructures the arguments of \code{powDPM()}
#' in a list, such that it can be used by \code{run_condition_monteCarlo()}.
#'
#' @inheritParams powDPM
#'
#' @noRd
create_conditions_DPM <- function(
  target_power,
  sample_size,
  time_points,
  AF_proportion,
  AF_cor,
  lagged_effects,
  dynamics_cor,
  reliability,
  loadings = NULL,
  skewness,
  kurtosis,
  estimate_ME,
  significance_criterion,
  reps,
  seed,
  constraints,
  bounds,
  estimator,
  software = "lavaan"
) {
  model <- "DPM"
  constraints <- normalize_constraints_for_software(constraints, software)
  reliability_conditions <- ireliability_conditions(reliability)

  conditions <- expand.grid(
    sample_size = sample_size,
    time_points = time_points,
    AF_proportion = AF_proportion,
    AF_cor = AF_cor,
    dynamics_cor = dynamics_cor,
    skewness = skewness,
    kurtosis = kurtosis,
    significance_criterion = significance_criterion,
    estimate_ME = estimate_ME,
    model = model,
    bounds = bounds,
    estimator = estimator,
    stringsAsFactors = FALSE
  )
  conditions <- conditions[rep(seq_len(nrow(conditions)), each = nrow(reliability_conditions)), , drop = FALSE]
  conditions <- cbind(
    conditions,
    reliability_conditions[rep(seq_len(nrow(reliability_conditions)), times = nrow(conditions) / nrow(reliability_conditions)), , drop = FALSE]
  )
  rownames(conditions) <- NULL

  conditions$constraints <- I(replicate(nrow(conditions), constraints, simplify = FALSE))
  conditions$lagged_effects <- replicate(nrow(conditions), lagged_effects, simplify = FALSE)
  conditions$loadings <- I(lapply(conditions$time_points, function(time_points) {
    inormalize_loadings(loadings, time_points, model)
  }))
  conditions$reliability_matrix <- I(Map(function(reliability, time_points) {
    inormalize_reliability(reliability, time_points)
  }, conditions$reliability_spec, conditions$time_points))
  conditions$reliability_spec <- NULL

  conditions$condition_id <- 1:nrow(conditions)
  conditions$AF_var <- sapply(conditions$AF_proportion, compute_AF_var)
  conditions$AF_cov <- mapply(compute_AF_cov, conditions$AF_cor, conditions$AF_var)
  conditions$ME_var <- lapply(conditions$reliability_matrix, function(reliability_matrix) {
    compute_ME_var(0, reliability_matrix)
  })
  conditions$DPM_values <- Map(
    compute_DPM_values,
    lagged_effects = conditions$lagged_effects,
    dynamics_cor = conditions$dynamics_cor,
    AF_var = conditions$AF_var,
    AF_cov = conditions$AF_cov,
    loadings = conditions$loadings
  )
  conditions$misspecification <- Map(
    detect_DPM_misspecification_condition,
    reliability_matrix = conditions$reliability_matrix,
    loadings = conditions$loadings,
    constraints = conditions$constraints,
    MoreArgs = list(estimate_ME = estimate_ME)
  )

  conditions <- split(conditions, seq(nrow(conditions)))
  conditions <- lapply(conditions, function(condition) {
    condition <- as.list(condition)
    condition$constraints <- condition$constraints[[1]]
    condition$lagged_effects <- condition$lagged_effects[[1]]
    condition$loadings <- condition$loadings[[1]]
    condition$reliability_matrix <- condition$reliability_matrix[[1]]
    condition$ME_var <- condition$ME_var[[1]]
    condition$DPM_values <- condition$DPM_values[[1]]
    condition$misspecification <- condition$misspecification[[1]]
    condition
  })

  lapply(conditions, create_lavaan_DPM)
}
