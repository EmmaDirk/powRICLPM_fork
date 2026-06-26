#' Set Up \code{powRICLPM} Analysis
#'
#' \code{create_conditions()} restructures the arguments of \code{powRICLPM()} in a list, such that it can be used by \code{run_condition_MonteCarlo()}, performing a Monte Carlo power analysis for each experimental condition.
#'
#' @inheritParams powRICLPM
#'
#' @noRd
create_conditions <- function(
  model = "RICLPM",
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
  if (identical(model, "DPM")) {
    return(create_conditions_DPM(
      model = model,
      target_power = target_power,
      sample_size = sample_size,
      time_points = time_points,
      AF_proportion = intraclass_correlation,
      AF_cor = RI_cor,
      lagged_effects = lagged_effects,
      dynamics_cor = within_cor,
      Psi = Psi,
      reliability = reliability,
      loadings = loadings,
      skewness = skewness,
      kurtosis = kurtosis,
      estimate_ME = estimate_ME,
      significance_criterion = significance_criterion,
      reps = reps,
      bootstrap_reps = bootstrap_reps,
      seed = seed,
      constraints = constraints,
      bounds = bounds,
      estimator = estimator,
      save_path = save_path,
      software = software
    ))
  }

  ICC <- intraclass_correlation
  constraints <- normalize_constraints_for_software(constraints, software)

  reliability_conditions <- ireliability_conditions(reliability)

  # Create data.frame with rows as experimental conditions
  conditions <- expand.grid(
    sample_size = sample_size,
    time_points = time_points,
    ICC = ICC,
    RI_cor = RI_cor,
    within_cor = within_cor,
    skewness = skewness,
    kurtosis = kurtosis,
    significance_criterion = significance_criterion,
    estimate_ME = estimate_ME,
    model = model,
    stringsAsFactors = FALSE
  )
  conditions <- conditions[rep(seq_len(nrow(conditions)), each = nrow(reliability_conditions)), , drop = FALSE]
  conditions <- cbind(
    conditions,
    reliability_conditions[rep(seq_len(nrow(reliability_conditions)), times = nrow(conditions) / nrow(reliability_conditions)), , drop = FALSE]
  )
  rownames(conditions) <- NULL

  # Add matrix input to conditions
  conditions$constraints <- I(replicate(nrow(conditions), constraints, simplify = FALSE))
  conditions$lagged_effects <- replicate(nrow(conditions), lagged_effects, simplify = FALSE)
  conditions$Psi <- replicate(nrow(conditions), Psi, simplify = FALSE)
  conditions$loadings <- I(lapply(conditions$time_points, function(time_points) {
    inormalize_loadings(loadings, time_points, model)
  }))
  conditions$reliability_matrix <- I(Map(function(reliability, time_points) {
    inormalize_reliability(reliability, time_points)
  }, conditions$reliability_spec, conditions$time_points))
  conditions$reliability_spec <- NULL

  # Compute and add additional parameters per condition
  conditions$condition_id <- 1:nrow(conditions)
  conditions$RI_var <- sapply(conditions$ICC, compute_RI_var)
  conditions$RI_cov <- mapply(compute_RI_cov, conditions$RI_cor, conditions$RI_var)
  conditions$ME_var <- Map(compute_ME_var, conditions$RI_var, conditions$reliability_matrix)

  # Create list of conditions
  conditions <- split(conditions, seq(nrow(conditions)))
  conditions <- lapply(conditions, function(condition) {
    condition <- as.list(condition)
    condition$constraints <- condition$constraints[[1]]
    condition$loadings <- condition$loadings[[1]]
    condition$reliability_matrix <- condition$reliability_matrix[[1]]
    condition$ME_var <- condition$ME_var[[1]]
    condition$misspecification <- detect_RICLPM_restrictive_misspecification(
      reliability_matrix = condition$reliability_matrix,
      estimate_ME = condition$estimate_ME,
      loadings = condition$loadings,
      constraints = condition$constraints
    )
    condition
  })

  # Create syntax per condition
  if (software == "lavaan") {
    conditions <- lapply(conditions, create_lavaan)
  } else if (software == "Mplus") {
    conditions <- lapply(conditions, create_Mplus, reps = reps , seed = seed)
    lapply(conditions, save_condition_Mplus, save_path)
  }

  return(conditions)

}


