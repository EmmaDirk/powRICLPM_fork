#' Power analysis for the DPM
#'
#' @description
#' Perform a Monte Carlo power analysis for the dynamic panel model (DPM). This
#' function uses the same lavaan simulation and summary machinery as
#' \code{\link{powRICLPM}}, but exposes the DPM-native arguments directly:
#' \code{AF_proportion}, \code{AF_cor}, and \code{dynamics_cor}.
#'
#' @param target_power A \code{double}, denoting the desired power.
#' @param search_lower,search_upper,search_step Optional sample-size search range.
#' @param sample_size An \code{integer} (vector), indicating sample size.
#' @param time_points An \code{integer} (vector), indicating number of time
#'   points. The DPM needs at least 3 waves, with the exact minimum depending
#'   on \code{estimate_ME}, \code{constraints}, and whether accumulating-factor
#'   loadings are freely estimated. If a specification is not identified with
#'   the supplied number of waves, use more waves or impose valid identifying
#'   constraints.
#' @param AF_proportion A \code{double} (vector) with elements between 0 and 1, denoting the proportion of DPM true-score variance attributed to the accumulating factors.
#' @param AF_cor A \code{double} between -1 and 1, denoting the correlation between accumulating factors in the DPM.
#' @param lagged_effects A 2 by 2 \code{matrix}, denoting the lagged effects.
#' @param dynamics_cor A \code{double} between -1 and 1, denoting the true-score dynamic-process correlation in the DPM.
#' @param reliability (optional) A \code{numeric} value, vector, or matrix with elements larger than 0.1 and at most 1, denoting the reliability of the observed variables in the DPM data-generating model.
#' @param estimate_ME (optional) A \code{logical}, denoting if measurement error variance should be estimated in the DPM.
#' @param loadings (optional) A \code{numeric} vector or matrix specifying accumulating-factor loadings for waves 2 through T in the DPM data-generating model. The wave-2 loading is fixed to 1.
#' @param skewness,kurtosis Numeric values used to generate nonnormal data.
#' @param significance_criterion Significance criterion used to compute power and coverage.
#' @param reps An integer number of Monte Carlo replications.
#' @param seed Optional random seed.
#' @param constraints Character vector of DPM constraints. Valid values are \code{"none"}, \code{"lagged"}, \code{"residuals"}, \code{"stationarity"}, \code{"ME"}, and \code{"AF_loadings_free"}.
#' @param bounds A \code{logical}, denoting if bounded lavaan estimation should be used.
#' @param estimator Lavaan estimator. Defaults to \code{"ML"} or \code{"MLR"} for nonnormal data.
#'
#' @details
#' The accumulating factors load on waves 2 through T and covary with the first
#' wave. The DPM uses \code{AF_proportion} as the accumulating-factor variance
#' and \code{dynamics_cor} as the target true-score process correlation. Residual
#' variances and covariances are computed from the DPM stationarity equations so
#' the DPM true-score variables have variance 1 across waves when the requested
#' values are admissible. Custom data-generating accumulating-factor loadings can
#' be supplied with \code{loadings}; the estimation model frees those loadings
#' only when \code{constraints = "AF_loadings_free"} is supplied. With custom
#' data-generating loadings, \code{AF_proportion} is not the accumulating-factor
#' proportion at every wave. When measurement error is included, observed variables
#' are single-indicator measurements of latent true-score variables named
#' \code{tA1}, \code{tA2}, ..., \code{tB1}, \code{tB2}, and so on.
#' The \code{reliability} argument controls the generated measurement error
#' variance. Setting \code{estimate_ME = TRUE} estimates those measurement error
#' variances in the DPM. As with the STARTS model, DPMs with estimated
#' measurement error may be prone to empirical under-identification,
#' non-convergence, or inadmissible solutions even when the model is formally
#' identified. Constraints such as \code{constraints = "ME"} or
#' \code{constraints = "stationarity"} can reduce model complexity.
#'
#' @return
#' An object with classes \code{powDPM} and \code{powRICLPM}, upon which
#' \code{summary()}, \code{print()}, and \code{plot()} can be used.
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{powRICLPM}}: Run power analysis for the RI-CLPM and STARTS model.
#'   \item \code{\link{summary.powRICLPM}}: Summarize the setup of a \code{powRICLPM} object.
#'   \item \code{\link{give}}: Extract information from \code{powRICLPM} objects.
#'   \item \code{\link{plot.powRICLPM}}: Visualize results from a \code{powRICLPM} object.
#' }
#'
#' @examples
#' lagged_effects <- matrix(c(.4, .1, .2, .3), ncol = 2, byrow = TRUE)
#'
#' \dontrun{
#' out <- powDPM(
#'   target_power = 0.8,
#'   sample_size = 500,
#'   time_points = 4,
#'   AF_proportion = 0.2,
#'   AF_cor = 0.3,
#'   lagged_effects = lagged_effects,
#'   dynamics_cor = 0.3,
#'   reps = 100,
#'   seed = 1234
#' )
#' }
#'
#' @export
powDPM <- function(
    target_power = 0.8,
    search_lower = NULL,
    search_upper = NULL,
    search_step = 20,
    sample_size = NULL,
    time_points,
    AF_proportion,
    AF_cor,
    lagged_effects,
    dynamics_cor,
    reliability = 1,
    estimate_ME = FALSE,
    loadings = NULL,
    skewness = 0,
    kurtosis = 0,
    significance_criterion = 0.05,
    reps = 20,
    seed = NA,
    constraints = "none",
    bounds = FALSE,
    estimator = NA
  ) {

  call_powDPM <- match.call()

  argument_names <- list(
    AF_proportion = "AF_proportion",
    AF_cor = "AF_cor",
    dynamics_cor = "dynamics_cor",
    lagged_effects = "lagged_effects"
  )

  time_start <- proc.time()

  cli::cli_h2("Checking Argument Input")

  reliability_expr <- call_powDPM$reliability
  if (!is.null(reliability_expr)) {
    icheck_reliability_matrix_call(
      expr = reliability_expr,
      env = parent.frame(),
      call = rlang::caller_env()
    )
  }

  icheck_target(target_power)
  icheck_T_DPM(time_points)
  icheck_AF_proportion(AF_proportion, arg = argument_names$AF_proportion)
  icheck_AF_cor(AF_cor, arg = argument_names$AF_cor)
  icheck_dynamics_cor(dynamics_cor, arg = argument_names$dynamics_cor)
  icheck_rel(reliability, time_points, "lavaan")
  icheck_ME(estimate_ME)
  icheck_DPM_lagged_effects(lagged_effects, arg = argument_names$lagged_effects)
  icheck_moment(skewness)
  icheck_moment(kurtosis)
  icheck_significance_criterion(significance_criterion)
  icheck_reps(reps)
  estimator <- icheck_estimator(estimator, skewness, kurtosis)
  icheck_DPM_constraints(constraints, estimate_ME = estimate_ME, time_points = time_points)
  icheck_loadings(loadings, time_points, "lavaan", constraints, model = "DPM")
  icheck_DPM_bounds(bounds)

  inote_ignored_input_names(sample_size, "sample_size", "Sample-size values")
  inote_ignored_input_names(time_points, "time_points", "Time-point values")
  inote_ignored_input_names(AF_proportion, argument_names$AF_proportion, "Accumulating-factor proportion values")
  inote_ignored_input_names(AF_cor, argument_names$AF_cor, "Correlation values")
  inote_ignored_input_names(dynamics_cor, argument_names$dynamics_cor, "Correlation values")
  inote_ignored_input_names(lagged_effects, argument_names$lagged_effects, "Lagged-effect values")
  inote_ignored_input_names(reliability, "reliability", "Reliability values")
  inote_ignored_input_names(loadings, "loadings", "Loading values")
  inote_ignored_input_names(constraints, "constraints", "Constraint values")

  if (is.null(sample_size)) {
    icheck_sample_size_search(search_lower, search_upper, search_step)
    sample_size <- seq(search_lower, search_upper, search_step)
  }
  icheck_DPM_identification(time_points, reliability, estimate_ME, constraints)
  icheck_N(sample_size, time_points, constraints, ME = estimate_ME, model = "DPM")
  if (length(seed) != 1L || !is.na(seed)) {
    seed <- icheck_seed(seed)
  }

  inote_custom_loadings_interpretation(loadings, time_points, model = "DPM")
  cli::cli_alert_success("Argument checking complete.")

  irun_power_analysis(
    model = "DPM",
    time_start = time_start,
    call = call_powDPM,
    argument_names = argument_names,
    target_power = target_power,
    sample_size = sample_size,
    time_points = time_points,
    intraclass_correlation = NULL,
    RI_cor = NULL,
    lagged_effects = lagged_effects,
    within_cor = NULL,
    AF_proportion = AF_proportion,
    AF_cor = AF_cor,
    dynamics_cor = dynamics_cor,
    Psi = NULL,
    reliability = reliability,
    loadings = loadings,
    skewness = skewness,
    kurtosis = kurtosis,
    estimate_ME = estimate_ME,
    significance_criterion = significance_criterion,
    reps = reps,
    bootstrap_reps = NULL,
    seed = seed,
    constraints = constraints,
    bounds = bounds,
    estimator = estimator,
    save_path = NULL,
    software = "lavaan"
  )
}
