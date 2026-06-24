#' Power analysis for the DPM
#'
#' @description
#' Perform a Monte Carlo power analysis for the dynamic panel model (DPM). This
#' function uses the same lavaan simulation and summary machinery as
#' \code{\link{powRICLPM}}, but exposes the DPM-native arguments directly:
#' \code{AF_proportion}, \code{AF_cor}, and \code{wave_cor}.
#'
#' @inheritParams powRICLPM
#' @param AF_proportion A \code{double} (vector) with elements between 0 and 1, denoting the proportion of observed variance attributed to the accumulating factors in the DPM.
#' @param AF_cor A \code{double} between -1 and 1, denoting the correlation between accumulating factors in the DPM.
#' @param wave_cor A \code{double} between -1 and 1, denoting the observed wave-level correlation in the DPM.
#' @param loadings (optional) A \code{numeric} vector or matrix specifying accumulating-factor loadings for the DPM in the \pkg{lavaan} data-generating model. DPM loadings must include one value per wave, with an explicit first-wave \code{NA}. The wave-2 loading is fixed to 1.
#' @param Phi Alternative name for \code{lagged_effects}.
#'
#' @details
#' The accumulating factors load on waves 2 through T and covary with the first
#' wave. The DPM uses \code{AF_proportion} as the accumulating-factor variance
#' and \code{wave_cor} as the target observed wave-level correlation. Residual
#' variances and covariances are computed from the DPM stationarity equations so
#' the observed variables have variance 1 across waves when the requested values
#' are admissible.
#'
#' @return
#' An object of class \code{powRICLPM}, upon which \code{summary()},
#' \code{print()}, and \code{plot()} can be used.
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
#'   time_points = 3,
#'   AF_proportion = 0.2,
#'   AF_cor = 0.3,
#'   lagged_effects = lagged_effects,
#'   wave_cor = 0.3,
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
    lagged_effects = NULL,
    wave_cor,
    loadings = NULL,
    skewness = 0,
    kurtosis = 0,
    significance_criterion = 0.05,
    alpha = NULL,
    reps = 20,
    bootstrap_reps = NULL,
    seed = NA,
    constraints = "none",
    estimator = NA,
    Phi = NULL
  ) {

  call_powDPM <- match.call()
  if (!is.null(Phi) && !is.null(lagged_effects)) {
    iabort_renamed_argument_conflict("lagged_effects", "Phi")
  }

  argument_names <- list(
    intraclass_correlation = "AF_proportion",
    RI_cor = "AF_cor",
    within_cor = "wave_cor",
    lagged_effects = if (!is.null(Phi)) "Phi" else "lagged_effects"
  )
  if (!is.null(Phi)) {
    lagged_effects <- Phi
  }

  time_start <- proc.time()

  cli::cli_h2("Checking Argument Input")

  if (!is.null(alpha)) {
    significance_criterion <- icheck_alpha(alpha)
  }

  icheck_target(target_power)
  icheck_T_DPM(time_points)
  icheck_ICC(AF_proportion, arg = argument_names$intraclass_correlation)
  icheck_cor(AF_cor, arg = argument_names$RI_cor)
  icheck_cor(wave_cor, arg = argument_names$within_cor)
  icheck_lagged_effects(lagged_effects, arg = argument_names$lagged_effects)
  icheck_moment(skewness)
  icheck_moment(kurtosis)
  icheck_significance_criterion(significance_criterion)
  icheck_reps(reps)
  estimator <- icheck_estimator(estimator, skewness, kurtosis)
  if (is.character(constraints)) {
    icheck_DPM_constraints(constraints)
  }
  icheck_constraints(constraints, ME = FALSE)
  icheck_loadings(loadings, time_points, "lavaan", constraints, model = "DPM")

  if (!is.null(bootstrap_reps)) {
    cli::cli_alert_warning("The argument {.arg bootstrap_reps} is superseded. Uncertainty regarding simulation estimates is now computed analytically based on Morris et al. (2017).")
  }

  lapply(time_points, function(time_points) {
    loadings_i <- inormalize_loadings(loadings, time_points, model = "DPM")
    lapply(AF_proportion, function(AF_proportion_i) {
      AF_var <- compute_AF_var(AF_proportion_i)
      dpm_values <- compute_DPM_values(
        lagged_effects = lagged_effects,
        wave_cor = wave_cor,
        AF_var = AF_var,
        AF_cov = compute_AF_cov(AF_cor, AF_var),
        loadings = loadings_i
      )
      lapply(seq_len(dim(dpm_values$Psi)[3]), function(i) {
        icheck_Psi(dpm_values$Psi[, , i], model = "DPM")
      })
    })
  })

  if (is.null(sample_size)) {
    icheck_sample_size_search(search_lower, search_upper, search_step)
    sample_size <- seq(search_lower, search_upper, search_step)
  }
  icheck_N(sample_size, time_points, constraints, ME = FALSE, model = "DPM")
  seed <- icheck_seed(seed)

  cli::cli_alert_success("Argument checking complete.")

  irun_power_analysis(
    model = "DPM",
    time_start = time_start,
    call = call_powDPM,
    argument_names = argument_names,
    target_power = target_power,
    sample_size = sample_size,
    time_points = time_points,
    intraclass_correlation = AF_proportion,
    RI_cor = AF_cor,
    lagged_effects = lagged_effects,
    within_cor = wave_cor,
    Psi = NULL,
    reliability = 1,
    loadings = loadings,
    skewness = skewness,
    kurtosis = kurtosis,
    estimate_ME = FALSE,
    significance_criterion = significance_criterion,
    reps = reps,
    bootstrap_reps = bootstrap_reps,
    seed = seed,
    constraints = constraints,
    bounds = FALSE,
    estimator = estimator,
    save_path = NULL,
    software = "lavaan"
  )
}
