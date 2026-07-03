#' Power analysis for the RI-CLPM and STARTS model
#'
#' @description
#' Perform a Monte Carlo power analysis for the random intercept cross-lagged panel model (RI-CLPM) and the stable trait autoregressive trait state model (STARTS). This function computes performance metrics such as bias, mean square error, coverage, power, etc, for all model parameters, and can perform power analyses across multiple experimental conditions simultaneously. Conditions are defined in terms of sample size, number of time points, intraclass correlation, and indicator reliability when available. See "Details" for information on (a) internal data simulation, (b) internal model estimation, (c) `powRICLPM`'s naming conventions of parameters, (d) parallel execution capabilities for speeding up the analysis, and (e) various extensions, such as the option to include measurement errors for data generation and estimation (i.e., the STARTS model), imposing various constraints over time, and many more. Use \code{\link{powDPM}} for power analysis of the dynamic panel model.
#'
#' @param target_power A numeric value between 0 and 1, denoting the targeted power level.
#' @param search_lower A positive \code{integer}, denoting the lower bound of a range of sample sizes.
#' @param search_upper A positive \code{integer}, denoting the upper bound of a range of sample sizes.
#' @param search_step A positive \code{integer}, denoting an increment in sample size.
#' @param sample_size (optional) An \code{integer} (vector), indicating specific sample sizes at which to evaluate power, rather than specifying a range using the \code{search_*} arguments.
#' @param time_points An \code{integer} (vector), indicating number of time points. The RI-CLPM and STARTS specifications need at least 3 waves, with the exact minimum depending on \code{estimate_ME}, \code{constraints}, and whether random-intercept loadings are freely estimated.
#' @param intraclass_correlation A \code{double} (vector) with elements between 0 and 1, denoting the proportion of (true score) variance at the between-unit level in the RI-CLPM. When measurement error is included in the data generating model, the intraclass correlation is computed as the variance of the random intercept factor divided by the true score variance (i.e., controlled for measurement error). When time-varying random-intercept loadings are supplied through \code{loadings}, this value remains the ICC at the first wave because the first loading is fixed to 1. Later loadings are interpreted relative to those first-wave between-unit differences: for example, a loading of 0.9 means the between-unit differences are 0.9 times as large on the score scale, and contribute \eqn{0.9^2} times the first-wave random-intercept variance.
#' @param RI_cor A \code{double} between 0 and 1, denoting the correlation between random intercepts.
#' @param lagged_effects A matrix, with standardized autoregressive effects (on the diagonal) and cross-lagged effects (off-diagonal) in the population. Columns represent predictors and rows represent outcomes.
#' @param within_cor A \code{double} between 0 and 1, denoting the correlation between the within-unit components in the RI-CLPM.
#' @param reliability (optional) A \code{numeric} value, vector, or matrix with elements larger than 0.1 and at most 1, denoting the reliability of the variables (see "Details").
#' @param loadings (optional) A \code{numeric} vector or matrix specifying random-intercept loadings for the RI-CLPM in the \pkg{lavaan} data-generating model (see "Details").
#' @param skewness (optional) A \code{numeric}, denoting the skewness values for the observed variables (see \code{\link[lavaan]{simulateData}}).
#' @param kurtosis (optional) A \code{numeric} value, denoting the excess kurtosis values (i.e., compared to the kurtosis of a normal distribution) for the observed variables (see \code{\link[lavaan]{simulateData}}).
#' @param estimate_ME (optional) A \code{logical}, denoting if measurement error variance should be estimated in the RI-CLPM (see "Details").
#' @param significance_criterion (optional) A \code{double}, denoting the significance criterion.
#' @param alpha (don't use) Deprecated, use `significance_criterion` instead.
#' @param reps A positive \code{integer}, denoting the number of Monte Carlo replications to be used during simulations.
#' @param bootstrap_reps (superseded) Uncertainty regarding simulation estimates is now computed analytically based on Morris et al. (2017). This argument is not used anymore.
#' @param seed An \code{integer} of length 1. If multiple cores are used, a seed will be used to generate a full L'Ecuyer-CMRG seed for all cores.
#' @param constraints (optional) A \code{character} vector, specifying the type of constraints that should be imposed on the estimation model (see "Details").
#' @param bounds (optional) A \code{logical}, denoting if bounded estimation should be used for the latent variable variances in the model (see "Details").
#' @param estimator (optional) A \code{character} string of length 1, denoting the estimator to be used. The default \code{NA} uses \code{ML} for normally generated data and \code{MLR} when \code{skewness} or \code{kurtosis} is nonzero (see "Details").
#' @param save_path A \code{character} string of length 1, naming the directory to save (data) files to (used for validation purposes of this package). Variables are saved in alphabetical and numerical order.
#' @param software A \code{character} string of length, naming which software to use for simulations; either "lavaan" or "Mplus" (see "Details").
#' @param ICC Alternative name for \code{intraclass_correlation}.
#' @param Phi Alternative name for \code{lagged_effects}.
#'
#' @details A rationale for the power analysis strategy implemented in this package can be found in Mulder (2023).
#'
#' \subsection{Data Generation}{Data are generated using \code{\link[lavaan]{simulateData}} from the \pkg{lavaan} package. Based on \code{lagged_effects} and \code{within_cor}, the residual variances and covariances for the within-components at wave 2 and later are computed, such that the within-components themselves have a variance of 1. This implies that the lagged effects in \code{lagged_effects} can be interpreted as standardized effects. By default, all random-intercept loadings in the data-generating model are fixed to 1. The \code{loadings} argument can be used to specify time-varying random-intercept loadings for lavaan data generation, with the first loading fixed to 1 and later loadings interpreted relative to the first occasion. The estimation model frees the corresponding random-intercept loadings only when \code{constraints = "RI_loadings_free"} is supplied. If the data-generating model uses varying random-intercept loadings while the estimation model keeps them fixed, \code{powRICLPM()} prints an informational message before running. With time-varying loadings, \code{intraclass_correlation} still gives the ICC at the first wave because the first loading is fixed to 1. Later loadings describe how strongly those first-wave between-unit differences carry into each later wave; for example, a loading of 0.9 means the between-unit differences are 0.9 times as large on the score scale, and contribute \eqn{0.9^2} times the first-wave random-intercept variance. Consequently, the supplied value is not a wave-invariant ICC.}
#'
#' \subsection{Model Estimation using lavaan}{When \code{software = "lavaan"} (default), generated data are analyzed using \code{\link[lavaan]{lavaan}} from the \pkg{lavaan} package. With the default \code{estimator = NA}, the estimator is maximum likelihood (\code{ML}) for normally generated data and robust maximum likelihood (\code{MLR}) when skewed or kurtosed data are generated (using the \code{skewness} and \code{kurtosis} arguments). Other maximum likelihood based estimators implemented in \href{https://lavaan.ugent.be/tutorial/est.html}{\pkg{lavaan}} can be specified as well. The population parameter values are used as starting values.
#'
#' Parameter estimates from non-converged model solutions are discarded from the results. When \code{bounds = FALSE}, inadmissible parameter estimates from converged solutions (e.g., a negative random intercept variance) are discarded. When \code{bounds = TRUE}, inadmissible parameter estimates are retained following advice by De Jonckere and Rosseel (2022). The results include the minimum estimates for all parameters across replications to diagnose which parameter(s) might be the cause of the inadmissible solution.}
#'
#' \subsection{Using Mplus}{When \code{software = "Mplus"}, Mplus input files will be generated and saved into \code{save_path}. Note that it is not possible to generate skewed or kurtosed data in Mplus via the `powRICLPM` package. Furthermore, bounded estimation and random-intercept loadings are not available in Mplus. The \code{"ME"} constraint can be combined with other Mplus-supported constraints when \code{estimate_ME = TRUE}. }
#'
#' \subsection{Naming Conventions Observed and Latent Variables}{The observed variables in the RI-CLPM are given default names, namely capital letters in alphabetical order, with numbers denoting the measurement occasion. For example, for a bivariate RICLPM with 3 time points, we observe \code{A1}, \code{A2}, \code{A3}, \code{B1}, \code{B2}, and \code{B3}. Their within-components are denoted by \code{wA1}, \code{wA2}, ..., \code{wB3}, respectively. The between-components have \code{RI_} prepended to the variable name, resulting in \code{RI_A} and \code{RI_B}.
#'
#' Parameters are denoted using \pkg{lavaan} model syntax (see \href{https://lavaan.ugent.be/tutorial/syntax1.html}{the \pkg{lavaan} website}). For example, the random intercept variances are denoted by \code{RI_A~~RI_A} and \code{RI_B~~RI_B}, the cross-lagged effects at the first wave as \code{wB2~wA1} and \code{wA2~wB1}, and the autoregressive effects as \code{wA2~wA1} and \code{wB2~wB1}. When factor loadings are freely estimated using \code{constraints = "RI_loadings_free"}, the freed loadings are denoted by the corresponding \pkg{lavaan} loading syntax, for example \code{RI_A=~A2}. Use \code{give(object, "names")} to extract parameter names from the \code{powRICLPM} object.}
#'
#' \subsection{Parallel Processing and Progress Bar}{To speed up the analysis, power analysis for multiple experimental conditions can be executed in parallel. This has been implemented using \pkg{future}. By default the analysis is executed sequentially (i.e., single-core). Parallel execution (i.e., multicore) can be setup using \code{\link[future]{plan}}, for example \code{plan(multisession)}. For more information and options, see \url{https://future.futureverse.org/articles/future-1-overview.html#controlling-how-futures-are-resolved}.
#'
#' A progress bar displaying the status of the power analysis has been implemented using \pkg{progressr}. By default, a simple progress bar will be shown. For more information on how to control this progress bar and several other notification options (e.g., auditory notifications), see \url{https://progressr.futureverse.org}.}
#'
#' \subsection{Extension: Measurement Errors (STARTS model)}{Including measurement error to the RI-CLPM makes the model equivalent to the bivariate STARTS model by Kenny and Zautra (2001) without constraints over time. Measurement error can be added to the generated data through the \code{reliability} argument. Setting \code{reliability = 0.8} implies that 80 percent is true score variance, and 20 percent is measurement error variance. A vector specifies multiple reliability conditions, each applied to both variables. A matrix with two rows specifies reliability conditions separately for variables A and B, with each column defining one experimental condition. \code{intraclass_correlation} then denotes the proportion of \emph{true score variance} captured by the random intercept factors. Estimating measurement errors (i.e., the STARTS model) is done by setting \code{estimate_ME = TRUE}. If generated measurement error is ignored by the estimation model, \code{powRICLPM()} prints an informational message before running.}
#'
#' \subsection{Extension: Imposing Constraints}{The following options can be supplied to the estimation model using the \code{constraints} argument. By default, \code{constraints = "none"}, and no equality or time-invariance constraints are imposed.
#'
#' \itemize{
#'   \item \code{"lagged"}: Time-invariant autoregressive and cross-lagged effects.
#'   \item \code{"residuals"}: Time-invariant within-unit residual variances and covariances.
#'   \item \code{"within"}: Time-invariant lagged effects and within-unit residual variances and covariances. This is equivalent to \code{constraints = c("lagged", "residuals")}.
#'   \item \code{"stationarity"}: Constraints such that at the within-unit level a stationary process is estimated. This includes time-invariant lagged effects and constraints on the residual variances. It cannot be combined with \code{"lagged"} or \code{"residuals"}.
#'   \item \code{"ME"}: Time-invariant measurement error variances. Only possible when \code{estimate_ME = TRUE}.
#'   \item \code{"RI_loadings_free"}: Freely estimated random-intercept factor loadings.
#' }
#'
#' The minimum number of waves is checked against the fitted RI-CLPM or STARTS specification:
#'
#' \tabular{llll}{
#'   \strong{Estimated ME} \tab \strong{Constraints} \tab \strong{Free RI loadings} \tab \strong{Minimum waves}\cr
#'   FALSE \tab none \tab FALSE \tab 3\cr
#'   FALSE \tab lagged \tab FALSE \tab 3\cr
#'   FALSE \tab residuals \tab FALSE \tab 3\cr
#'   FALSE \tab lagged + residuals/within \tab FALSE \tab 3\cr
#'   FALSE \tab stationarity \tab FALSE \tab 3\cr
#'   FALSE \tab none \tab TRUE \tab 4\cr
#'   FALSE \tab lagged \tab TRUE \tab 3\cr
#'   FALSE \tab residuals \tab TRUE \tab 3\cr
#'   FALSE \tab lagged + residuals/within \tab TRUE \tab 3\cr
#'   FALSE \tab stationarity \tab TRUE \tab 3\cr
#'   TRUE \tab none \tab FALSE \tab 4\cr
#'   TRUE \tab none + ME \tab FALSE \tab 4\cr
#'   TRUE \tab lagged \tab FALSE \tab 4\cr
#'   TRUE \tab lagged + ME \tab FALSE \tab 3\cr
#'   TRUE \tab residuals \tab FALSE \tab 4\cr
#'   TRUE \tab residuals + ME \tab FALSE \tab 3\cr
#'   TRUE \tab lagged + residuals/within \tab FALSE \tab 3\cr
#'   TRUE \tab lagged + residuals/within + ME \tab FALSE \tab 3\cr
#'   TRUE \tab stationarity \tab FALSE \tab 3\cr
#'   TRUE \tab stationarity + ME \tab FALSE \tab 3\cr
#'   TRUE \tab none \tab TRUE \tab 5\cr
#'   TRUE \tab none + ME \tab TRUE \tab 4\cr
#'   TRUE \tab lagged \tab TRUE \tab 4\cr
#'   TRUE \tab lagged + ME \tab TRUE \tab 4\cr
#'   TRUE \tab residuals \tab TRUE \tab 4\cr
#'   TRUE \tab residuals + ME \tab TRUE \tab 4\cr
#'   TRUE \tab lagged + residuals/within \tab TRUE \tab 4\cr
#'   TRUE \tab lagged + residuals/within + ME \tab TRUE \tab 3\cr
#'   TRUE \tab stationarity \tab TRUE \tab 4\cr
#'   TRUE \tab stationarity + ME \tab TRUE \tab 3
#' }
#'
#' The factor loadings are a special case. In the standard RI-CLPM specification, random-intercept loadings are fixed to 1 so that each observed variable loads equally on its random intercept. The option \code{constraints = "RI_loadings_free"} relaxes this assumption in the estimation model by freely estimating the random-intercept factor loadings from the second wave onward.
#'
#' Compatible options can be combined in a character vector. For example, \code{constraints = c("lagged", "RI_loadings_free")} constrains the lagged effects over time and freely estimates random-intercept factor loadings. Similarly, \code{constraints = c("lagged", "residuals")} imposes the same RI-CLPM constraints as \code{constraints = "within"}. Some options cannot be combined because they impose overlapping or incompatible restrictions. For example, \code{constraints = "stationarity"} already imposes restrictions on the residual variances, so it cannot be combined with \code{constraints = "residuals"}. Similarly, because \code{constraints = "within"} already combines lagged and residual constraints, it cannot be combined with \code{"lagged"} or \code{"residuals"} in the same constraint vector.
#'
#'}
#'
#' \subsection{Extension: Bounded Estimation}{Bounded estimation is useful to avoid nonconvergence in small samples. Here, automatic wide bounds are used as advised by De Jonckere and Rosseel (2022), see \code{optim.bounds} in \code{\link[lavaan]{lavOptions}}. This option can be used with lavaan specifications, including specifications with equality or time-invariance constraints. Bounded estimation is not available for \code{software = "Mplus"}.}
#'
#' @return
#' An object of class `powRICLPM`, upon which \code{summary()}, \code{print()}, and \code{plot()} can be used. The returned object is a \code{list} with a \code{conditions} and \code{session} element. \code{condition} itself is a \code{list} of experimental conditions, where each element is again a \code{list} containing the input and output of the power analysis for that particular experimental condition. \code{session} is a \code{list} containing information common to all experimental conditions.
#'
#' @author Jeroen D. Mulder \email{j.d.mulder@@uu.nl}
#'
#' @references
#' De Jonckere, J., & Rosseel, Y. (2022). Using bounded estimation to avoid nonconvergence in small sample structural equation modeling. \emph{Structural Equation Modeling}, \emph{29}(3), 412-427. \doi{10.1080/10705511.2021.1982716}
#'
#' Kenny, D. A., & Zautra, A. (2001). Trait–state models for longitudinal data. \emph{New methods for the analysis of change} (pp. 243–263). American Psychological Association. \doi{10.1037/10409-008}
#'
#' Mulder, J. D. (2023). Power analysis for the random intercept cross-lagged panel model using the \emph{powRICLPM} R-package. \emph{Structural Equation Modeling: A Multidisciplinary Journal}, \emph{30}(4), 645-658. \doi{10.1080/10705511.2022.2122467}
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{powDPM}}: Run power analysis for the dynamic panel model.
#'   \item \code{\link{summary.powRICLPM}}: Summarize the setup of \code{powRICLPM} object.
#'   \item \code{\link{give}}: Extract information from \code{powRICLPM} objects.
#'   \item \code{\link{plot.powRICLPM}}: Visualize results \code{powRICLPM} object for a specific parameter.
#' }
#'
#' @examples
#' # Define population parameters for lagged effects
#' lagged_effects <- matrix(c(.4, .1, .2, .3), ncol = 2, byrow = TRUE)
#'
#' \dontrun{
#' # (optional) Set up parallel computing (i.e., multicore, speeding up the analysis)
#' library(future)
#' library(progressr)
#'
#' local({
#' old_plan <- future::plan()
#' on.exit(future::plan(old_plan), add = TRUE)
#' future::plan(future::multisession)
#'
#' # Run analysis (`reps` is small, because this is an example)
#' with_progress({
#'   out_preliminary <- powRICLPM(
#'     target_power = 0.8,
#'     search_lower = 500,
#'     search_upper = 700,
#'     search_step = 100,
#'     time_points = c(3, 4),
#'     intraclass_correlation = c(0.4, 0.6),
#'     reliability = 0.8,
#'     RI_cor = 0.3,
#'     lagged_effects = lagged_effects,
#'     within_cor = 0.3,
#'     reps = 100,
#'     seed = 1234
#'   )
#' })
#' })
#' }
#'
#' @importFrom future.apply future_lapply
#' @importFrom future plan multisession sequential
#' @export
powRICLPM <- function(
    target_power = 0.8,
    search_lower = NULL,
    search_upper = NULL,
    search_step = 20,
    sample_size = NULL,
    time_points,
    intraclass_correlation = NULL,
    RI_cor = NULL,
    lagged_effects = NULL,
    within_cor = NULL,
    reliability = 1,
    loadings = NULL,
    skewness = 0,
    kurtosis = 0,
    estimate_ME = FALSE,
    significance_criterion = 0.05,
    alpha = NULL,
    reps = 20,
    bootstrap_reps = NULL,
    seed = NA,
    constraints = "none",
    bounds = FALSE,
    estimator = NA,
    save_path = NULL,
    software = "lavaan",
    ICC = NULL,
    Phi = NULL
  ) {

  # Get call
  call_powRICLPM <- match.call()
  reliability_expr <- call_powRICLPM$reliability
  if (!is.null(ICC) && !is.null(intraclass_correlation)) {
    iabort_renamed_argument_conflict("intraclass_correlation", "ICC")
  }
  if (!is.null(Phi) && !is.null(lagged_effects)) {
    iabort_renamed_argument_conflict("lagged_effects", "Phi")
  }

  argument_names <- list(
    intraclass_correlation = if (!is.null(ICC)) {
      "ICC"
    } else {
      "intraclass_correlation"
    },
    RI_cor = "RI_cor",
    within_cor = "within_cor",
    lagged_effects = if (!is.null(Phi)) "Phi" else "lagged_effects"
  )
  if (!is.null(ICC)) {
    intraclass_correlation <- ICC
  }
  if (!is.null(Phi)) {
    lagged_effects <- Phi
  }

  ICC <- intraclass_correlation

  # Start time
  time_start <- proc.time()

  # Inform user that input checking is starting
  cli::cli_h2("Checking Argument Input")

  # Input checker (deprecated)
  if (! is.null(alpha)) {
    significance_criterion <- icheck_alpha(alpha)
  }

  # Input checkers (icheck) I
  icheck_target(target_power)
  icheck_T(time_points, estimate_ME)
  icheck_ICC(ICC, arg = argument_names$intraclass_correlation)
  icheck_cor(RI_cor, arg = argument_names$RI_cor)
  icheck_cor(within_cor, arg = argument_names$within_cor)
  icheck_lagged_effects(lagged_effects, arg = argument_names$lagged_effects)
  icheck_moment(skewness)
  icheck_moment(kurtosis)
  icheck_significance_criterion(significance_criterion)
  icheck_ME(estimate_ME)
  icheck_reps(reps)
  estimator <- icheck_estimator(estimator, skewness, kurtosis)
  save_path <- icheck_path(save_path, software)
  icheck_software(software, skewness, kurtosis)
  if (!is.null(reliability_expr)) {
    icheck_reliability_matrix_call(
      expr = reliability_expr,
      env = parent.frame(),
      call = rlang::caller_env()
    )
  }
  icheck_rel(reliability, time_points, software)
  icheck_constraints(constraints, estimate_ME)
  icheck_bounds(bounds, constraints, software)
  icheck_loadings(loadings, time_points, software, constraints, model = "RICLPM")
  icheck_constraints_software(constraints, software)
  icheck_RICLPM_identification(time_points, estimate_ME, constraints)

  inote_ignored_input_names(sample_size, "sample_size", "Sample-size values")
  inote_ignored_input_names(time_points, "time_points", "Time-point values")
  inote_ignored_input_names(ICC, argument_names$intraclass_correlation, "Intraclass-correlation values")
  inote_ignored_input_names(RI_cor, argument_names$RI_cor, "Correlation values")
  inote_ignored_input_names(within_cor, argument_names$within_cor, "Correlation values")
  inote_ignored_input_names(lagged_effects, argument_names$lagged_effects, "Lagged-effect values")
  inote_ignored_input_names(reliability, "reliability", "Reliability values")
  inote_ignored_input_names(loadings, "loadings", "Loading values")
  inote_ignored_input_names(constraints, "constraints", "Constraint values")

  # powRICLPM updates
  if (!is.null(bootstrap_reps)) {
    cli::cli_alert_warning("The argument {.arg bootstrap_reps} is superseded. Uncertainty regarding simulation estimates is now computed analytically based on Morris et al. (2017).")
  }

  # Compute population parameter values for data generation
  Psi <- compute_Psi(lagged_effects, within_cor)
  icheck_Psi(Psi)

  # Get candidate sample sizes
  if (is.null(sample_size)) {
    icheck_sample_size_search(search_lower, search_upper, search_step)
    sample_size <- seq(search_lower, search_upper, search_step)
  }
  icheck_N(sample_size, time_points, constraints, estimate_ME, model = "RICLPM")
  seed <- icheck_seed(seed)

  # Inform user that input check is complete
  inote_custom_loadings_interpretation(
    loadings,
    time_points,
    model = "RICLPM",
    proportion_label = argument_names$intraclass_correlation
  )
  cli::cli_alert_success("Argument checking complete.")

  irun_power_analysis(
    model = "RICLPM",
    time_start = time_start,
    call = call_powRICLPM,
    argument_names = argument_names,
    target_power = target_power,
    sample_size = sample_size,
    time_points = time_points,
    intraclass_correlation = ICC,
    RI_cor = RI_cor,
    lagged_effects = lagged_effects,
    within_cor = within_cor,
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
  )
}


irun_power_analysis <- function(
  model,
  time_start,
  call,
  argument_names,
  target_power,
  sample_size,
  time_points,
  intraclass_correlation,
  RI_cor,
  lagged_effects,
  within_cor,
  AF_proportion = NULL,
  AF_cor = NULL,
  dynamics_cor = NULL,
  Psi,
  reliability,
  loadings,
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
  conditions <- if (identical(model, "DPM")) {
    create_conditions_DPM(
      target_power = target_power,
      sample_size = sample_size,
      time_points = time_points,
      AF_proportion = AF_proportion,
      AF_cor = AF_cor,
      lagged_effects = lagged_effects,
      dynamics_cor = dynamics_cor,
      reliability = reliability,
      loadings = loadings,
      skewness = skewness,
      kurtosis = kurtosis,
      estimate_ME = estimate_ME,
      significance_criterion = significance_criterion,
      reps = reps,
      seed = seed,
      constraints = constraints,
      bounds = bounds,
      estimator = estimator,
      software = software
    )
  } else {
    create_conditions(
      model = model,
      target_power = target_power,
      sample_size = sample_size,
      time_points = time_points,
      intraclass_correlation = intraclass_correlation,
      RI_cor = RI_cor,
      lagged_effects = lagged_effects,
      within_cor = within_cor,
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
    )
  }

  misspecification <- if (identical(model, "DPM")) {
    combine_DPM_misspecification(lapply(conditions, function(condition) {
      condition$misspecification
    }))
  } else if (identical(model, "RICLPM")) {
    combine_RICLPM_misspecification(lapply(conditions, function(condition) {
      condition$misspecification
    }))
  } else {
    NULL
  }
  if (!is.null(misspecification)) {
    confirm_restrictive_misspecification(misspecification)
  }

  if (software == "lavaan") {
    seed <- icheck_seed(seed)

    # Inform user that simulations in lavaan have started
    cli::cli_h2("\nPerforming Simulations")

    # Prepare progress bar
    p <- progressr::progressor(steps = (length(conditions) * reps))

    # Run Monte Carlo simulation for each condition
    conditions <- future.apply::future_lapply(
      conditions,
      FUN = function(x) {
        run_condition_monteCarlo(
          condition = x,
          p = p,
          bounds = bounds,
          estimator = estimator,
          reps = reps,
          save_path = save_path
        )
      },
      future.seed = seed,
      future.globals = TRUE
    )

    # Inform user that are complete
    cli::cli_alert_success("Simulations complete.")

    # End time
    time_end <- proc.time()
    time_taken <- time_end - time_start

    # Create powRICLPM object, combining conditions and general session info
    out <- list(
      conditions = conditions,
      session = list(
        estimate_ME = estimate_ME,
        model = model,
        reps = reps,
        target_power = target_power,
        constraints = constraints,
        bounds = bounds,
        estimator = estimator,
        software = software,
        save_path = save_path,
        misspecified_restrictive = if (is.null(misspecification)) FALSE else misspecification$restrictive,
        misspecified_general = if (is.null(misspecification)) FALSE else misspecification$general,
        misspecification_reasons = if (is.null(misspecification)) character() else misspecification$reasons,
        misspecification_restrictive_reasons = if (is.null(misspecification)) character() else misspecification$restrictive_reasons,
        misspecification_general_reasons = if (is.null(misspecification)) character() else misspecification$general_reasons,
        time_taken = time_taken,
        version = utils::packageVersion("powRICLPM"),
        call = call,
        argument_names = argument_names
      )
    )

    if (identical(model, "DPM")) {
      class(out) <- c("powDPM", "powRICLPM", class(out))
    } else {
      class(out) <- c("powRICLPM", class(out))
    }

    return(out)
  } else if (software == "Mplus") {
    # Inform user that simulations in Mplus have started
    cli::cli_h2("\nPerforming Simulations")

    # Inform user of results
    print.powRICLPM.Mplus(
      conditions,
      save_path = save_path,
      icc_label = if (identical(argument_names$intraclass_correlation, "ICC")) {
        "ICC"
      } else {
        "Intraclass correlation"
      }
    )

    invisible()
  }
}


#' (superseded) Create Mplus Syntax for RI-CLPM Power Analysis
#'
#' @description
#' The function `powRICLPM_Mplus` has been superseded and is no longer available in the current version of the package.
#'
#' @details
#' The functionality previously provided by `powRICLPM_Mplus` has been replaced by `powRICLPM`, where you can now set `software = "Mplus"`.
#'
#' @param ... (don't use)
#'
#' @seealso [powRICLPM()]
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#'   lagged_effects <- matrix(c(.4, .1, .2, .3), ncol = 2, byrow = TRUE)
#'
#'   # Use `software = "Mplus"` to setup power analysis for Mplus
#'   out_preliminary <- powRICLPM(
#'     target_power = 0.8,
#'     search_lower = 500,
#'     search_upper = 700,
#'     search_step = 100,
#'     time_points = c(3, 4),
#'     intraclass_correlation = c(0.4, 0.6),
#'     reliability = 0.8,
#'     RI_cor = 0.3,
#'     lagged_effects = lagged_effects,
#'     within_cor = 0.3,
#'     reps = 1000,
#'     seed = 1234,
#'     software = "Mplus"
#'   )
#' }
powRICLPM_Mplus <- function(...) {
  .Deprecated(msg = "The function 'powRICLPM_Mplus' is deprecated and is removed since version 0.2.0. Please use 'powRICLPM()' instead, and set the `software` argument to 'Mplus'.")
}

