#' Power analysis for the RI-CLPM, STARTS model, and DPM
#'
#' @description
#' Perform a Monte Carlo power analysis for the random intercept cross-lagged panel model (RI-CLPM), the stable trait autoregressive trait state model (STARTS), and the dynamic panel model (DPM). This function computes performance metrics such as bias, mean square error, coverage, power, etc, for all model parameters, and can perform power analyses across multiple experimental conditions simultaneously. Conditions are defined in terms of sample size, number of time points, stable-factor variance proportion (`intraclass_correlation` for the RI-CLPM or `AF_proportion` for the DPM), and indicator reliability when available. See "Details" for information on (a) internal data simulation, (b) internal model estimation, (c) `powRICLPM`'s naming conventions of parameters, (d) parallel execution capabilities for speeding up the analysis, and (e) various extensions, such as the option to include measurement errors for data generation and estimation (i.e., the STARTS model), imposing various constraints over time, and many more.
#'
#' @param model A \code{character} string of length 1, denoting the model to use. The default \code{"RICLPM"} keeps the random-intercept cross-lagged panel model. Use \code{"DPM"} for the dynamic panel model.
#' @param target_power A numeric value between 0 and 1, denoting the targeted power level.
#' @param search_lower A positive \code{integer}, denoting the lower bound of a range of sample sizes.
#' @param search_upper A positive \code{integer}, denoting the upper bound of a range of sample sizes.
#' @param search_step A positive \code{integer}, denoting an increment in sample size.
#' @param sample_size (optional) An \code{integer} (vector), indicating specific sample sizes at which to evaluate power, rather than specifying a range using the \code{search_*} arguments.
#' @param time_points An \code{integer} (vector), indicating number of time points. The RI-CLPM needs at least 3 waves, STARTS estimation needs at least 4 waves, and the DPM needs at least 2 waves.
#' @param intraclass_correlation A \code{double} (vector) with elements between 0 and 1, denoting the proportion of (true score) variance at the between-unit level in the RI-CLPM. When measurement error is included in the data generating model, the intraclass correlation is computed as the variance of the random intercept factor divided by the true score variance (i.e., controlled for measurement error). When time-varying random-intercept loadings are supplied through \code{loadings}, this value remains the ICC at the first wave because the first loading is fixed to 1. Later loadings are interpreted relative to those first-wave between-unit differences: for example, a loading of 0.9 means the between-unit differences are 0.9 times as large on the score scale, and contribute \eqn{0.9^2} times the first-wave random-intercept variance.
#' @param AF_proportion A \code{double} (vector) with elements between 0 and 1, denoting the proportion of observed variance attributed to the accumulating factors in the DPM. This is the DPM counterpart of \code{intraclass_correlation}.
#' @param RI_cor A \code{double} between 0 and 1, denoting the correlation between random intercepts.
#' @param AF_cor A \code{double} between 0 and 1, denoting the correlation between accumulating factors in the DPM.
#' @param lagged_effects A matrix, with standardized autoregressive effects (on the diagonal) and cross-lagged effects (off-diagonal) in the population. Columns represent predictors and rows represent outcomes.
#' @param within_cor A \code{double} between 0 and 1, denoting the correlation between the within-unit components in the RI-CLPM.
#' @param wave_cor A \code{double} between 0 and 1, denoting the observed wave-level correlation in the DPM. This is the DPM counterpart of \code{within_cor}.
#' @param reliability (optional) A \code{numeric} value, vector, or matrix with elements larger than 0.1 and at most 1, denoting the reliability of the variables (see "Details"). This extension is not available for the DPM.
#' @param loadings (optional) A \code{numeric} vector or matrix specifying random-intercept loadings for the RI-CLPM or accumulating-factor loadings for the DPM in the \pkg{lavaan} data-generating model (see "Details").
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
#' @param software A \code{character} string of length, naming which software to use for simulations; either "lavaan" or "Mplus" (see "Details"). The DPM is currently available only with \code{software = "lavaan"}.
#' @param ICC Alternative name for \code{intraclass_correlation}.
#' @param Phi Alternative name for \code{lagged_effects}.
#'
#' @details A rationale for the power analysis strategy implemented in this package can be found in Mulder (2023).
#'
#' \subsection{Data Generation}{Data are generated using \code{\link[lavaan]{simulateData}} from the \pkg{lavaan} package. Based on \code{lagged_effects} and \code{within_cor}, the residual variances and covariances for the within-components at wave 2 and later are computed, such that the within-components themselves have a variance of 1. This implies that the lagged effects in \code{lagged_effects} can be interpreted as standardized effects. By default, all random-intercept loadings in the data-generating model are fixed to 1. The \code{loadings} argument can be used to specify time-varying random-intercept loadings for lavaan data generation, with the first loading fixed to 1 and later loadings interpreted relative to the first occasion. Supplying \code{loadings} requires \code{constraints = "RI_loadings_free"} or \code{constraints = "loadings_free"} so that the estimation model also frees the corresponding random-intercept loadings. With time-varying loadings, \code{intraclass_correlation} still gives the ICC at the first wave because the first loading is fixed to 1. Later loadings describe how strongly those first-wave between-unit differences carry into each later wave; for example, a loading of 0.9 means the between-unit differences are 0.9 times as large on the score scale, and contribute \eqn{0.9^2} times the first-wave random-intercept variance. Consequently, the supplied value is not a wave-invariant ICC.
#'
#' For \code{model = "DPM"}, the accumulating factors load on waves 2 through T and covary with the first wave. The DPM uses \code{AF_proportion} as the accumulating-factor variance and \code{wave_cor} as the target observed wave-level correlation. Residual variances and covariances are computed from the DPM stationarity equations so the observed variables have variance 1 across waves when the requested values are admissible. DPM \code{loadings} must include one value per wave, with an explicit first-wave \code{NA}. The wave-2 loading is fixed to 1.}
#'
#' \subsection{Model Estimation using lavaan}{When \code{software = "lavaan"} (default), generated data are analyzed using \code{\link[lavaan]{lavaan}} from the \pkg{lavaan} package. With the default \code{estimator = NA}, the estimator is maximum likelihood (\code{ML}) for normally generated data and robust maximum likelihood (\code{MLR}) when skewed or kurtosed data are generated (using the \code{skewness} and \code{kurtosis} arguments). Other maximum likelihood based estimators implemented in \href{https://lavaan.ugent.be/tutorial/est.html}{\pkg{lavaan}} can be specified as well. The population parameter values are used as starting values.
#'
#' Parameter estimates from non-converged model solutions are discarded from the results. When \code{bounds = FALSE}, inadmissible parameter estimates from converged solutions (e.g., a negative random intercept variance) are discarded. When \code{bounds = TRUE}, inadmissible parameter estimates are retained following advice by De Jonckere and Rosseel (2022). The results include the minimum estimates for all parameters across replications to diagnose which parameter(s) might be the cause of the inadmissible solution.}
#'
#' \subsection{Using Mplus}{When \code{software = "Mplus"}, Mplus input files will be generated and saved into \code{save_path}. Note that it is not possible to generate skewed or kurtosed data in Mplus via the `powRICLPM` package. Furthermore, bounded estimation, time-varying reliability, time-varying random-intercept loadings (both for data generation and estimation), and the DPM are not available in Mplus. The \code{"ME"} constraint can be combined with other Mplus-supported constraints when \code{estimate_ME = TRUE}. }
#'
#' \subsection{Naming Conventions Observed and Latent Variables}{The observed variables in the RI-CLPM are given default names, namely capital letters in alphabetical order, with numbers denoting the measurement occasion. For example, for a bivariate RICLPM with 3 time points, we observe \code{A1}, \code{A2}, \code{A3}, \code{B1}, \code{B2}, and \code{B3}. Their within-components are denoted by \code{wA1}, \code{wA2}, ..., \code{wB3}, respectively. The between-components have \code{RI_} prepended to the variable name, resulting in \code{RI_A} and \code{RI_B}.
#'
#' Parameters are denoted using \pkg{lavaan} model syntax (see \href{https://lavaan.ugent.be/tutorial/syntax1.html}{the \pkg{lavaan} website}). For example, the random intercept variances are denoted by \code{RI_A~~RI_A} and \code{RI_B~~RI_B}, the cross-lagged effects at the first wave as \code{wB2~wA1} and \code{wA2~wB1}, and the autoregressive effects as \code{wA2~wA1} and \code{wB2~wB1}. DPM accumulating factors are denoted \code{AF_A} and \code{AF_B}, and DPM lagged effects are observed-level regressions such as \code{B2~A1}. When factor loadings are freely estimated using \code{constraints = "loadings_free"} or \code{constraints = "RI_loadings_free"}, the freed loadings are denoted by the corresponding \pkg{lavaan} loading syntax, for example \code{RI_A=~A2} in the RI-CLPM or \code{AF_A=~A3} in the DPM. Use \code{give(object, "names")} to extract parameter names from the \code{powRICLPM} object.}
#'
#' \subsection{Parallel Processing and Progress Bar}{To speed up the analysis, power analysis for multiple experimental conditions can be executed in parallel. This has been implemented using \pkg{future}. By default the analysis is executed sequentially (i.e., single-core). Parallel execution (i.e., multicore) can be setup using \code{\link[future]{plan}}, for example \code{plan(multisession, workers = 4)}. For more information and options, see \url{https://future.futureverse.org/articles/future-1-overview.html#controlling-how-futures-are-resolved}.
#'
#' A progress bar displaying the status of the power analysis has been implemented using \pkg{progressr}. By default, a simple progress bar will be shown. For more information on how to control this progress bar and several other notification options (e.g., auditory notifications), see \url{https://progressr.futureverse.org}.}
#'
#' \subsection{Extension: Measurement Errors (STARTS model)}{Including measurement error to the RI-CLPM makes the model equivalent to the bivariate STARTS model by Kenny and Zautra (2001) without constraints over time. Measurement error can be added to the generated data through the \code{reliability} argument. Setting \code{reliability = 0.8} implies that 80 percent is true score variance, and 20 percent is measurement error variance. A vector of length \code{time_points} specifies reliabilities that vary over time and are the same for both variables. A matrix with two rows and one column per time point specifies reliabilities separately for variables A and B. Vector and matrix reliability specifications are only available with \code{software = "lavaan"} and one value of \code{time_points}. To compare multiple reliability patterns, use multiple calls to \code{powRICLPM()}. \code{intraclass_correlation} then denotes the proportion of \emph{true score variance} captured by the random intercept factors. Estimating measurement errors (i.e., the STARTS model) is done by setting \code{estimate_ME = TRUE}.}
#'
#' \subsection{Extension: Imposing Constraints}{The following options can be supplied to the estimation model using the \code{constraints} argument. By default, \code{constraints = "none"}, and no equality or time-invariance constraints are imposed.
#'
#' \itemize{
#'   \item \code{"lagged"}: Time-invariant autoregressive and cross-lagged effects.
#'   \item \code{"residuals"}: Time-invariant within-unit residual variances and covariances.
#'   \item \code{"within"}: Time-invariant lagged effects and within-unit residual variances and covariances. This is equivalent to \code{constraints = c("lagged", "residuals")}.
#'   \item \code{"stationarity"}: Constraints such that at the within-unit level a stationary process is estimated. This includes time-invariant lagged effects and constraints on the residual variances. It cannot be combined with \code{"lagged"} or \code{"residuals"}.
#'   \item \code{"ME"}: Time-invariant measurement error variances. Only possible when \code{estimate_ME = TRUE}; not available for the DPM.
#'   \item \code{"loadings_free"}: Freely estimated factor loadings. For backwards compatibility, \code{"RI_loadings_free"} remains available for RI-CLPM random-intercept loadings.
#' }
#'
#' The factor loadings are a special case. In the standard RI-CLPM specification, random-intercept loadings are fixed to 1 so that each observed variable loads equally on its random intercept. The option \code{constraints = "loadings_free"} or \code{constraints = "RI_loadings_free"} relaxes this assumption in the estimation model by freely estimating the random-intercept factor loadings from the second wave onward. For the DPM, \code{constraints = "loadings_free"} frees accumulating-factor loadings from wave 3 onward while keeping the wave-2 loading fixed to 1.
#'
#' Compatible options can be combined in a character vector. For example, \code{constraints = c("lagged", "loadings_free")} constrains the lagged effects over time and freely estimates factor loadings. Similarly, \code{constraints = c("lagged", "residuals")} imposes the same RI-CLPM constraints as \code{constraints = "within"}. Some options cannot be combined because they impose overlapping or incompatible restrictions. For example, \code{constraints = "stationarity"} already imposes restrictions on the residual variances, so it cannot be combined with \code{constraints = "residuals"}. Similarly, because \code{constraints = "within"} already combines lagged and residual constraints, it cannot be combined with \code{"lagged"} or \code{"residuals"} in the same constraint vector. The shorthand \code{"within"} is RI-CLPM terminology and is not available for the DPM; use \code{c("lagged", "residuals")} instead.
#'
#'}
#'
#' \subsection{Extension: Bounded Estimation}{Bounded estimation is useful to avoid nonconvergence in small samples. Here, automatic wide bounds are used as advised by De Jonckere and Rosseel (2022), see \code{optim.bounds} in \code{\link[lavaan]{lavOptions}}. This option can only be used when no equality or time-invariance constraints are imposed on the estimation model, except for \code{constraints = "loadings_free"} or \code{constraints = "RI_loadings_free"} since these free factor loadings rather than constraining parameters. Bounded estimation is not yet available for the DPM.}
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
#'   \item \code{\link{summary.powRICLPM}}: Summarize the setup of \code{powRICLPM} object.
#'   \item \code{\link{give}}: Extract information from \code{powRICLPM} objects.
#'   \item \code{\link{plot.powRICLPM}}: Visualize results \code{powRICLPM} object for a specific parameter.
#' }
#'
#' @examples
#' # Define population parameters for lagged effects
#' lagged_effects <- matrix(c(.4, .1, .2, .3), ncol = 2, byrow = TRUE)
#'
#' # (optional) Set up parallel computing (i.e., multicore, speeding up the analysis)
#' library(future)
#' library(progressr)
#' future::plan(multisession, workers = 6)
#'
#' \dontrun{
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
#' }
#'
#' \dontshow{
#' ## Shut down parallel workers (done for sake of example, normally not needed)
#' future::plan("sequential")
#' }
#'
#' @importFrom future.apply future_lapply
#' @importFrom future plan multisession sequential
#' @export
powRICLPM <- function(
    model = "RICLPM",
    target_power = 0.8,
    search_lower = NULL,
    search_upper = NULL,
    search_step = 20,
    sample_size = NULL,
    time_points,
    intraclass_correlation = NULL,
    AF_proportion = NULL,
    RI_cor = NULL,
    AF_cor = NULL,
    lagged_effects = NULL,
    within_cor = NULL,
    wave_cor = NULL,
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
  model <- icheck_model(model)

  if (model == "DPM") {
    icheck_DPM_aliases(AF_proportion, intraclass_correlation, ICC,
                       RI_cor, within_cor)
  }
  if (!is.null(ICC) && !is.null(intraclass_correlation)) {
    iabort_renamed_argument_conflict("intraclass_correlation", "ICC")
  }
  if (!is.null(AF_proportion) &&
      (!is.null(intraclass_correlation) || !is.null(ICC))) {
    iabort_renamed_argument_conflict("AF_proportion", "intraclass_correlation/ICC")
  }
  if (!is.null(AF_cor) && !is.null(RI_cor)) {
    iabort_renamed_argument_conflict("AF_cor", "RI_cor")
  }
  if (!is.null(wave_cor) && !is.null(within_cor)) {
    iabort_renamed_argument_conflict("wave_cor", "within_cor")
  }
  if (!is.null(Phi) && !is.null(lagged_effects)) {
    iabort_renamed_argument_conflict("lagged_effects", "Phi")
  }

  if (model == "RICLPM" && (!is.null(AF_proportion) || !is.null(AF_cor) || !is.null(wave_cor))) {
    cli::cli_abort(
      c(
        "These arguments look like a DPM setup, but `model` defaults to `'RICLPM'`:",
        i = "Use `model = 'DPM'` with `AF_proportion`, `AF_cor`, and `wave_cor`.",
        i = "For the RI-CLPM, use `intraclass_correlation`, `RI_cor`, and `within_cor`.",
        x = paste0(
          "DPM argument(s) supplied: ",
          paste(c("AF_proportion", "AF_cor", "wave_cor")[
            !vapply(list(AF_proportion, AF_cor, wave_cor), is.null, logical(1))
          ], collapse = ", "),
          "."
        )
      )
    )
  }

  argument_names <- list(
    intraclass_correlation = if (model == "DPM") {
      "AF_proportion"
    } else if (!is.null(ICC)) {
      "ICC"
    } else {
      "intraclass_correlation"
    },
    RI_cor = if (model == "DPM") "AF_cor" else "RI_cor",
    within_cor = if (model == "DPM") "wave_cor" else "within_cor",
    lagged_effects = if (!is.null(Phi)) "Phi" else "lagged_effects"
  )
  if (!is.null(AF_proportion)) {
    intraclass_correlation <- AF_proportion
  }
  if (!is.null(ICC)) {
    intraclass_correlation <- ICC
  }
  if (!is.null(AF_cor)) {
    RI_cor <- AF_cor
  }
  if (!is.null(wave_cor)) {
    within_cor <- wave_cor
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
  icheck_T_model(time_points, estimate_ME, model)
  icheck_ICC(ICC, arg = argument_names$intraclass_correlation)
  icheck_cor(RI_cor, arg = argument_names$RI_cor)
  icheck_cor(within_cor, arg = argument_names$within_cor)
  icheck_lagged_effects(lagged_effects, arg = argument_names$lagged_effects)
  icheck_moment(skewness)
  icheck_moment(kurtosis)
  icheck_significance_criterion(significance_criterion)
  icheck_ME(estimate_ME)
  icheck_reps(reps)
  seed <- icheck_seed(seed)
  icheck_constraints(constraints, estimate_ME)
  icheck_bounds(bounds, constraints, software)
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
  icheck_DPM_compatibility(model, reliability, estimate_ME, software, constraints, bounds)
  icheck_loadings(loadings, time_points, software, constraints, model = model)
  icheck_constraints_software(constraints, software)

  # powRICLPM updates
  if (!is.null(bootstrap_reps)) {
    cli::cli_alert_warning("The argument {.arg bootstrap_reps} is superseded. Uncertainty regarding simulation estimates is now computed analytically based on Morris et al. (2017).")
  }

  # Compute population parameter values for data generation
  if (model == "RICLPM") {
    Psi <- compute_Psi(lagged_effects, within_cor)
    icheck_Psi(Psi)
  } else {
    Psi <- NULL
    lapply(time_points, function(time_points) {
      loadings_i <- inormalize_loadings(loadings, time_points, model = model)
      dpm_values <- compute_DPM_values(
        lagged_effects = lagged_effects,
        wave_cor = within_cor,
        AF_var = compute_AF_var(ICC),
        AF_cov = compute_AF_cov(RI_cor, compute_AF_var(ICC)),
        loadings = loadings_i
      )
      lapply(seq_len(dim(dpm_values$Psi)[3]), function(i) {
        icheck_Psi(dpm_values$Psi[, , i])
      })
    })
  }

  # Get candidate sample sizes
  if (is.null(sample_size)) {
    icheck_sample_size_search(search_lower, search_upper, search_step)
    sample_size <- seq(search_lower, search_upper, search_step)
  }
  icheck_N(sample_size, time_points, constraints, estimate_ME, model = model)

  # Inform user that input check is complete
  cli::cli_alert_success("Argument checking complete.")

  # Setup analysis
  conditions <- create_conditions(
    model = model,
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
    save_path,
    software
  )

  if (software == "lavaan") {
    # Inform user that simulations in lavaan have started
    cli::cli_h2("\nPerforming Simulations Using lavaan")

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
        save_path = save_path,
        time_taken = time_taken,
        version = utils::packageVersion("powRICLPM"),
        call = call_powRICLPM,
        argument_names = argument_names
      )
    )

    class(out) <- c("powRICLPM", class(out))

    return(out)
  } else if (software == "Mplus") {
    # Inform user that simulations in Mplus have started
    cli::cli_h2("\nPerforming Simulations Using Mplus")

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

