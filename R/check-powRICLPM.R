#' Check \code{target} Argument
#'
#' \code{icheck_target()} tests if the \code{target} argument is an integer.
#'
#' @noRd
icheck_target <- function(x, arg = rlang::caller_arg(x), call = rlang::caller_env()) {
  if (!is.double(x)) {
    cli::cli_abort(
      c(
        "{.arg {arg}} must be a numeric:",
        x = paste0("Your {.arg {arg}} is a ", typeof(x), ".")
      ),
      call = call
    )
  }
  if (x >= 1 || x <= 0) {
    cli::cli_abort(
      c(
        "{.arg {arg}} must be between 0 and 1:",
        x = paste0("Your {.arg {arg}} is ", x, ".")
      )
    )
  }
}

#' Check \code{time_points} Argument
#'
#' \code{icheck_T()} checks if the \code{time_points} argument represents (a) valid number(s) of repeated measures.
#'
#' @noRd
icheck_T <- function(x, ME, arg = rlang::caller_arg(x), call = rlang::caller_env()) {
  if (!all(is.numeric(x))) {
    cli::cli_abort(
      c(
        "{.arg {arg}} must be a vector of integers:",
        x = "Not all elements are numeric."
      )
    )
  }
  if (!all(x %% 1 == 0)) {
    cli::cli_abort(
      c(
        "{.arg {arg}} must be a vector of integers:",
        x = "Not all elements are integers."
      )
    )
  }
  if (any(x < 3) && !ME) {
    cli::cli_abort(
      c(
        "Elements in {.arg {arg}} should be larger than 2:",
        i = "The RI-CLPM is not identified with fewer than 3 time points.",
        x = "You've supplied a number of time points smaller than 3."
      )
    )
  }
  if (any(x < 4) && ME) {
    cli::cli_abort(
      c(
        "Elements in {.arg {arg}} should be larger than 3:",
        i = "If you want to estimate measurement errors in the RI-CLPM, you should have at least 4 waves of data.",
        x = "You've supplied a number of time points smaller than 4."
      )
    )
  }
  if (any(x > 15)) {
    cli::cli_warn(
      c(
        "You've supplied a large number of time points:",
        i = "This can lead to computational problems in estimation. You might want to consider methods for intensive longitudinal data."
      )
    )
  }
}

#' Check \code{ICC} Argument
#'
#' \code{icheck_ICC()} checks if the \code{ICC} argument represents (a) valid intraclass correlation(s).
#'
#' @noRd
icheck_ICC <- function(x, arg = rlang::caller_arg(x), call = rlang::caller_env()) {
  if (!all(is.double(x))) {
    cli::cli_abort(
      c(
        "{.arg {arg}} should be a double (vector):",
        x = "Not all elements in {.arg {arg}} are a double."
      )
    )
  }
  if (!all(x > 0 & x < 1)) {
    cli::cli_abort(
      c(
        "Elements in {.arg {arg}} must be between 0 and 1:",
        x = "Some elements in {.arg {arg}} are smaller than 0 or larger than 1."
      )
    )
  }
  if (!all(x < .99)) {
    cli::cli_abort(
      c(
        "Elements in {.arg {arg}} are very close to 1:",
        i = "This can lead to problems with estimation. Use a maximum of .99."
      )
    )
  }
}

#' Check Correlation Arguments
#'
#' \code{icheck_cor()} tests if an argument represents a valid correlation.
#'
#' @noRd
icheck_cor <- function(x, arg = rlang::caller_arg(x), call = rlang::caller_env()) {
  if (!is.double(x)) {
    cli::cli_abort(
      c(
        "{.arg {arg}} must be a double:",
        x = paste0("You've provided a ", typeof(x), ".")
      )
    )
  }
  if (x < -1 || x > 1) {
    cli::cli_abort(
      c(
        "A correlation must be between -1 and 1:",
        x = paste0("Your {.arg {arg}} was ", x, ".")
      )
    )
  }
  if (length(x) > 1) {
    cli::cli_abort(
      c(
        "{.arg {arg}} should a double of length 1:",
        x = "You've provided multiple values for {.arg {arg}}."
      )
    )
  }
}

#' Check \code{lagged_effects} Argument
#'
#' \code{icheck_lagged_effects()} tests if \code{lagged_effects} represents a valid regression matrix for the within-components of the RI-CLPM.
#'
#' @noRd
icheck_lagged_effects <- function(x, arg = rlang::caller_arg(x), call = rlang::caller_env()) {
  if (!is.matrix(x)) {
    cli::cli_abort(
      c(
        "{.arg {arg}} must be a matrix:",
        x = paste0("Your {.arg {arg}} is a ", typeof(x), ".")
      )
    )
  }
  if (!is_unit(x)) {
    cli::cli_abort(
      c(
        "{.arg {arg}} must specify a stationary process:",
        i = "This is checked by testing if the eigenvalues of {.arg {arg}} lie within unit circle.",
        x = "Eigenvalues of {.arg {arg}} are not within unit circle. Try out smaller lagged effects?"
      )
    )
  }
}

#' Check \code{reliability} Argument
#'
#' \code{icheck_rel()} checks if the \code{reliability} argument is a valid reliability coefficient.
#'
#' @noRd
icheck_rel <- function(x, arg = rlang::caller_arg(x), call = rlang::caller_env()) {
  if (!all(is.numeric(x))) {
    cli::cli_abort(
      c(
        "{.arg {arg}} must be a numeric (vector):",
        x = paste0("Your {.arg {arg}} is a ", class(x), ".")
      )
    )
  }
  if (!all(x <= 1 & x > 0)) {
    cli::cli_abort(
      c(
        "Elements in {.arg {arg}} must be between 0 and 1:",
        x = paste0("Some elements in {.arg {arg}} are smaller than 0 or larger than 1.")
      )
    )
  }
  if (!all(x > 0.1)) {
    cli::cli_abort(
      c(
        "Some elements in {.arg {arg}} are close to 0:",
        x = paste0("A low reliability can lead to problems with estimation. Use a minimum of 0.1.")
      )
    )
  }
}

#' Check \code{loadings} Argument
#'
#' \code{icheck_loadings()} checks if the \code{loadings} argument is a valid
#' random-intercept loading specification for the lavaan data-generating model.
#'
#' @noRd
icheck_loadings <- function(x, time_points, software, constraints = "none",
                            arg = rlang::caller_arg(x),
                            t_arg = rlang::caller_arg(time_points),
                            con_arg = rlang::caller_arg(constraints),
                            call = rlang::caller_env()) {
  if (is.null(x)) {
    return(invisible(NULL))
  }
  arg <- format_argument_label(arg, "loadings")
  t_arg <- format_argument_label(t_arg, "time_points")
  con_arg <- format_argument_label(con_arg, "constraints")
  if (software == "Mplus") {
    cli::cli_abort(
      c(
        "{.arg {arg}} can only be used with `software = 'lavaan'`:",
        x = "Time-varying random-intercept loadings are not available for Mplus."
      ),
      call = call
    )
  }
  if (length(time_points) != 1) {
    cli::cli_abort(
      c(
        "{.arg {arg}} can only be used with one value of {.arg {t_arg}}:",
        i = "If you want to compare power for multiple values of {.arg {t_arg}} while using free loadings, use multiple {.fun powRICLPM} calls.",
        x = paste0("length({.arg {t_arg}}) = ", length(time_points), ".")
      ),
      call = call
    )
  }
  if (!is.numeric(x)) {
    cli::cli_abort(
      c(
        "{.arg {arg}} must be a numeric vector or matrix:",
        x = paste0("Your {.arg {arg}} is a ", typeof(x), ".")
      ),
      call = call
    )
  }
  if (is.matrix(x)) {
    if (!identical(dim(x), c(2L, as.integer(time_points)))) {
      cli::cli_abort(
        c(
          "{.arg {arg}} must have 2 rows and one column per time point:",
          i = "Rows correspond to variables A and B; columns correspond to time points.",
          x = paste0("Your {.arg {arg}} has dimensions ", paste(dim(x), collapse = " x "),
                     ", implying ", ncol(x), " time points, but {.arg {t_arg}} = ",
                     time_points, ".")
        ),
        call = call
      )
    }
    first_loadings <- x[, 1]
  } else {
    if (length(x) != time_points) {
      cli::cli_abort(
        c(
          "{.arg {arg}} must have one value per time point:",
          x = paste0("Your {.arg {arg}} has length ", length(x),
                     ", but {.arg {t_arg}} = ", time_points, ".")
        ),
        call = call
      )
    }
    first_loadings <- x[1]
  }
  if (!all(is.finite(x))) {
    cli::cli_abort(
      c(
        "{.arg {arg}} must contain only finite values:",
        x = paste0("Your {.arg {arg}} contains: ", format_loadings(x), ".")
      ),
      call = call
    )
  }
  if (!all(first_loadings == 1)) {
    cli::cli_abort(
      c(
        "The first random-intercept loading must be fixed to 1:",
        i = "{.arg {arg}} is specified relative to the first occasion.",
        x = "Set the first loading = 1 for each variable."
      ),
      call = call
    )
  }
  if (!has_constraint(constraints, "RI_loadings_free")) {
    cli::cli_abort(
      c(
        "{.arg {arg}} can only be used when random-intercept loadings are freely estimated:",
        i = "Use `constraints = 'RI_loadings_free'` or include `'RI_loadings_free'` in a constraint vector.",
        x = paste0("You supplied {.arg {arg}}, but {.arg {con_arg}} = ", format_constraints(constraints), ".")
      ),
      call = call
    )
  }
  invisible(NULL)
}

inormalize_loadings <- function(loadings, time_points) {
  if (is.null(loadings)) {
    return(matrix(1, nrow = 2, ncol = time_points))
  }
  if (is.matrix(loadings)) {
    return(loadings)
  }
  rbind(loadings, loadings)
}

#' Check \code{moment} Arguments
#'
#' \code{icheck_moment()} tests if a \code{moment} argument (e.g., skewness, kurtosis) is a numeric value of length 1.
#'
#' @noRd
icheck_moment <- function(x, arg = rlang::caller_arg(x), call = rlang::caller_env()) {
  if (!is.numeric(x)) {
    cli::cli_abort(
      c(
        "{.arg {arg}} must be a numeric (vector):",
        x = paste0("Your {.arg {arg}} is a ", typeof(x), ".")
      )
    )
  }
  if (length(x) != 1) {
    cli::cli_abort(
      c(
        "{.arg {arg}} should be of length 1:",
        i = "This version of `powRICLPM` only accepts a single value for this factor.",
        x = paste0("Your {.arg {arg}} is of length ", length(x), ".")
      )
    )
  }
}

#' Check \code{significance_criterion} Argument
#'
#' \code{icheck_significance_criterion()} tests if the \code{alpha} argument is numeric value representing a valid significance criterion.
#'
#' @param x A double.
#'
#' @noRd
icheck_significance_criterion <- function(x, arg = rlang::caller_arg(x), call = rlang::caller_env()) {
  if (!is.double(x)) {
    cli::cli_abort(
      c(
        "{.arg {arg}} must be a double:",
        x = paste0("Your {.arg {arg}} is a ", typeof(x), ".")
      )
    )
  }
  if (1 < x || x < 0) {
    cli::cli_abort(
      c(
        "{.arg {arg}} must be between 0 and 1:",
        x = paste0("Your {.arg {arg}} is ", x, ".")

      )
    )
  }
}

#' Check \code{estimate_ME} Argument
#'
#' \code{icheck_ME()} checks if \code{estimate_ME} is a logical of length 1.
#'
#' @noRd
icheck_ME <- function(x, arg = rlang::caller_arg(x), call = rlang::caller_env()) {
  if (!is.logical(x)) {
    cli::cli_abort(
      c(
        "{.arg {arg}} must be a logical:",
        x = paste0("Your {.arg {arg}} is a ", class(x), ".")
      )
    )
  }
  if (length(x) > 1) {
    cli::cli_abort(
      c(
        "{.arg {arg}} should be of length 1:",
        x = paste0("Your {.arg {arg}} is of length ", length(x), ".")
      )
    )
  }
}

#' Check \code{reps} Argument
#'
#' \code{icheck_reps()} tests if provided argument is a valid number of replications.
#'
#' @noRd
icheck_reps <- function(x, arg = rlang::caller_arg(x), call = rlang::caller_env()) {
  if (!is.numeric(x)) {
    cli::cli_abort(
      c(
        "{.arg {arg}} must be an integer:",
        x = paste0("Your {.arg {arg}} is a ", typeof(x), ".")
      )
    )
  }
  if (x %% 1 != 0) {
    cli::cli_abort(
      c(
        "{.arg {arg}} must be an integer:",
        x = paste0("Your {.arg {arg}} is not a 'whole' number.")
      )
    )
  }
  if (x < 0) {
    cli::cli_abort(
      c(
        "{.arg {arg}} must be a positive integer:",
        x = paste0("Your {.arg {arg}} is: ", x, ".")
      )
    )
  }
}


#' Check \code{seed} Argument
#'
#' \code{icheck_seed()} checks if a seed is specified. If yes, then it tests provided argument is a valid seed. If not, a randomly generates is returned.
#'
#' @noRd
icheck_seed <- function(x, arg = rlang::caller_arg(x), call = rlang::caller_env()) {
  if (is.na(x)) {
    x <- floor(stats::runif(1, 0, 1000000))
    cli::cli_warn(
      c("No seed was specified. A random seed was generated: {.value {x}}")
    )
    return(x)
  }
  if (!is.numeric(x)) {
    cli::cli_abort(
      c(
        "{.arg {arg}} must be a numeric:",
        x = paste0("Your {.arg {arg}} is a {typeof(x)}.")
      )
    )
  } else if (x %% 1 != 0) {
    cli::cli_abort(
      c(
        "{.arg {arg}} should be an integer:",
        x = "Your {.arg {arg}} is not a 'whole' number."
      )
    )
  }
  return(x)
}

#' Check \code{constraints} Argument
#'
#' \code{check_constraints()} checks if the \code{constraints} arguments refers to a valid set of constraints.
#'
#' @noRd
icheck_constraints <- function(x, ME, arg = rlang::caller_arg(x), call = rlang::caller_env()) {
  if (!is.character(x)) {
    cli::cli_abort(
      c(
        "{.arg {arg}} must be a character vector:",
        x = paste0("Your {.arg {arg}} is ", typeof(x), ".")
      )
    )
  }
  if (length(x) == 0) {
    cli::cli_abort(
      c(
        "{.arg {arg}} must contain at least one constraint option:",
        i = "Use 'none' when no constraints should be imposed."
      )
    )
  }

  valid_constraints <- c(
    "none", "lagged", "residuals", "within", "stationarity", "ME",
    "RI_loadings_free"
  )
  if (!all(x %in% valid_constraints)) {
    cli::cli_abort(
      c(
        "{.arg {arg}} contains invalid constraints:",
        i = paste0("Valid constraints are: ", paste(valid_constraints, collapse = ", "), "."),
        x = paste0("Your {.arg {arg}} is ", format_constraints(x), ".")
      )
    )
  }
  constraints <- normalize_constraints(x)
  if (has_constraint(constraints, "none") && length(constraints) > 1) {
    cli::cli_abort(
      c(
        "{.arg {arg}} can only include 'none' by itself:",
        x = paste0("Your {.arg {arg}} is ", format_constraints(constraints), ".")
      )
    )
  }
  if ("within" %in% constraints && any(constraints %in% c("lagged", "residuals"))) {
    cli::cli_abort(
      c(
        "{.arg {arg}} cannot combine 'within' with 'lagged' or 'residuals':",
        i = "'within' already constrains both lagged effects and residual variances.",
        x = paste0("Your {.arg {arg}} is ", format_constraints(constraints), ".")
      )
    )
  }
  if (has_constraint(constraints, "stationarity") &&
      any(constraints %in% c("lagged", "residuals", "within"))) {
    cli::cli_abort(
      c(
        "{.arg {arg}} cannot combine 'stationarity' with lagged or residual constraints:",
        i = "'stationarity' already imposes its own lagged-effect and residual-variance constraints.",
        x = paste0("Your {.arg {arg}} is ", format_constraints(constraints), ".")
      )
    )
  }
  if (has_constraint(constraints, "ME") && !ME) {
    cli::cli_abort(
      c(
        "{.arg {arg}} can only be set to 'ME' when `estimate_ME = TRUE`:",
        x = "Your {.arg {arg}} is FALSE."
      )
    )
  }
}

#' Check \code{estimator} Argument
#'
#' \code{icheck_estimator()} checks if the \code{estimator} argument refers to a valid maximum likelihood estimators implemented in \pkg{lavaan}.
#'
#' @noRd
icheck_estimator <- function(x, skewness, kurtosis, arg = rlang::caller_arg(x), call = rlang::caller_env()) {
  if (is.na(x)) {
    if (skewness != 0 || kurtosis != 0) {
      cli::cli_alert(
        c(
          i = "Estimator defaulting to `MLR` for skewed and/or kurtosed data."
        )
      )
      return("MLR")
    } else {
      return("ML")
    }
  }
  if (!any(x == c("ML", "MLR", "MLM", "MLMVS", "MLMV", "MLF"))) {
    cli::cli_abort(
      c(
        "{.arg {arg}} must be 'ML', 'MLR', 'MLM', 'MLMVS', 'MLMV', or 'MLF':",
        x = paste0("Your {.arg {arg}} is '", x, "'.")
      )
    )
  }
}

#' Check \code{save_path} Argument
#'
#' \code{icheck_path()} tests if the provided argument is a valid path argument that points to an existing folder.
#'
#' @noRd
icheck_path <- function(x, software, arg = rlang::caller_arg(x), call = rlang::caller_env()) {
  if (!is.null(x)) {
    if (!is.character(x)) {
      cli::cli_abort(
        c(
          "{.arg {arg}} must be a character string:",
          x = paste0("Your {.arg {arg}} is a ", typeof(x), ".")
        )
      )
    }
    if (!dir.exists(x)) {
      cli::cli_abort(
        c(
          "{.arg {arg}} must point to an existing folder:",
          x = "The directory in {.arg {arg}} does not exist. Create this folder or change {.arg {arg}}."
        )
      )
    }
  } else if (is.null(x) && software == "Mplus") {
    path <- getwd()
    return(path)
  }
  return(x)
}

#' Check computed \code{Psi} matrix
#'
#' \code{icheck_Psi()} checks if \code{Psi} represents a valid variance-covariance matrix.
#'
#' @noRd
icheck_Psi <- function(x, arg = rlang::caller_arg(x), call = rlang::caller_env()) {
  if (!is_PD(x)) {
    cli::cli_abort(
      c(
        "The residual variance-covariance matrix for within-components (Psi) must be positive definite:",
        i = "It is computed from the specified `lagged_effects` and `within_cor` arguments.",
        x = "Psi is not positive definite. Try smaller values for `lagged_effects` and `within_cor`?"
      )
    )
  }
}

#' Check \code{sample_size} Argument
#'
#' \code{icheck_N()} checks if the \code{sample_size} argument represents valid sample size(s) given the specified constraints, number of time points, and estimation of measurement error.
#'
#' @noRd
icheck_N <- function(x, t, constraints = "none", ME = FALSE, arg = rlang::caller_arg(x), call = rlang::caller_env()) {
  if (!all(is.numeric(x))) {
    cli::cli_abort(
      c(
        "{.arg {arg}} must be an integer vector:",
        x = "Not all elements are numeric."
      )
    )
  }
  if (!all(x %% 1 == 0)) {
    cli::cli_abort(
      c(
        "{.arg {arg}} must be an integer vector:",
        x = "Not all elements are integers."
      )
    )
  }
  if (!all(x > 0)) {
    cli::cli_abort(
      c(
        "{.arg {arg}} must only contain positive integers:",
        x = "Your {.arg {arg}} contains negative numbers."
      )
    )
  }
  n_parameters <- count_parameters(2, t, constraints, ME)
  if (!all(x > n_parameters)) {
    cli::cli_abort(
      c(
        "The sample size must be larger than the number of parameters estimated (for all experimental conditions):",
        i = paste0("The highest number of estimated parameters in an experimental condition is ", n_parameters, "."),
        x = paste0("The smallest sample size you have specified is ", min(x), ".")
      )
    )
  }
}


#' Check sample size search arguments
#'
#' \code{icheck_sample_size_search()} checks if the \code{search_*}
#' arguments define a valid range when \code{sample_size} is not supplied.
#'
#' @noRd
icheck_sample_size_search <- function(search_lower, search_upper, search_step, call = rlang::caller_env()) {
  search_args <- list(
    search_lower = search_lower,
    search_upper = search_upper,
    search_step = search_step
  )
  missing_args <- names(search_args)[vapply(search_args, is.null, logical(1))]

  if (length(missing_args) > 0) {
    cli::cli_abort(
      c(
        "Either {.arg sample_size} or a complete sample-size search range must be supplied:",
        i = "Use {.arg sample_size} to evaluate specific sample sizes.",
        i = "Or use {.arg search_lower}, {.arg search_upper}, and {.arg search_step} to generate a range.",
        x = paste0("Missing search argument", if (length(missing_args) > 1) "s" else "", ": ", paste(missing_args, collapse = ", "), ".")
      ),
      call = call
    )
  }

  invalid_type <- names(search_args)[!vapply(search_args, function(x) {
    is.numeric(x) && length(x) == 1 && !is.na(x)
  }, logical(1))]

  if (length(invalid_type) > 0) {
    cli::cli_abort(
      c(
        "The sample-size search range must use single numeric values:",
        x = paste0("Check: ", paste(invalid_type, collapse = ", "), ".")
      ),
      call = call
    )
  }

  invalid_integer <- names(search_args)[!vapply(search_args, function(x) x %% 1 == 0, logical(1))]

  if (length(invalid_integer) > 0) {
    cli::cli_abort(
      c(
        "The sample-size search range must use integer values:",
        x = paste0("Check: ", paste(invalid_integer, collapse = ", "), ".")
      ),
      call = call
    )
  }

  invalid_positive <- names(search_args)[!vapply(search_args, function(x) x > 0, logical(1))]

  if (length(invalid_positive) > 0) {
    cli::cli_abort(
      c(
        "The sample-size search range must use positive integers:",
        x = paste0("Check: ", paste(invalid_positive, collapse = ", "), ".")
      ),
      call = call
    )
  }

  if (search_upper < search_lower) {
    cli::cli_abort(
      c(
        "{.arg search_upper} must be greater than or equal to {.arg search_lower}:",
        x = paste0("You supplied search_lower = ", search_lower, " and search_upper = ", search_upper, ".")
      ),
      call = call
    )
  }
}

#' Check \code{bounds} Argument
#'
#' \code{check_bounds()} tests if \code{bounds} is a logical, and whether bounded estimation can be used.
#'
#' @inheritParams powRICLPM
#'
#' @noRd
icheck_bounds <- function(bounds, constraints, software,
                          con = rlang::caller_arg(constraints), soft = rlang::caller_arg(software)) {
  if (!is.logical(bounds)) {
    cli::cli_abort(
      c(
        "`bounds` should be a `logical`:",
        x = paste0("Your `bounds` is a `", typeof(bounds), "`.")
      )
    )
  }
  constraints_for_bounds <- setdiff(
    normalize_constraints(constraints),
    c("none", "RI_loadings_free")
  )
  if (bounds && length(constraints_for_bounds) > 0) {
    cli::cli_abort(
      c(
        "Bounded estimation can only be used without equality or time-invariance constraints on the estimation model:",
        x = paste0("You've placed the following constraints on the estimation model: ", format_constraints(constraints), ".")
      )
    )
  }
  if (bounds && software == "Mplus") {
    cli::cli_abort(
      c(
        "Bounded estimation can only be used with lavaan:",
        x = "You've set the `software` argument to: {.arg {soft}}"
      )
    )
  }
}

icheck_constraints_software <- function(constraints, software) {
  constraints <- normalize_constraints(constraints)
  if (software == "Mplus" && has_constraint(constraints, "RI_loadings_free")) {
    cli::cli_abort(
      c(
        "Free random-intercept loadings are not available for Mplus:",
        x = paste0(
          "`constraints = ", format_constraints(constraints),
          "` can only be used with `software = 'lavaan'`."
        )
      )
    )
  }
  if (software == "Mplus" &&
      length(constraints) > 1 &&
      !has_constraint(constraints, "ME") &&
      !is_lagged_residuals_constraint(constraints)) {
    cli::cli_abort(
      c(
        "This vector-valued constraint combination is not supported for Mplus:",
        i = "`constraints = c('lagged', 'residuals')` is supported for Mplus, where it is treated as `constraints = 'within'`. The `ME` constraint can also be combined with other Mplus-supported constraints.",
        x = paste0("Your `constraints` argument is ", format_constraints(constraints), ".")
      )
    )
  }
}


normalize_constraints <- function(x) {
  if (is.list(x) && length(x) == 1) {
    x <- x[[1]]
  }
  unique(x)
}


is_lagged_residuals_constraint <- function(constraints) {
  constraints <- normalize_constraints(constraints)
  length(constraints) == 2 && setequal(constraints, c("lagged", "residuals"))
}


normalize_constraints_for_software <- function(constraints, software) {
  constraints <- normalize_constraints(constraints)
  if (software == "Mplus" && is_lagged_residuals_constraint(constraints)) {
    return("within")
  }
  constraints
}


has_constraint <- function(constraints, constraint) {
  constraints <- normalize_constraints(constraints)
  if (constraint == "lagged") {
    return(any(constraints %in% c("lagged", "within", "stationarity")))
  }
  if (constraint == "residuals") {
    return(any(constraints %in% c("residuals", "within")))
  }
  any(constraints == constraint)
}


format_constraints <- function(constraints) {
  constraints <- normalize_constraints(constraints)
  if (length(constraints) == 1) {
    return(constraints)
  }
  paste0("c(", paste0("'", constraints, "'", collapse = ", "), ")")
}


format_loadings <- function(loadings) {
  values <- paste(as.character(loadings), collapse = ", ")
  if (is.matrix(loadings)) {
    return(paste0(
      "matrix(c(", values, "), nrow = ", nrow(loadings),
      ", ncol = ", ncol(loadings), ")"
    ))
  }
  paste0("c(", values, ")")
}


format_argument_label <- function(arg, fallback) {
  if (grepl("^[.A-Za-z][.A-Za-z0-9_]*$", arg)) {
    return(arg)
  }
  fallback
}


#' Check `software` Argument
#'
#' `icheck_software()` tests if the `software` argument is a valid option, either "lavaan" or "Mplus".
#'
#' @noRd
icheck_software <- function(
    x, skewness, kurtosis,
    skew = rlang::caller_arg(skewness), kurt = rlang::caller_arg(kurtosis)
  ) {
  if (!any(x == c("lavaan", "Mplus"))) {
    cli::cli_abort(
      c(
        "The `software` argument does not specify a valid option:",
        x = "Your `software` is: {.value {x}}",
        i = "Please specify either `lavaan` or `Mplus`."
      )
    )
  }
  if (x == "Mplus" && !(skewness == 0 && kurtosis == 0)) {
    cli::cli_abort(
      c(
        "When using Mplus, it is not possible to generate skewed or kurtosed data:",
        x = "You've set `skewness` to {.arg {skew}} and `kurtosis` to {.arg {kurt}}."
      )
    )
  }
}


icheck_alpha <- function(x) {
  cli::cli_warn(c("`alpha` is deprecated and will be removed in a future version. Please use 'significance_criterion' instead."))
  return(x)
}


iabort_renamed_argument_conflict <- function(new, old, call = rlang::caller_env()) {
  cli::cli_abort(
    c(
      "Both {.arg {new}} and {.arg {old}} were supplied.",
      x = "Use only one of these arguments in the same call."
    ),
    call = call
  )
}


normalize_intraclass_correlation_value <- function(x) {
  if (identical(x, "ICC")) {
    return("intraclass_correlation")
  }
  x
}


iargument_display_name <- function(object = NULL, call = NULL, primary, alternate) {
  call_names <- character()
  if (!is.null(call)) {
    call_names <- names(as.list(call))
  }
  if (alternate %in% call_names && !(primary %in% call_names)) {
    return(alternate)
  }
  if (primary %in% call_names && !(alternate %in% call_names)) {
    return(primary)
  }

  session_argument <- NULL
  if (!is.null(object) && !is.null(object$session$argument_names[[primary]])) {
    session_argument <- object$session$argument_names[[primary]]
  }
  if (identical(session_argument, alternate)) {
    return(alternate)
  }

  if (!is.null(object) && !is.null(object$session$call)) {
    session_call_names <- names(as.list(object$session$call))
    if (alternate %in% session_call_names && !(primary %in% session_call_names)) {
      return(alternate)
    }
  }

  primary
}


iicc_value_name <- function(object = NULL, call = NULL) {
  iargument_display_name(
    object = object,
    call = call,
    primary = "intraclass_correlation",
    alternate = "ICC"
  )
}


iicc_table_name <- function(object = NULL, call = NULL) {
  if (identical(iicc_value_name(object = object, call = call), "ICC")) {
    return("ICC")
  }
  "Intraclass correlation"
}


ilabel_icc_column <- function(x, label) {
  if (is.data.frame(x) && "ICC" %in% names(x) && !identical(label, "ICC")) {
    names(x)[names(x) == "ICC"] <- label
  }
  x
}



