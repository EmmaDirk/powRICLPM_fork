#' Extract Information From \code{powRICLPM} Object
#'
#' Extract information stored within a \code{powRICLPM} object (internally used by \code{\link{print.powRICLPM}} and \code{\link{summary.powRICLPM}}). See "Details" for which pieces of information can be extracted. The information is presented by condition (i.e., sample size, number of time points, intraclass correlation or accumulating-factor proportion, and reliability when applicable).
#'
#' @param from A \code{powRICLPM} object
#' @param what A character string, denoting the information to extract, such as "conditions", "estimation_problems", "results", or "names" (see "Details").
#' @param parameter (optional) When \code{what = "results"}, a character string denoting the parameter to extract the results for.
#'
#' @details
#' The following information can be extracted from the \code{powRICLPM} object:
#'
#' \itemize{
#'   \item \code{conditions}: A \code{data.frame} with the different experimental conditions per row, where each condition is defined by a unique combination of sample size, number of time points, intraclass correlation or DPM accumulating-factor proportion, and reliability when applicable. Matrix reliability specifications are shown with a compact variable-specific label.
#'   \item \code{sample_size}, \code{time_points}, \code{intraclass_correlation}, \code{ICC}, \code{AF_proportion}, or \code{reliability}: The same conditions \code{data.frame}. \code{reliability} is available for DPM objects when measurement error is part of the DPM conditions.
#'   \item \code{estimation_problems}: The proportion of fatal errors, inadmissible values, or non-converged estimations (columns) per experimental conditions (row).
#'   \item \code{results}: The average estimate (\code{average}), signed difference between the average estimate and population value (\code{bias}), minimum estimate (\code{minimum}), empirical standard error of parameter estimates (\code{EmpSE}), the average standard error (\code{SEAvg}), the mean square error (\code{MSE}), the average width of the confidence interval (\code{accuracy}), the coverage rate (\code{coverage}), and the proportion of times the \emph{p}-value was lower than the significance criterion (\code{power}). It requires setting the \code{parameter = "..."} argument.
#'   \item \code{names}: Parameter names available in every experimental condition.
#'   \item \code{uncertainty}: Monte Carlo standard errors for a specific parameter. It requires setting the \code{parameter = "..."} argument.
#' }
#'
#' @return A \code{data.frame}.
#' @export
#'
#' @examples
#' \dontshow{
#' load(system.file("extdata", "out_preliminary.RData", package = "powRICLPM"))
#' }
#' # Return data frame with number of estimation problems per experimental condition
#' give(out_preliminary, "estimation_problems")
#'
#' # Return data frame with performance measures for "wB2~wA1" per experimental condition
#' give(out_preliminary, "results", parameter = "wB2~wA1")
#'
#' # Return character vector with parameter names
#' give(out_preliminary, "names")
give <- function(from, what, parameter = NULL) {

  # Input checking
  icheck_object_summary(from)
  raw_what <- what
  icc_column_label <- if (identical(what, "ICC")) {
    "ICC"
  } else if (identical(what, "AF_proportion")) {
    "AF_proportion"
  } else if (identical(what, "intraclass_correlation")) {
    "intraclass_correlation"
  } else {
    iicc_value_name(object = from)
  }
  what <- normalize_intraclass_correlation_value(what)
  icheck_what_give(what)
  if (iis_DPM_object(from) &&
      any(identical(raw_what, "intraclass_correlation"), identical(raw_what, "ICC"))) {
    cli::cli_abort(
      c(
        "`intraclass_correlation` and `ICC` are not available for DPM objects:",
        i = "Use `AF_proportion` for the DPM accumulating-factor proportion."
      )
    )
  }
  if (!iis_DPM_object(from) && identical(raw_what, "AF_proportion")) {
    cli::cli_abort(
      c(
        "`AF_proportion` is only available for DPM objects:",
        i = "Use `intraclass_correlation` or `ICC` for RI-CLPM and STARTS objects."
      )
    )
  }
  if (iis_DPM_object(from) && identical(what, "reliability") &&
      !iDPM_has_reliability_conditions(from)) {
    cli::cli_abort(
      c(
        "`reliability` is not available for DPM objects:",
        i = "This DPM object does not include measurement error, so reliability is not part of its simulation conditions.",
        i = "Use `give(object, 'conditions')` or `give(object, 'AF_proportion')` to inspect DPM conditions."
      )
    )
  }

  # Call to relevant give_*() based on `what` argument
  if (what == "conditions" || what == "sample_size" || what == "time_points" ||
      what == "intraclass_correlation" || what == "reliability") {
    out <- give_powRICLPM_conditions(object = from)
  } else if (what == "estimation_problems") {
    out <- give_powRICLPM_estimation_problems(object = from)
  } else if (what == "results") {
      out <- give_powRICLPM_results(object = from, parameter = parameter)
  } else if (what == "names") {
      out <- give_powRICLPM_parameter_names(object = from)
  } else if (what == "uncertainty") {
      out <- give_powRICLPM_MCSE_parameter(object = from, parameter = parameter)
  }

  out <- idrop_DPM_reliability_column(from, out)
  ilabel_icc_column(out, icc_column_label)
}

give_powRICLPM_conditions <- function(object) {

  # Combine sample sizes and simulated power across conditions
  d <- do.call(rbind, lapply(object$conditions, function(condition) {
    data.frame(
      sample_size = condition$sample_size,
      time_points = condition$time_points,
      ICC = icondition_proportion(condition),
      reliability = condition$reliability,
      stringsAsFactors = FALSE
    )
  }))
  return(d)
}

iis_DPM_object <- function(object) {
  identical(object$session$model, "DPM")
}

iDPM_has_reliability_conditions <- function(object) {
  if (!iis_DPM_object(object)) {
    return(TRUE)
  }
  isTRUE(object$session$estimate_ME) ||
    any(vapply(object$conditions, function(condition) {
      !is.null(condition$ME_var) && any(condition$ME_var != 0)
    }, logical(1)))
}

imodel_display_name <- function(object) {
  if (iis_DPM_object(object)) {
    return("Dynamic Panel Model (DPM)")
  }
  if (isTRUE(object$session$estimate_ME)) {
    return("STARTS model")
  }
  "RI-CLPM"
}

imodel_display_short <- function(object) {
  if (iis_DPM_object(object)) {
    return("DPM")
  }
  if (isTRUE(object$session$estimate_ME)) {
    return("STARTS")
  }
  "RI-CLPM"
}

ipowRICLPM_version <- function(object) {
  if (!is.null(object$session$powRICLPM_version)) {
    return(object$session$powRICLPM_version)
  }
  if (!is.null(object$session$version)) {
    return(object$session$version)
  }
  utils::packageVersion("powRICLPM")
}

iexperimental_condition_label <- function(n) {
  if (identical(n, 1L)) {
    return("experimental condition")
  }
  "experimental conditions"
}

idrop_DPM_reliability_column <- function(object, x) {
  if (!iis_DPM_object(object) || !is.data.frame(x) || !"reliability" %in% names(x)) {
    return(x)
  }
  if (iDPM_has_reliability_conditions(object)) {
    return(x)
  }
  x$reliability <- NULL
  x
}

ihas_free_loadings <- function(condition) {
  if (is.null(condition$constraints)) {
    return(FALSE)
  }
  has_constraint(condition$constraints, "loadings_free")
}

inote_condition_loading_anchor <- function(object, output = NULL) {
  if (!inherits(object, "powRICLPM")) {
    return(invisible(NULL))
  }
  if (is.character(output)) {
    return(invisible(NULL))
  }
  uses_free_loadings <- any(vapply(object$conditions, ihas_free_loadings, logical(1)))
  if (!uses_free_loadings) {
    return(invisible(NULL))
  }
  if (identical(object$session$model, "DPM")) {
    cli::cli_alert_info(
      "With free DPM loadings, `AF_proportion` refers to wave 2, where the accumulating-factor loading is fixed to 1."
    )
  } else {
    cli::cli_alert_info(
      "With free random-intercept loadings, `intraclass_correlation`/`ICC` is anchored at wave 1, where the random-intercept loading is fixed to 1."
    )
  }
  invisible(NULL)
}

give_powRICLPM_estimation_problems <- function(object) {

  # Combine sample sizes and simulated power across conditions
  d <- do.call(rbind, lapply(object$conditions, function(condition) {
    data.frame(
      sample_size = condition$sample_size,
      time_points = condition$time_points,
      ICC = icondition_proportion(condition),
      reliability = condition$reliability,
      errors = condition$estimation_information$n_error,
      not_converged = condition$estimation_information$n_nonconvergence,
      inadmissible = condition$estimation_information$n_inadmissible,
      stringsAsFactors = FALSE
    )
  }))
  return(d)
}

give_powRICLPM_results <- function(object, parameter = NULL) {
  if (is.null(parameter)) {
    cli::cli_abort(
      c(
        "No {.arg parameter} was specified:",
        i = "{.code give(object, 'results')} needs to know which parameter to extract.",
        i = "Use {.code give(object, 'names')} to see available parameter names."
      )
    )
  }

  # Combine simulation results across experimental conditions
  d <- do.call(rbind, lapply(object$conditions, function(condition) {

    # Extract and round estimates per condition
    estimates <- if (is.data.frame(condition$estimates)) {
      condition$estimates[condition$estimates$parameter == parameter, -1, drop = FALSE]
    } else {
      data.frame()
    }
    if (nrow(estimates) == 0L) {
      if (parameter %in% give_powRICLPM_parameter_names_any(object)) {
        cli::cli_abort(
          c(
            "No usable estimates are available for {.arg parameter} `{parameter}` in at least one experimental condition:",
            i = "All replications for that condition failed, did not converge, or were inadmissible.",
            i = "Use {.code give(object, 'estimation_problems')} to inspect estimation problems."
          )
        )
      }
      cli::cli_abort(
        c(
          "The requested {.arg parameter} was not found:",
          x = paste0("No results are available for parameter `", parameter, "`."),
          i = "Use {.code give(object, 'names')} to see available parameter names."
        )
      )
    }
    estimates <- round(estimates, digits = 3)

    # Combine extracted info in data frame
    data.frame(
      sample_size = condition$sample_size,
      time_points = condition$time_points,
      ICC = icondition_proportion(condition),
      reliability = condition$reliability,
      estimates
    )
  }))
  return(d)
}

give_powRICLPM_parameter_names <- function(object) {
  Reduce(intersect, lapply(object$conditions, function(condition) {
    icondition_parameter_names(condition)
  }))
}

give_powRICLPM_parameter_names_any <- function(object) {
  unique(unlist(lapply(object$conditions, function(condition) {
    icondition_parameter_names(condition)
  }), use.names = FALSE))
}

icondition_parameter_names <- function(condition) {
  if (is.data.frame(condition$estimates) && "parameter" %in% names(condition$estimates)) {
    return(condition$estimates$parameter)
  }
  if (is.data.frame(condition$est_tab) && all(c("lhs", "op", "rhs", "free") %in% names(condition$est_tab))) {
    return(paste0(
      condition$est_tab$lhs,
      condition$est_tab$op,
      condition$est_tab$rhs
    )[condition$est_tab$free])
  }
  character()
}


give_powRICLPM_MCSE_parameter <- function(object, parameter) {

  # Combine data frames across conditions
  d <- do.call(rbind, lapply(object$conditions, function(condition) {

    # Extract rows from uncertainty where parameter matches
    uncertainty_filtered <- condition$MCSEs[which(condition$estimates$parameter == parameter), ]

    # Create data frame
    data.frame(
      sample_size = condition$sample_size,
      time_points = condition$time_points,
      ICC = icondition_proportion(condition),
      reliability = condition$reliability,
      uncertainty_filtered,
      stringsAsFactors = FALSE
    )
  }))
  return(d)
}

icondition_proportion <- function(condition) {
  if (!is.null(condition$AF_proportion)) {
    return(condition$AF_proportion)
  }
  condition$ICC
}
