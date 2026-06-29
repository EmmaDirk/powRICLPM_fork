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
#'   \item \code{estimation_problems}: The number of fatal errors, inadmissible solutions, or non-converged estimations across replications for each experimental condition.
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
    iabort_condition_selector_unavailable("intraclass_correlation", from, "give")
  }
  if (!iis_DPM_object(from) && identical(raw_what, "AF_proportion")) {
    iabort_condition_selector_unavailable("AF_proportion", from, "give")
  }
  if (iis_DPM_object(from) && identical(what, "reliability") &&
      !iDPM_has_reliability_conditions(from)) {
    iabort_condition_selector_unavailable("reliability", from, "give")
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
    cbind(
      data.frame(
        sample_size = condition$sample_size,
        time_points = condition$time_points,
        ICC = icondition_proportion(condition),
        stringsAsFactors = FALSE
      ),
      icondition_reliability_columns(condition),
      data.frame(loadings = iloading_status(condition), stringsAsFactors = FALSE)
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
  reliability_cols <- intersect(c("reliability", "reliability_A", "reliability_B"), names(x))
  if (!iis_DPM_object(object) || !is.data.frame(x) || length(reliability_cols) == 0L) {
    return(x)
  }
  if (iDPM_has_reliability_conditions(object)) {
    return(x)
  }
  x[reliability_cols] <- NULL
  x
}

iabort_condition_selector_unavailable <- function(selector, object, context,
                                                  call = rlang::caller_env()) {
  context_label <- switch(
    context,
    summary = "summaries",
    plot = "plots",
    give = "objects",
    "objects"
  )
  if (identical(selector, "intraclass_correlation") || identical(selector, "ICC")) {
    cli::cli_abort(
      c(
        "`intraclass_correlation` and `ICC` are not available for DPM objects:",
        i = paste0("Use `AF_proportion` for DPM ", context_label, ".")
      ),
      call = call
    )
  }
  if (identical(selector, "AF_proportion")) {
    cli::cli_abort(
      c(
        "`AF_proportion` is only available for DPM objects:",
        i = paste0("Use `intraclass_correlation` or `ICC` for RI-CLPM and STARTS ", context_label, ".")
      ),
      call = call
    )
  }
  if (identical(selector, "reliability")) {
    hint <- switch(
      context,
      summary = "Select DPM conditions with `sample_size`, `time_points`, and `AF_proportion`.",
      plot = "Use `shape_by = 'none'` or map aesthetics to another simulation condition.",
      give = "Use `give(object, 'conditions')` or `give(object, 'AF_proportion')` to inspect DPM conditions.",
      "Use another simulation condition."
    )
    cli::cli_abort(
      c(
        "`reliability` is not available for DPM objects:",
        i = "This DPM object does not include measurement error, so reliability is not part of its simulation conditions.",
        i = hint
      ),
      call = call
    )
  }
}

icondition_reliability_columns <- function(condition) {
  reliability_matrix <- condition$reliability_matrix
  if (!is.null(reliability_matrix) &&
      !identical(reliability_matrix[1, ], reliability_matrix[2, ])) {
    return(data.frame(
      reliability_A = format_reliability_value(reliability_matrix[1, 1]),
      reliability_B = format_reliability_value(reliability_matrix[2, 1]),
      stringsAsFactors = FALSE
    ))
  }
  data.frame(reliability = condition$reliability, stringsAsFactors = FALSE)
}

ireliability_columns <- function(x) {
  intersect(c("reliability", "reliability_A", "reliability_B"), names(x))
}

icondition_key_columns <- function(x) {
  c("sample_size", "time_points", "ICC", ireliability_columns(x))
}

icondition_table_names <- function(object, x, icc_table_label) {
  col_names <- c("Sample size", "Time points", icc_table_label)
  if ("reliability" %in% names(x)) {
    col_names <- c(col_names, "Reliability")
  }
  if ("reliability_A" %in% names(x)) {
    col_names <- c(col_names, "Reliability A")
  }
  if ("reliability_B" %in% names(x)) {
    col_names <- c(col_names, "Reliability B")
  }
  if ("loadings" %in% names(x)) {
    col_names <- c(col_names, "Loadings")
  }
  col_names
}

iloading_status <- function(condition) {
  fitted <- if (ihas_free_loadings(condition)) "free" else "fixed"
  custom <- !is.null(condition$loadings) && any(condition$loadings != 1)
  paste0(fitted, if (custom) "*" else "")
}

ihas_custom_loadings <- function(object) {
  any(vapply(object$conditions, function(condition) {
    !is.null(condition$loadings) && any(condition$loadings != 1)
  }, logical(1)))
}

iprint_custom_loadings_note <- function(object) {
  if (!ihas_custom_loadings(object)) {
    return(invisible(NULL))
  }
  index <- if (length(object$conditions) == 1L) "1" else "i"
  cat(
    "\n* Custom data-generating loadings supplied; inspect with x$conditions[[",
    index,
    "]]$loadings.\n",
    sep = ""
  )
  invisible(NULL)
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
    cbind(
      data.frame(
        sample_size = condition$sample_size,
        time_points = condition$time_points,
        ICC = icondition_proportion(condition),
        stringsAsFactors = FALSE
      ),
      icondition_reliability_columns(condition),
      data.frame(
        loadings = iloading_status(condition),
        errors = condition$estimation_information$n_error,
        not_converged = condition$estimation_information$n_nonconvergence,
        inadmissible = condition$estimation_information$n_inadmissible,
        stringsAsFactors = FALSE
      )
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
  if (!is.character(parameter) || length(parameter) != 1L) {
    cli::cli_abort(
      c(
        "{.arg parameter} must be a character string of length 1:",
        x = paste0("Your {.arg parameter} is ", format_object_type(parameter), " of length ", length(parameter), ".")
      )
    )
  }
  parameter_names_any <- give_powRICLPM_parameter_names_any(object)

  # Combine simulation results across experimental conditions
  d <- do.call(rbind, lapply(object$conditions, function(condition) {

    # Extract and round estimates per condition
    estimates <- if (is.data.frame(condition$estimates)) {
      condition$estimates[condition$estimates$parameter == parameter, -1, drop = FALSE]
    } else {
      data.frame()
    }
    if (nrow(estimates) == 0L) {
      if (parameter %in% parameter_names_any) {
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
    cbind(
      data.frame(
        sample_size = condition$sample_size,
        time_points = condition$time_points,
        ICC = icondition_proportion(condition),
        stringsAsFactors = FALSE
      ),
      icondition_reliability_columns(condition),
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
  parameter_names <- give_powRICLPM_parameter_names(object)
  icheck_give_parameter(parameter, object, what = "uncertainty", parameter_names = parameter_names)

  # Combine data frames across conditions
  d <- do.call(rbind, lapply(object$conditions, function(condition) {

    # Extract rows from uncertainty where parameter matches
    uncertainty_filtered <- condition$MCSEs[which(condition$estimates$parameter == parameter), , drop = FALSE]

    # Create data frame
    cbind(
      data.frame(
        sample_size = condition$sample_size,
        time_points = condition$time_points,
        ICC = icondition_proportion(condition),
        stringsAsFactors = FALSE
      ),
      icondition_reliability_columns(condition),
      uncertainty_filtered
    )
  }))
  return(d)
}

icheck_give_parameter <- function(parameter, object, what = "results",
                                  parameter_names = NULL,
                                  call = rlang::caller_env()) {
  if (is.null(parameter)) {
    cli::cli_abort(
      c(
        "No {.arg parameter} was specified:",
        i = paste0("{.code give(object, '", what, "')} needs to know which parameter to extract."),
        i = "Use {.code give(object, 'names')} to see available parameter names."
      ),
      call = call
    )
  }
  if (!is.character(parameter) || length(parameter) != 1L) {
    cli::cli_abort(
      c(
        "{.arg parameter} must be a character string of length 1:",
        x = paste0("Your {.arg parameter} is ", format_object_type(parameter), " of length ", length(parameter), ".")
      ),
      call = call
    )
  }
  if (is.null(parameter_names)) {
    parameter_names <- give_powRICLPM_parameter_names(object)
  }
  if (!parameter %in% parameter_names) {
    cli::cli_abort(
      c(
        "The requested {.arg parameter} was not found across all experimental conditions:",
        x = paste0("No ", what, " is available for parameter `", parameter, "`."),
        i = "Use {.code give(object, 'names')} to see available parameter names."
      ),
      call = call
    )
  }
}

icondition_proportion <- function(condition) {
  if (!is.null(condition$AF_proportion)) {
    return(condition$AF_proportion)
  }
  condition$ICC
}
