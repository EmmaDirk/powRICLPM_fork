#' Check \code{parameter} argument
#'
#' \code{check_parameter()} tests if a parameter was specified.
#'
#' @param x A character string.
#'
#' @noRd
icheck_plot_parameter <- function(parameter, object, arg = rlang::caller_arg(parameter), call = rlang::caller_env()) {
  if (is.null(parameter)) {
    cli::cli_abort(
      c(
        "No `parameter` was specified:",
        i = "`plot()` needs to know which specific parameter to create a plot for."
      )
    )
  }

  if (length(parameter) > 1) {
    cli::cli_abort(
      c(
        "{.arg {arg}} must be a character vector of size 1:",
        "x" = "Your {.arg {arg}} is of length {length(parameter)}."
      ),
      call = call
    )
  }

  if (!is.character(parameter)) {
    cli::cli_abort(
      c(
        "{.arg {arg}} must be a character string:",
        x = paste0("Your {.arg {arg}} is of type ", typeof(parameter), ".")
      )
    )
  }

  parameter_names <- give_powRICLPM_parameter_names(object)

  if (!any(parameter == parameter_names)) {
    cli::cli_abort(
      c(
        x = "Your {.arg {arg}} is not available across all experimental conditions.",
        i = "Perhaps use `give(object, what = 'names')` to get an overview of parameter names in the `powRICLPM` object."
      )
    )
  }
}

#' Check \code{y} argument
#'
#' \code{check_y()} tests if a \code{y} represents a valid outcome for plotting \code{powRICLPM} results.
#'
#' @param y Character string, specifying which outcome is plotted on the y-axis.
#'
#' @noRd
icheck_y <- function(x, arg = rlang::caller_arg(x), call = rlang::caller_env()) {
  if (length(x) > 1) {
    cli::cli_abort(
      c(
        "You can only plot a single outcome on the y-axis:",
        x = paste0("Your {.arg {arg}} contains {length(x)} outcomes.")
      )
    )
  }
  if (identical(x, "minimum")) {
    cli::cli_abort(
      c(
        "{.arg {arg}} cannot be 'minimum':",
        i = "Monte Carlo standard errors are not computed for the minimum estimate."
      ),
      call = call
    )
  }
  if (!any(x == c("power", "coverage", "accuracy", "MSE", "bias", "average", "EmpSE", "SD", "SEAvg"))) {
    cli::cli_abort(
      c(
        "{.arg {arg}} must be 'power', 'coverage', 'accuracy', 'MSE', 'bias', 'average', 'EmpSE', or 'SEAvg':",
        i = "'SD' is still accepted as an alias for 'EmpSE'.",
        x = paste0("Your {.arg {arg}} is '", x, "'.")
      )
    )
  }
}


icheck_plot_options <- function(x, arg = rlang::caller_arg(x), call = rlang::caller_env()) {
  if (!is.character(x)) {
    cli::cli_abort(
      c(
        "{.arg {arg}} must be a character string:",
        x = "Your {.arg {arg}} is a {.cls {typeof(x)}}."
      ),
      call = call
    )
  }

  if (length(x) > 1) {
    cli::cli_abort(
      c(
        "You can only set a single variable for {.arg {arg}}:",
        x = paste0("Your {.arg {arg}} contains {length(x)} variables")
      )
    )
  }
  x <- normalize_intraclass_correlation_value(x)
  if (!any(x == c("time_points", "intraclass_correlation", "AF_proportion", "reliability", "none"))) {
    cli::cli_abort(
      c(
        "{.arg {arg}} must be 'time_points', 'intraclass_correlation', 'ICC', 'AF_proportion', 'reliability', or 'none':",
        x = "Your {.arg {arg}} is {.val {x}}."
      )
    )
  }
}

