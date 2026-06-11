#' Check \code{what} Argument
#'
#' \code{check_give} checks the \code{what} argument of \code{give()} by testing if it is a character string of length 1.
#'
#' @param x A character string of length 1.
#'
#' @noRd
icheck_what_give <- function(x, arg = rlang::caller_arg(x), call = rlang::caller_env()) {
  if (!is.character(x)) {
    cli::cli_abort(
      c(
        "{.arg {arg}} must be a character string:",
        "x" = "Your {.arg {arg}} is a {.cls {typeof(x)}}."
      ),
      call = call
    )
  }

  if (length(x) > 1) {
    cli::cli_abort(
      c(
        "{.arg {arg}} must be of length 1:",
        "x" = "Your {.arg {arg}} is of length {length(x)}."
      )
    )
  }

  x <- normalize_intraclass_correlation_value(x)
  if (!any(x == c(
    "conditions", "estimation_problems", "results", "names",
    "sample_size", "time_points", "intraclass_correlation",
    "reliability", "uncertainty"
  ))) {
    cli::cli_abort(
      c(
        "{.arg {arg}} is not an accepted input. Please change it to `conditions`, `estimation_problems`, `results`, `names`, `sample_size`, `time_points`, `intraclass_correlation`, `reliability`, or `uncertainty`.",
        "x" = "Your {.arg {arg}} is {.val {x}}."
      ),
      call = call
    )
  }
}


