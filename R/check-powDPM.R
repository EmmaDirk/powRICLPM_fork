#' Check \code{time_points} Argument for the DPM
#'
#' @noRd
icheck_T_DPM <- function(x, arg = rlang::caller_arg(x), call = rlang::caller_env()) {
  if (!all(is.numeric(x))) {
    cli::cli_abort(
      c(
        "{.arg {arg}} must be a vector of integers:",
        x = "Not all elements are numeric."
      ),
      call = call
    )
  }
  if (!all(is.finite(x))) {
    cli::cli_abort(
      c(
        "{.arg {arg}} must be a vector of finite integers:",
        x = "Some elements are `NA`, `NaN`, `Inf`, or `-Inf`."
      ),
      call = call
    )
  }
  if (!all(x %% 1 == 0)) {
    cli::cli_abort(
      c(
        "{.arg {arg}} must be a vector of integers:",
        x = "Not all elements are integers."
      ),
      call = call
    )
  }
  if (any(x < 2)) {
    cli::cli_abort(
      c(
        "Elements in {.arg {arg}} should be larger than 1:",
        i = "The DPM needs at least 2 time points to be identified.",
        x = "You've supplied a number of time points smaller than 2."
      ),
      call = call
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
  invisible(NULL)
}

icheck_DPM_constraints <- function(constraints, call = rlang::caller_env()) {
  constraints <- normalize_constraints(constraints)
  if (has_constraint(constraints, "ME")) {
    cli::cli_abort(
      c(
        "`constraints = 'ME'` cannot be used with `powDPM()`:",
        i = "The DPM does not separate measurement error."
      ),
      call = call
    )
  }
  if ("within" %in% constraints) {
    cli::cli_abort(
      c(
        "`constraints = 'within'` cannot be used with `powDPM()`:",
        i = "Use `constraints = c('lagged', 'residuals')` if you want both DPM lagged effects and observed-level residuals constrained over time."
      ),
      call = call
    )
  }

  invisible(NULL)
}

icheck_DPM_loadings <- function(x, time_points, constraints, arg, t_arg, con_arg, call) {
  is_numeric_vector <- is.numeric(x) && is.null(dim(x))
  is_numeric_matrix <- is.numeric(x) && is.matrix(x)
  if (!is_numeric_vector && !is_numeric_matrix) {
    cli::cli_abort(
      c(
        "{.arg {arg}} must be a numeric vector or numeric matrix:",
        x = paste0("Your {.arg {arg}} is ", format_object_type(x), ".")
      ),
      call = call
    )
  }
  if (is_numeric_vector && length(x) == 0L) {
    cli::cli_abort(
      c(
        "{.arg {arg}} must contain at least one value:",
        x = paste0("Your {.arg {arg}} has length 0.")
      ),
      call = call
    )
  }
  if (!all(is.finite(x) | is.na(x))) {
    cli::cli_abort(
      c(
        "{.arg {arg}} must contain only finite values, except for the required first-wave NA:",
        x = paste0("Your {.arg {arg}} contains: ", format_loadings(x), ".")
      ),
      call = call
    )
  }

  if (is.matrix(x)) {
    if (nrow(x) != 2L) {
      row_label <- if (nrow(x) == 1L) "row" else "rows"
      cli::cli_abort(
        c(
          "{.arg {arg}} must have 2 rows:",
          i = "Rows correspond to the accumulating factors for variables A and B.",
          x = paste0("Your {.arg {arg}} has ", nrow(x), " ", row_label, ".")
        ),
        call = call
      )
    }
    if (ncol(x) != time_points) {
      cli::cli_abort(
        c(
          "{.arg {arg}} must have one column per time point for `model = 'DPM'`:",
          i = "The first column must be `NA` because accumulating factors do not load on wave 1.",
          x = paste0(
            "Your {.arg {arg}} has ", ncol(x),
            " columns, but {.arg {t_arg}} = ", time_points, "."
          )
        ),
        call = call
      )
    }
    if (!all(is.na(x[, 1]) & !is.nan(x[, 1]))) {
      cli::cli_abort(
        c(
          "The first DPM loading column must be `NA`:",
          i = "Accumulating factors load on waves 2 through T, not on wave 1.",
          x = paste0("The first column is ", format_loadings(x[, 1]), ".")
        ),
        call = call
      )
    }
    x <- x[, -1, drop = FALSE]
    first_loadings <- x[, 1]
  } else {
    if (length(x) != time_points) {
      cli::cli_abort(
        c(
          "{.arg {arg}} must have one value per time point for `model = 'DPM'`:",
          i = "The first value must be `NA` because accumulating factors do not load on wave 1.",
          x = paste0("Your {.arg {arg}} has length ", length(x),
                     ", but {.arg {t_arg}} = ", time_points, ".")
        ),
        call = call
      )
    }
    if (!is.na(x[1]) || is.nan(x[1])) {
      cli::cli_abort(
        c(
          "The first DPM loading must be `NA`:",
          i = "Accumulating factors load on waves 2 through T, not on wave 1.",
          x = paste0("The first value is ", x[1], ".")
        ),
        call = call
      )
    }
    x <- x[-1]
    first_loadings <- x[1]
  }
  if (!all(is.finite(x))) {
    cli::cli_abort(
      c(
        "{.arg {arg}} must contain only finite values after the required first-wave NA:",
        x = paste0("Your {.arg {arg}} contains: ", format_loadings(x), ".")
      ),
      call = call
    )
  }
  if (!all(first_loadings == 1)) {
    cli::cli_abort(
      c(
        "The second DPM loading must be 1:",
        i = "The first loading must be `NA`; the second loading is the wave-2 loading.",
        x = paste0("The second loading is ", paste(first_loadings, collapse = ", "), ".")
      ),
      call = call
    )
  }
  if (!has_constraint(constraints, "loadings_free")) {
    cli::cli_abort(
      c(
        "{.arg {arg}} can only be used when accumulating-factor loadings are freely estimated:",
        i = "Use `constraints = 'loadings_free'` or include `'loadings_free'` in a constraint vector.",
        x = paste0("You supplied {.arg {arg}}, but {.arg {con_arg}} = ", format_constraints(constraints), ".")
      ),
      call = call
    )
  }
  invisible(NULL)
}

inormalize_DPM_loadings <- function(loadings, time_points) {
  if (is.null(loadings)) {
    return(matrix(1, nrow = 2, ncol = time_points - 1L))
  }
  if (is.matrix(loadings)) {
    if (ncol(loadings) == time_points) {
      return(unname(loadings[, -1, drop = FALSE]))
    }
    return(unname(loadings))
  }
  if (length(loadings) == time_points) {
    loadings <- loadings[-1]
  }
  unname(rbind(loadings, loadings))
}
