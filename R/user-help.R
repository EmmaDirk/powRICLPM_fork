#' Check Interpretation of Lagged Effects
#'
#' Write a textual interpretation of the values in `lagged_effects`. This can be used to check if `lagged_effects` has been correctly specified.
#'
#' @inheritParams powRICLPM
#' @param ... Additional arguments are not allowed.
#'
#' @return Invisibly returns \code{NULL}. The function is used for printing an
#'   interpretation of \code{lagged_effects}.
#' @export
#'
#' @examples
#' # Correctly specified lagged effects
#' lagged_effects1 <- matrix(c(.4, .1, .2, .3), ncol = 2, byrow = TRUE)
#' check_lagged_effects(lagged_effects1)
#'
#' # Lagged effects with too large standardized effects
#' lagged_effects2 <- matrix(c(.6, .5, .4, .7), ncol = 2, byrow = TRUE)
#' lagged_effects2 <- check_lagged_effects(lagged_effects2)
check_lagged_effects <- function(lagged_effects = NULL, Phi = NULL, ...) {
  argument_name <- "lagged_effects"
  if (!is.null(Phi) && !is.null(lagged_effects)) {
    iabort_renamed_argument_conflict("lagged_effects", "Phi")
  }
  if (!is.null(Phi)) {
    lagged_effects <- Phi
    argument_name <- "Phi"
  }

  dots <- list(...)
  if (length(dots) > 0) {
    dot_names <- names(dots)
    dot_names[dot_names == ""] <- "<unnamed>"
    cli::cli_abort(
      c(
        "Unexpected argument in {.fn check_lagged_effects}:",
        x = paste0("Unknown argument(s): ", paste0(dot_names, collapse = ", "), ".")
      )
    )
  }

  iwrite_lagged_effects_check(lagged_effects, argument_name = argument_name)
}

iwrite_lagged_effects_check <- function(lagged_effects, argument_name,
                                        call = rlang::caller_env()) {
  # Check argument type
  if (!is.matrix(lagged_effects)) {
    cli::cli_abort(
      c(
        paste0("`", argument_name, "` must be a matrix:"),
        x = paste0("Your `", argument_name, "` is a `", typeof(lagged_effects), "`.")
      ),
      call = call
    )
  }

  # Interpretation cross-lagged effects
  writeLines(
    rlang::format_error_bullets(c(
      paste0("According to `", argument_name, "`, the lagged effects are:"),
      "*" = paste0("Autoregressive effect of A: ", lagged_effects[1, 1]),
      "*" = paste0("Autoregressive effect of B: ", lagged_effects[2, 2]),
      "*" = paste0("Cross-lagged effect of A -> B: ", lagged_effects[2, 1]),
      "*" = paste0("Cross-lagged effect of B -> A: ", lagged_effects[1, 2])
    ))
  )

  # Check if positive definite
  if (!is_unit(lagged_effects)) {
    writeLines(
        rlang::format_error_bullets(c(
        paste0("\nHowever, `", argument_name, "` must specify a stationary process:"),
        i = paste0("This is checked by testing if the eigenvalues of `", argument_name, "` lie within unit circle."),
        x = paste0("The eigenvalues of `", argument_name, "` are not within unit circle. Try out smaller lagged effects?")
      ))
    )
  }
}

#' @rdname check_lagged_effects
#'
#' @inheritParams powRICLPM
#' @param Phi Alternative name for \code{lagged_effects}.
#'
#' @details `check_Phi()` is retained for users who prefer the `Phi` notation.
#' @export
#'
#' @examples
#' lagged_effects <- matrix(c(.4, .1, .2, .3), ncol = 2, byrow = TRUE)
#' check_Phi(lagged_effects)
check_Phi <- function(lagged_effects = NULL, Phi = NULL) {
  argument_name <- "Phi"
  if (!is.null(Phi) && !is.null(lagged_effects)) {
    iabort_renamed_argument_conflict("lagged_effects", "Phi")
  }
  if (is.null(Phi)) {
    Phi <- lagged_effects
    if (!is.null(lagged_effects)) {
      argument_name <- "lagged_effects"
    }
  }

  iwrite_lagged_effects_check(Phi, argument_name = argument_name)
}

#' Check Interpretation of Factor Loadings
#'
#' Write a textual interpretation of the values in `loadings`. This can be used
#' to check if time-varying random-intercept loadings have been correctly
#' specified for lavaan data generation, or if accumulating-factor loadings have
#' been correctly specified for the DPM. Time-varying loadings are not supported
#' for Mplus.
#'
#' @param loadings A \code{numeric} vector or matrix specifying
#'   random-intercept loadings in the lavaan data-generating model. A vector of
#'   length \code{time_points} is applied to both variables. A matrix must have
#'   two rows, one for each variable, and one column per time point. The first
#'   loading must be 1, so later values are interpreted relative to the first
#'   occasion. For \code{model = "DPM"}, vectors or matrices specify waves 2
#'   through T only; the first supplied loading is the wave-2 loading and must
#'   be 1.
#' @param time_points (optional) A single \code{integer} indicating the number
#'   of time points. If omitted, this is inferred from \code{loadings}.
#' @param model A \code{character} string, either \code{"RICLPM"} or
#'   \code{"DPM"}, indicating how \code{loadings} should be interpreted.
#' @param ... Additional arguments are not allowed.
#'
#' @return Invisibly returns \code{NULL}. The function is used for printing an
#'   interpretation of \code{loadings}.
#' @export
#'
#' @examples
#' # Same random-intercept loadings for A and B
#' loadings1 <- c(1, .8, 1.1, .9)
#' check_loadings(loadings1)
#'
#' # Different random-intercept loadings for A and B
#' loadings2 <- matrix(c(1, .8, 1.1, .9, 1, 1.2, .7, 1), nrow = 2, byrow = TRUE)
#' check_loadings(loadings2, time_points = 4)
check_loadings <- function(loadings = NULL, time_points = NULL, model = "RICLPM", ...) {
  model <- icheck_model(model)
  dots <- list(...)
  if (length(dots) > 0) {
    dot_names <- names(dots)
    dot_names[dot_names == ""] <- "<unnamed>"
    cli::cli_abort(
      c(
        "Unexpected argument in {.fn check_loadings}:",
        x = paste0("Unknown argument(s): ", paste0(dot_names, collapse = ", "), ".")
      )
    )
  }

  if (is.null(loadings)) {
    cli::cli_abort(
      c(
        "`loadings` must be supplied:",
        x = "Your `loadings` is `NULL`."
      )
    )
  }
  if (is.numeric(loadings) && is.null(dim(loadings)) && length(loadings) == 0L) {
    cli::cli_abort(
      c(
        "`loadings` must contain at least one value:",
        x = "Your `loadings` has length 0."
      )
    )
  }

  if (missing(time_points) || is.null(time_points)) {
    time_points <- iinfer_loadings_time_points(loadings, model)
  }

  icheck_loadings(
    loadings,
    time_points,
    software = "lavaan",
    constraints = if (model == "DPM") "AF_loadings_free" else "RI_loadings_free",
    model = model
  )
  iwrite_loadings_check(loadings, time_points, model)
  invisible(NULL)
}

iinfer_loadings_time_points <- function(loadings, model = "RICLPM") {
  if (is.matrix(loadings)) {
    if (model == "DPM") {
      return(ncol(loadings) + 1L)
    }
    return(ncol(loadings))
  }
  if (model == "DPM") {
    return(length(loadings) + 1L)
  }
  length(loadings)
}

iwrite_loadings_check <- function(loadings, time_points, model = "RICLPM") {
  loadings_matrix <- inormalize_loadings(loadings, time_points, model)
  same_loadings <- identical(loadings_matrix[1, ], loadings_matrix[2, ])
  factor_label <- if (model == "DPM") "accumulating-factor" else "random-intercept"
  latent_names <- if (model == "DPM") c("AF_A", "AF_B") else c("RI_A", "RI_B")
  waves <- if (model == "DPM") 2:time_points else seq_len(time_points)

  if (same_loadings) {
    loading_lines <- paste0(
      latent_names[1], " loads on A",
      waves,
      " and ", latent_names[2], " loads on B",
      waves,
      " with ",
      loadings_matrix[1, ],
      "."
    )

    writeLines(
      rlang::format_error_bullets(c(
        paste0("According to `loadings`, the ", factor_label, " loadings in the data-generating model are:"),
        stats::setNames(loading_lines, rep("*", length(loading_lines)))
      ))
    )
    return(invisible(NULL))
  }

  A_lines <- paste0(latent_names[1], " loads on A", waves, " with ", loadings_matrix[1, ], ".")
  B_lines <- paste0(latent_names[2], " loads on B", waves, " with ", loadings_matrix[2, ], ".")

  writeLines(
    rlang::format_error_bullets(c(
      paste0("According to `loadings`, the ", factor_label, " loadings in the data-generating model are:"),
      stats::setNames(A_lines, rep("*", length(A_lines)))
    ))
  )
  writeLines("")
  writeLines(rlang::format_error_bullets(
    stats::setNames(B_lines, rep("*", length(B_lines)))
  ))
  invisible(NULL)
}

#' Check Interpretation of Reliability
#'
#' Write a textual interpretation of the values in `reliability`. This can be
#' used to check if reliability conditions have been correctly specified for
#' data generation.
#'
#' @param reliability A \code{numeric} value, vector, or matrix specifying
#'   reliability in the lavaan data-generating model. A single value is applied
#'   to both variables at every time point. A vector supplies multiple
#'   reliability conditions, with each value applied to both variables at every
#'   time point. A matrix must have two rows, one for each variable, and one
#'   column per reliability condition.
#' @param time_points (optional) A single \code{integer} indicating the number
#'   of time points. Reliability values are applied to every time point.
#' @param ... Additional arguments are not allowed.
#'
#' @return Invisibly returns \code{NULL}. The function is used for printing an
#'   interpretation of \code{reliability}.
#' @export
#'
#' @examples
#' # Compare multiple reliabilities for A and B
#' reliability1 <- c(.8, .7, .9)
#' check_reliability(reliability1)
#'
#' # Compare reliability conditions that differ between A and B
#' reliability2 <- matrix(c(.8, .9, .7, .85), nrow = 2, byrow = TRUE)
#' check_reliability(reliability2, time_points = 4)
check_reliability <- function(reliability = NULL, time_points = NULL, ...) {
  reliability_expr <- substitute(reliability)
  time_points_missing <- missing(time_points)
  dots <- list(...)
  if (length(dots) > 0) {
    dot_names <- names(dots)
    dot_names[dot_names == ""] <- "<unnamed>"
    cli::cli_abort(
      c(
        "Unexpected argument in {.fn check_reliability}:",
        x = paste0("Unknown argument(s): ", paste0(dot_names, collapse = ", "), ".")
      )
    )
  }

  if (is.null(reliability)) {
    cli::cli_abort(
      c(
        "`reliability` must be supplied:",
        x = "Your `reliability` is `NULL`."
      )
    )
  }

  icheck_reliability_matrix_call(
    expr = reliability_expr,
    env = parent.frame(),
    call = rlang::caller_env()
  )

  if (missing(time_points) || is.null(time_points)) {
    if (!time_points_missing && is.null(time_points)) {
      cli::cli_alert_info("`time_points = NULL` was supplied, so reliability is interpreted as applying to every time point. You can omit `time_points` for the same behavior.")
    }
    time_points <- 1L
  } else if (length(time_points) != 1L) {
    cli::cli_abort(
      c(
        "{.arg time_points} must be one positive whole number:",
        x = paste0("length({.arg time_points}) = ", length(time_points), ".")
      )
    )
  } else if (!is.numeric(time_points) ||
             is.na(time_points) ||
             !is.finite(time_points) ||
             time_points %% 1 != 0 ||
             time_points < 1) {
    cli::cli_abort(
      c(
        "{.arg time_points} must be one positive whole number:",
        x = paste0("{.arg time_points} = ", format_scalar_value(time_points), ".")
      )
    )
  }

  icheck_rel(reliability, time_points, software = "lavaan")
  iwrite_reliability_check(reliability, time_points)
  invisible(NULL)
}

iwrite_reliability_check <- function(reliability, time_points) {
  reliability_conditions <- ireliability_conditions(reliability)

  reliability_lines <- vapply(seq_len(nrow(reliability_conditions)), function(i) {
    reliability_matrix <- inormalize_reliability(
      reliability_conditions$reliability_spec[[i]],
      time_points[1]
    )

    if (identical(reliability_matrix[1, ], reliability_matrix[2, ])) {
      paste0(
        "Condition ", i, ": Variables A and B have reliability ",
        reliability_matrix[1, 1],
        " at every time point."
      )
    } else {
      paste0(
        "Condition ", i, ": ",
        ireliability_condition_label(reliability_matrix),
        "."
      )
    }
  }, character(1))

  writeLines(
    rlang::format_error_bullets(c(
      "According to `reliability`, the reliabilities in the data-generating model are:",
      stats::setNames(reliability_lines, rep("*", length(reliability_lines)))
    ))
  )
  invisible(NULL)
}
