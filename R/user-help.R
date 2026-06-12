#' Check Interpretation of Lagged Effects
#'
#' Write a textual interpretation of the values in `lagged_effects`. This can be used to check if `lagged_effects` has been correctly specified.
#'
#' @inheritParams powRICLPM
#' @param ... Not used.
#'
#' @return No return value, called for side effects.
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

  # Check argument type
  if (!is.matrix(lagged_effects)) {
    stop(rlang::format_error_bullets(c(
      paste0("`", argument_name, "` must be a matrix:"),
      x = paste0("Your `", argument_name, "` is a `", typeof(lagged_effects), "`.")
    )))
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
  if (!is.null(Phi) && !is.null(lagged_effects)) {
    iabort_renamed_argument_conflict("lagged_effects", "Phi")
  }
  if (!is.null(Phi)) {
    return(check_lagged_effects(Phi = Phi))
  }

  check_lagged_effects(lagged_effects)
}
