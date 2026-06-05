#' Check Interpretation of Lagged Effects
#'
#' Write a textual interpretation of the values in `lagged_effects`. This can be used to check if `lagged_effects` has been correctly specified.
#'
#' @inheritParams powRICLPM
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
check_lagged_effects <- function(lagged_effects) {
  # Check argument type
  if (!is.matrix(lagged_effects)) {
    stop(rlang::format_error_bullets(c(
      "`lagged_effects` must be a matrix:",
      x = paste0("Your `lagged_effects` is a `", typeof(lagged_effects), "`.")
    )))
  }

  # Interpretation cross-lagged effects
  writeLines(
    rlang::format_error_bullets(c(
      "According to `lagged_effects`, the lagged effects are:",
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
        "\nHowever, `lagged_effects` must specify a stationary process:",
        i = "This is checked by testing if the eigenvalues of `lagged_effects` lie within unit circle.",
        x = "The eigenvalues of `lagged_effects` are not within unit circle. Try out smaller lagged effects?"
      ))
    )
  }
}

#' @rdname check_lagged_effects
#'
#' @inheritParams powRICLPM
#' @param Phi Deprecated. Use \code{lagged_effects} instead.
#'
#' @details `check_Phi()` is deprecated. Use [check_lagged_effects()] instead.
#' @export
#'
#' @examples
#' # Deprecated; use check_lagged_effects() instead.
#' lagged_effects <- matrix(c(.4, .1, .2, .3), ncol = 2, byrow = TRUE)
#' check_Phi(lagged_effects)
check_Phi <- function(lagged_effects = NULL, Phi = NULL) {
  if (!is.null(Phi) && !is.null(lagged_effects)) {
    iabort_renamed_argument_conflict("lagged_effects", "Phi")
  }
  if (!is.null(Phi)) {
    lagged_effects <- Phi
  }

  check_lagged_effects(lagged_effects)
}
