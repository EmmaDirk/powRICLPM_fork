icheck_object_summary <- function(object, arg = rlang::caller_arg(object), call = rlang::caller_env()) {
  if (!inherits(object, "powRICLPM")) {
    cli::cli_abort(
      c(
        "{.arg {arg}} must be of class 'powRICLPM'."
      ),
      call = call
    )
  }
}


icheck_sample_size_summary <- function(sample_size, object, arg = rlang::caller_arg(sample_size), call = rlang::caller_env()) {

  if (length(sample_size) > 1) {
    cli::cli_abort(
      c(
        "{.arg {arg}} must be a single number:",
        "x" = "Your {.arg {arg}} is of length {length(sample_size)}."
      ),
      call = call
    )
  }

  sample_sizes <- lapply(object$conditions, function(x) {x$sample_size})

  if (!any(sample_size == sample_sizes)) {
    cli::cli_abort(
      c(
        "{.arg {arg}} must refer to an experimental condition in the {.cls {class(object)}} object with that sample size:",
        "i" = "The sample size you've indicated is not included in any experimental condition.",
        "x" = "Perhaps you meant any of the following sample sizes?",
        paste(unique(sample_sizes), collapse = ", ")
      ),
      call = call
    )
  }
}


icheck_parameter_summary <- function(x, object, arg = rlang::caller_arg(x), call = rlang::caller_env()) {

  if (length(x) > 1) {
    cli::cli_abort(
      c(
        "{.arg {arg}} must be a character string of length 1:",
        "x" = "Your {.arg {arg}} is of length {length(x)}."
      ),
      call = call
    )
  }

  if (!is.character(x)) {
    cli::cli_abort(
      c(
        "{.arg {arg}} must be a character string:",
        "x" = "Your {.arg {arg}} is of type {typeof(x)}."
      ),
      call = call
    )
  }

  names_parameters <- give_powRICLPM_parameter_names(object)

  if (!any(x == names_parameters)) {
    cli::cli_abort(
      c(
        "The requested {.arg {arg}} was not found across all experimental conditions:",
        "x" = "No summary is available for parameter `{x}`.",
        "i" = "Use `give(object, what = 'names')` to see available parameter names."
      ),
      call = call
    )
  }
}



icheck_time_points_summary <- function(time_points, object, arg = rlang::caller_arg(time_points), call = rlang::caller_env()) {

  if (length(time_points) > 1) {
    cli::cli_abort(
      c(
        "{.arg {arg}} must be a single number:",
        "x" = "Your {.arg {arg}} is of length {length(time_points)}."
      ),
      call = call
    )
  }

  time_points_object <- lapply(object$conditions, function(x) {x$time_points})

  if (!any(time_points == time_points_object)) {
    cli::cli_abort(
      c(
        "{.arg {arg}} must refer to an experimental condition in the {.cls {class(object)}} object with that number of time points:",
        "i" = "The {.arg {arg}} you've indicated is not included in any experimental condition.",
        "x" = "Perhaps you meant any of the following number of time points?",
        paste(unique(time_points_object), collapse = ", ")
      ),
      call = call
    )
  }
}


icheck_ICC_summary <- function(ICC, object, arg = rlang::caller_arg(ICC), call = rlang::caller_env()) {
  if (length(ICC) > 1) {
    cli::cli_abort(
      c(
        "{.arg {arg}} must be a single number:",
        "x" = "Your {.arg {arg}} is of length {length(ICC)}."
      ),
      call = call
    )
  }

  ICCs <- lapply(object$conditions, icondition_proportion)

  if (!any(ICC == ICCs)) {
    cli::cli_abort(
      c(
        "{.arg {arg}} must refer to an experimental condition in the {.cls {class(object)}} object with that ICC:",
        "i" = "The {.arg {arg}} you've indicated is not included in any experimental condition.",
        "x" = "Perhaps you meant any of the following ICCs?",
        paste(unique(ICCs), collapse = ", ")
      ),
      call = call
    )
  }
}


icheck_reliability_summary <- function(reliability, object, arg = rlang::caller_arg(reliability), call = rlang::caller_env()) {

  if (is.character(reliability)) {
    cli::cli_abort(
      c(
        "{.arg {arg}} must be either a single numeric value or a named numeric vector:",
        i = "For variable-specific reliability, use `reliability = c(A = 0.8, B = 0.7)`.",
        i = "Character strings such as \"A = 0.8, B = 0.7\" are not supported."
      ),
      call = call
    )
  }
  if (!is.numeric(reliability) || !is.null(dim(reliability))) {
    cli::cli_abort(
      c(
        "{.arg {arg}} must be either a single numeric value or a named numeric vector:",
        i = "For variable-specific reliability, use `reliability = c(A = 0.8, B = 0.7)`.",
        x = paste0("Your {.arg {arg}} is ", format_object_type(reliability), ".")
      ),
      call = call
    )
  }
  if (length(reliability) == 0L) {
    cli::cli_abort(
      c(
        "{.arg {arg}} must contain one reliability selector:",
        i = "Use a single numeric value or `reliability = c(A = 0.8, B = 0.7)`."
      ),
      call = call
    )
  }
  if (!all(is.finite(reliability))) {
    cli::cli_abort(
      c(
        "{.arg {arg}} must contain only finite, non-missing reliability values:",
        x = paste0("Your {.arg {arg}} contains: ", format_reliability(reliability), ".")
      ),
      call = call
    )
  }
  if (length(reliability) > 1L && !is_combined_reliability_selector(reliability)) {
    cli::cli_abort(
      c(
        "For variable-specific reliability, {.arg {arg}} must specify exactly variables {.val A} and {.val B}:",
        i = "Use `reliability = c(A = 0.8, B = 0.7)`.",
        x = paste0("Your {.arg {arg}} has names: ", format_reliability_names(names(reliability)), ".")
      ),
      call = call
    )
  }

  reliabilities <- vapply(object$conditions, function(x) {as.character(x$reliability)}, character(1))

  if (!any(ireliability_matches(reliability, reliabilities))) {
    cli::cli_abort(
      c(
        "{.arg {arg}} must refer to an experimental condition in the {.cls {class(object)}} object with that reliability:",
        "i" = "The reliability you've indicated is not included in any experimental condition.",
        "x" = "Perhaps you meant one of the following reliability conditions?",
        paste(ireliability_available_display_labels(object), collapse = "; ")
      ),
      call = call
    )
  }
}

is_combined_reliability_selector <- function(reliability) {
  is.numeric(reliability) &&
    is.null(dim(reliability)) &&
    length(reliability) == 2L &&
    !is.null(names(reliability)) &&
    all(nzchar(names(reliability))) &&
    !anyDuplicated(names(reliability)) &&
    setequal(names(reliability), c("A", "B"))
}

format_reliability_names <- function(x) {
  if (is.null(x)) {
    return("<none>")
  }
  x[!nzchar(x)] <- "<empty>"
  paste(x, collapse = ", ")
}


imatch_condition_summary <- function(object, sample_size, time_points, ICC, reliability = NULL,
                                     call = rlang::caller_env()) {
  matches <- Filter(function(x) {
    x$sample_size == sample_size &&
      x$time_points == time_points &&
      icondition_proportion(x) == ICC &&
      (is.null(reliability) || ireliability_matches(reliability, x$reliability))
  }, object$conditions)

  if (length(matches) == 0) {
    supplied_reliability <- if (is.null(reliability)) {
      "not supplied"
    } else {
      reliability
    }
    cli::cli_abort(
      c(
        "No experimental condition matches the supplied condition arguments:",
        i = "Check the combination of {.arg sample_size}, {.arg time_points}, {.arg intraclass_correlation}, and {.arg reliability}.",
        x = paste0(
          "Supplied values: sample_size = ", sample_size,
          ", time_points = ", time_points,
          ", intraclass_correlation = ", ICC,
          ", reliability = ", supplied_reliability, "."
        )
      ),
      call = call
    )
  }

  matching_reliabilities <- unique(vapply(matches, function(x) as.character(x$reliability), character(1)))
  if (is.null(reliability) && length(matching_reliabilities) > 1) {
    cli::cli_abort(
      c(
        "Multiple experimental conditions match the supplied condition arguments:",
        i = "{.arg reliability} is needed to select one condition.",
        i = "For variable-specific matrix reliability, use `reliability = c(A = ..., B = ...)`.",
        x = paste0(
          "Matching reliability conditions are: ",
          paste(vapply(matching_reliabilities, ireliability_display_label, character(1)), collapse = "; "),
          "."
        )
      ),
      call = call
    )
  }

  matches[[1]]
}

