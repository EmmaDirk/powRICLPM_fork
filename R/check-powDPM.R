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
  if (any(x < 3)) {
    cli::cli_abort(
      c(
        "Elements in {.arg {arg}} should be larger than 2:",
        i = "The DPM needs at least 3 time points before identification constraints can be evaluated.",
        x = "You've supplied a number of time points smaller than 3."
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

DPM_identification_table <- data.frame(
  estimate_ME = c(
    rep(FALSE, 10),
    rep(TRUE, 20)
  ),
  constraints_key = c(
    "none", "lagged", "residuals", "lagged+residuals", "stationarity",
    "none", "lagged", "residuals", "lagged+residuals", "stationarity",
    "none", "none+ME", "lagged", "lagged+ME", "residuals", "residuals+ME",
    "lagged+residuals", "lagged+residuals+ME", "stationarity", "stationarity+ME",
    "none", "none+ME", "lagged", "lagged+ME", "residuals", "residuals+ME",
    "lagged+residuals", "lagged+residuals+ME", "stationarity", "stationarity+ME"
  ),
  AF_loadings_free = c(
    rep(FALSE, 5), rep(TRUE, 5),
    rep(FALSE, 10), rep(TRUE, 10)
  ),
  min_waves = c(
    4, 3, 3, 3, 3,
    4, 4, 4, 3, 3,
    5, 4, 4, 4, 4, 4, 4, 3, 3, 3,
    5, 5, 4, 4, 5, 4, 4, 3, 3, 3
  ),
  stringsAsFactors = FALSE
)

classify_DPM_spec <- function(estimate_ME, constraints) {
  constraints <- normalize_constraints(constraints)
  AF_loadings_free <- has_constraint(constraints, "AF_loadings_free")
  base_constraints <- constraints[!constraints %in% c("AF_loadings_free", "ME")]
  base_key <- if (length(base_constraints) == 0L || has_constraint(base_constraints, "none")) {
    "none"
  } else if (setequal(base_constraints, c("lagged", "residuals"))) {
    "lagged+residuals"
  } else if (length(base_constraints) == 1L) {
    base_constraints
  } else {
    paste(sort(base_constraints), collapse = "+")
  }
  if (has_constraint(constraints, "ME")) {
    base_key <- paste0(base_key, "+ME")
  }
  list(
    estimate_ME = isTRUE(estimate_ME),
    constraints_key = base_key,
    AF_loadings_free = AF_loadings_free
  )
}

lookup_DPM_min_waves <- function(estimate_ME, constraints) {
  spec <- classify_DPM_spec(estimate_ME, constraints)
  row <- DPM_identification_table[
    DPM_identification_table$estimate_ME == spec$estimate_ME &
      DPM_identification_table$constraints_key == spec$constraints_key &
      DPM_identification_table$AF_loadings_free == spec$AF_loadings_free,
    ,
    drop = FALSE
  ]
  if (nrow(row) != 1L) {
    cli::cli_abort(
      c(
        "No DPM identification rule is available for this specification:",
        x = paste0(
          "estimate_ME = ", spec$estimate_ME,
          ", constraints = ", spec$constraints_key,
          ", AF_loadings_free = ", spec$AF_loadings_free, "."
        )
      )
    )
  }
  row$min_waves
}

icheck_DPM_identification <- function(time_points, reliability, estimate_ME, constraints,
                                      call = rlang::caller_env()) {
  min_waves <- lookup_DPM_min_waves(estimate_ME, constraints)
  if (any(time_points < min_waves)) {
    supplied <- min(time_points[time_points < min_waves])
    cli::cli_abort(
      c(
        paste0(
          "The requested DPM specification is not identified with ",
          supplied, " measurement waves; this specification requires at least ",
          min_waves, " waves."
        ),
        i = "Use more waves or impose valid identifying constraints."
      ),
      call = call
    )
  }
  invisible(NULL)
}

icheck_AF_proportion <- function(x, arg = rlang::caller_arg(x), call = rlang::caller_env()) {
  if (!is.numeric(x) || !is.null(dim(x))) {
    cli::cli_abort(
      c(
        "{.arg {arg}} should be a numeric vector:",
        x = paste0("Your {.arg {arg}} is ", format_object_type(x), ".")
      ),
      call = call
    )
  }
  if (!all(is.finite(x))) {
    cli::cli_abort(
      c(
        "{.arg {arg}} must contain only finite values:",
        x = paste0("Your {.arg {arg}} contains: ", paste(as.character(x), collapse = ", "), ".")
      ),
      call = call
    )
  }
  if (!all(x > 0 & x < 1)) {
    cli::cli_abort(
      c(
        "Elements in {.arg {arg}} must be between 0 and 1:",
        x = "Some elements in {.arg {arg}} are smaller than 0 or larger than 1."
      ),
      call = call
    )
  }
  if (!all(x < .99)) {
    cli::cli_abort(
      c(
        "Elements in {.arg {arg}} are very close to 1:",
        i = "This can lead to problems with estimation. Use a maximum of .99."
      ),
      call = call
    )
  }
  invisible(NULL)
}

icheck_AF_cor <- function(x, arg = rlang::caller_arg(x), call = rlang::caller_env()) {
  icheck_DPM_cor(x, arg = arg, call = call)
}

icheck_dynamics_cor <- function(x, arg = rlang::caller_arg(x), call = rlang::caller_env()) {
  icheck_DPM_cor(x, arg = arg, call = call)
}

icheck_DPM_cor <- function(x, arg, call) {
  if (!is.numeric(x) || !is.null(dim(x))) {
    cli::cli_abort(
      c(
        "{.arg {arg}} must be a numeric value:",
        x = paste0("Your {.arg {arg}} is ", format_object_type(x), ".")
      ),
      call = call
    )
  }
  if (length(x) != 1L) {
    cli::cli_abort(
      c(
        "{.arg {arg}} should be of length 1:",
        x = paste0("You've provided ", length(x), " values.")
      ),
      call = call
    )
  }
  if (!is.finite(x)) {
    cli::cli_abort(
      c(
        "{.arg {arg}} must be finite:",
        x = paste0("Your {.arg {arg}} is ", format_scalar_value(x), ".")
      ),
      call = call
    )
  }
  if (x < -1 || x > 1) {
    cli::cli_abort(
      c(
        "A correlation must be between -1 and 1:",
        x = paste0("Your {.arg {arg}} was ", x, ".")
      ),
      call = call
    )
  }
  invisible(NULL)
}

icheck_DPM_lagged_effects <- function(x, arg = rlang::caller_arg(x), call = rlang::caller_env()) {
  if (!is.matrix(x)) {
    cli::cli_abort(
      c(
        "{.arg {arg}} must be a matrix:",
        x = paste0("Your {.arg {arg}} is ", format_object_type(x), ".")
      ),
      call = call
    )
  }
  if (!is.numeric(x)) {
    cli::cli_abort(
      c(
        "{.arg {arg}} must be a numeric matrix:",
        x = paste0("Your {.arg {arg}} is a matrix of ", typeof(x), ".")
      ),
      call = call
    )
  }
  if (!identical(dim(x), c(2L, 2L))) {
    cli::cli_abort(
      c(
        "{.arg {arg}} must be a 2 by 2 matrix:",
        x = paste0("Your {.arg {arg}} has dimensions ", paste(dim(x), collapse = " by "), ".")
      ),
      call = call
    )
  }
  if (!all(is.finite(x))) {
    cli::cli_abort(
      c(
        "{.arg {arg}} must contain only finite values.",
        x = paste0("Your {.arg {arg}} contains: ", paste(as.character(x), collapse = ", "), ".")
      ),
      call = call
    )
  }
  if (!is_unit(x)) {
    cli::cli_abort(
      c(
        "{.arg {arg}} must specify a stationary process:",
        i = "This is checked by testing if the eigenvalues of {.arg {arg}} lie within the unit circle.",
        x = "Eigenvalues of {.arg {arg}} are not within the unit circle. Try smaller lagged effects?"
      ),
      call = call
    )
  }
  invisible(NULL)
}

icheck_DPM_bounds <- function(bounds, arg = rlang::caller_arg(bounds), call = rlang::caller_env()) {
  if (!is.logical(bounds) || length(bounds) != 1L || is.na(bounds)) {
    cli::cli_abort(
      c(
        "{.arg {arg}} must be TRUE or FALSE:",
        x = paste0("Your {.arg {arg}} is ", format_object_type(bounds), " of length ", length(bounds), ".")
      ),
      call = call
    )
  }
  invisible(NULL)
}

icheck_DPM_constraints <- function(constraints, estimate_ME = FALSE, time_points = NULL,
                                   call = rlang::caller_env()) {
  valid_constraints <- c(
    "none", "lagged", "residuals", "stationarity", "ME", "AF_loadings_free"
  )
  if (!is.character(constraints)) {
    cli::cli_abort(
      c(
        "{.arg constraints} must be a character vector:",
        x = paste0("Your {.arg constraints} is ", typeof(constraints), ".")
      ),
      call = call
    )
  }
  if (length(constraints) == 0L) {
    cli::cli_abort(
      c(
        "{.arg constraints} must contain at least one constraint option:",
        i = "Use 'none' when no constraints should be imposed."
      ),
      call = call
    )
  }
  constraints <- normalize_constraints(constraints)
  if (!all(constraints %in% valid_constraints)) {
    cli::cli_abort(
      c(
        "{.arg constraints} contains invalid DPM constraints:",
        i = paste0("Valid DPM constraints are: ", paste(valid_constraints, collapse = ", "), "."),
        x = paste0("Your {.arg constraints} is ", format_constraints(constraints), ".")
      ),
      call = call
    )
  }
  if (has_constraint(constraints, "none") && length(constraints) > 1L) {
    cli::cli_abort(
      c(
        "{.arg constraints} can only include 'none' by itself:",
        x = paste0("Your {.arg constraints} is ", format_constraints(constraints), ".")
      ),
      call = call
    )
  }
  if (has_constraint(constraints, "ME") && !isTRUE(estimate_ME)) {
    cli::cli_abort(
      c(
        "{.arg constraints} can only be set to 'ME' when `estimate_ME = TRUE`:",
        x = "Your `estimate_ME` is FALSE."
      ),
      call = call
    )
  }
  if (has_constraint(constraints, "stationarity") &&
      any(constraints %in% c("lagged", "residuals"))) {
    cli::cli_abort(
      c(
        "{.arg constraints} cannot combine 'stationarity' with lagged or residual constraints:",
        i = "'stationarity' already imposes its own DPM lagged-effect and residual restrictions.",
        x = paste0("Your {.arg constraints} is ", format_constraints(constraints), ".")
      ),
      call = call
    )
  }
  if (!is.null(time_points) && has_constraint(constraints, "AF_loadings_free") &&
      any(time_points < 3L)) {
    cli::cli_abort(
      c(
        "`AF_loadings_free` requires at least 3 time points:",
        i = "The wave-2 accumulating-factor loading is fixed to 1, so free loadings start at wave 3."
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
  if (!all(is.finite(x))) {
    cli::cli_abort(
      c(
        "{.arg {arg}} must contain only finite values:",
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
    if (ncol(x) != time_points - 1L) {
      cli::cli_abort(
        c(
          "{.arg {arg}} must have one column for each wave from 2 through T:",
          i = "DPM accumulating-factor loadings are supplied for waves 2 through T only.",
          x = paste0(
            "{.arg {arg}} has ", ncol(x),
            " columns, but {.arg {t_arg}} = ", time_points,
            " requires ", time_points - 1L,
            " loading columns because DPM accumulating-factor loadings are specified for waves 2 through ",
            time_points, "."
          )
        ),
        call = call
      )
    }
    first_loadings <- x[, 1]
  } else {
    if (length(x) != time_points - 1L) {
      cli::cli_abort(
        c(
          "{.arg {arg}} must have one value for each wave from 2 through T:",
          i = "DPM accumulating-factor loadings are supplied for waves 2 through T only.",
          x = paste0("{.arg {arg}} has length ", length(x),
                     ", but {.arg {t_arg}} = ", time_points,
                     " requires ", time_points - 1L,
                     " loadings because DPM accumulating-factor loadings are specified for waves 2 through ",
                     time_points, ".")
        ),
        call = call
      )
    }
    first_loadings <- x[1]
  }
  if (!all(first_loadings == 1)) {
    cli::cli_abort(
      c(
        "The first DPM loading must be 1:",
        i = "The first supplied DPM loading is the wave-2 accumulating-factor loading.",
        x = paste0("The first supplied loading is ", paste(first_loadings, collapse = ", "), ".")
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
    return(unname(loadings))
  }
  unname(rbind(loadings, loadings))
}

detect_DPM_misspecification <- function(reliability, loadings, time_points, constraints,
                                        estimate_ME) {
  reliability_conditions <- ireliability_conditions(reliability)
  records <- list()
  for (i in seq_len(nrow(reliability_conditions))) {
    for (time_points_i in time_points) {
      reliability_matrix <- inormalize_reliability(
        reliability_conditions$reliability_spec[[i]],
        time_points_i
      )
      loadings_matrix <- inormalize_DPM_loadings(loadings, time_points_i)
      records[[length(records) + 1L]] <- detect_DPM_misspecification_condition(
        reliability_matrix = reliability_matrix,
        loadings = loadings_matrix,
        constraints = constraints,
        estimate_ME = estimate_ME
      )
    }
  }
  combine_DPM_misspecification(records)
}

detect_DPM_misspecification_condition <- function(reliability_matrix, loadings,
                                                 constraints, estimate_ME) {
  restrictive_reasons <- character()
  general_reasons <- character()

  if (any(reliability_matrix < 1) && !isTRUE(estimate_ME)) {
    restrictive_reasons <- c(
      restrictive_reasons,
      "Generated measurement error is ignored in the fitted model."
    )
  }
  if (!any(reliability_matrix < 1) && isTRUE(estimate_ME)) {
    general_reasons <- c(
      general_reasons,
      "Measurement error is estimated although the data-generating model has no measurement error."
    )
  }
  generated_AF_loadings_general <- any(loadings != 1)
  estimated_AF_loadings_free <- has_constraint(constraints, "AF_loadings_free")
  if (generated_AF_loadings_general && !estimated_AF_loadings_free) {
    restrictive_reasons <- c(
      restrictive_reasons,
      "Generated AF loadings vary but fitted AF loadings are fixed."
    )
  }
  if (!generated_AF_loadings_general && estimated_AF_loadings_free) {
    general_reasons <- c(
      general_reasons,
      "AF loadings are freely estimated although the data-generating model fixes them to 1."
    )
  }

  list(
    restrictive = length(restrictive_reasons) > 0L,
    general = length(general_reasons) > 0L,
    reasons = unique(c(restrictive_reasons, general_reasons)),
    restrictive_reasons = unique(restrictive_reasons),
    general_reasons = unique(general_reasons)
  )
}

combine_DPM_misspecification <- function(records) {
  restrictive_reasons <- unique(unlist(lapply(records, function(x) {
    x$restrictive_reasons
  }), use.names = FALSE))
  general_reasons <- unique(unlist(lapply(records, function(x) {
    x$general_reasons
  }), use.names = FALSE))
  list(
    restrictive = length(restrictive_reasons) > 0L,
    general = length(general_reasons) > 0L,
    reasons = unique(c(restrictive_reasons, general_reasons)),
    restrictive_reasons = restrictive_reasons,
    general_reasons = general_reasons
  )
}

confirm_restrictive_misspecification <- function(misspecification,
                                                 call = rlang::caller_env()) {
  if (!isTRUE(misspecification$restrictive)) {
    return(invisible(TRUE))
  }
  message <- paste(
    "The estimation model is more restrictive than the data-generating model; power may be overestimated and cross-lagged or autoregressive estimates may be biased.",
    "",
    sep = "\n"
  )
  if (!interactive()) {
    cli::cli_abort(
      c(
        message,
        x = "Restrictive DPM misspecification requires interactive confirmation.",
        i = "Run interactively and type `YES` to continue."
      ),
      call = call
    )
  }
  cat(message)
  response <- readline("Type YES to continue: ")
  if (!identical(response, "YES")) {
    cli::cli_abort("DPM simulation aborted because `YES` was not entered.", call = call)
  }
  invisible(TRUE)
}
