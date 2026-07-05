#' @title
#' Print \code{powRICLPM} Object
#'
#' @description
#' \code{print.powRICLPM} prints a table listing all experimental conditions contained in the \code{powRICLPM} object, as well as the frequency of the estimation problems that occurred in each.
#'
#' @param x A \code{powRICLPM} object.
#' @param ... (don't use)
#'
#' @return No return value, called for side effects.
#'
#' @method print powRICLPM
#' @export
print.powRICLPM <- function(x, ...) {
  icc_table_label <- iicc_table_name(object = x)

  # Collect condition table
  df_conditions <- give_powRICLPM_conditions(x)
  df_conditions <- cbind(condition = 1:length(x$conditions), df_conditions)
  df_conditions <- idrop_DPM_reliability_column(x, df_conditions)
  condition_col_names <- c("Condition", icondition_table_names(x, df_conditions, icc_table_label))

  # Print header
  n_conditions <- length(x$conditions)
  cat(
    "powRICLPM (", as.character(ipowRICLPM_version(x)), ") simulated power for the ",
    imodel_display_name(x), " for ", n_conditions, " ",
    iexperimental_condition_label(n_conditions), ":",
    sep = ""
  )
  iprint_misspecification_warning(x)

  # Format condition table
  print(
    knitr::kable(
      df_conditions,
      format = "simple",
      align = rep("r", times = length(colnames(df_conditions))),
      col.names = condition_col_names
    )
  )
  iprint_custom_loadings(x)
  invisible(x)
}



#' Print Summary Call powRICLPM
#'
#' @noRd
print.summary.powRICLPM <- function(x, ..., object) {
  n_conditions <- nrow(x)
  cat(
    "powRICLPM (", as.character(ipowRICLPM_version(object)), ") simulated power for the ",
    imodel_display_name(object), " for ", n_conditions, " ",
    iexperimental_condition_label(n_conditions), ".",
    sep = ""
  )
  iprint_misspecification_warning(object)
  cat("\n")
  print(
    knitr::kable(
      x,
      align = rep("r", times = length(colnames(x))),
      format = "simple",
      caption = "SUMMARY OF ANALYSIS PER EXPERIMENTAL CONDITION"
    )
  )
  iprint_custom_loadings(object)
}

#' Print Summary Condition Call powRICLPM
#'
#' @noRd
print.summary.powRICLPM.condition <- function(x, ..., object) {
  iprint_misspecification_warning(object)
  cat("\n")
  print(
    knitr::kable(
      x$replications,
      align = "lr",
      format = "simple",
      caption = "SUMMARY OF ANALYSIS"
    )
  )
  cat("\n")
  print(
    knitr::kable(
      x$summary_condition,
      align = "lr",
      format = "simple",
      caption = paste0("SUMMARY OF ", imodel_display_short(object), " SIMULATION CONDITION")
    )
  )
  cat("\n")
  print(
    knitr::kable(
      x$results,
      format = "simple",
      align = "lrrrrrrrr",
      escape = FALSE,
      caption = "SIMULATION RESULTS"
    )
  )
}

#' Print Summary Parameter Call powRICLPM
#'
#' @noRd
print.summary.powRICLPM.parameter <- function(x, ..., parameter, object) {
  iprint_misspecification_warning(object)
  print(
    knitr::kable(
      x,
      format = "simple",
      align = rep("r", times = length(colnames(x))),
      caption = paste0("SIMULATION RESULTS FOR ", parameter, " (", imodel_display_short(object), ")")
    )
  )
}

iprint_misspecification_warning <- function(object) {
  if (!isTRUE(object$session$misspecified_restrictive)) {
    return(invisible(NULL))
  }
  cat(
    "\n\n",
    imisspecification_warning_text(object$session$misspecification_restrictive_reasons, past = TRUE),
    sep = ""
  )
  invisible(NULL)
}

#' Print Mplus Call powRICLPM
#'
#' @noRd
print.powRICLPM.Mplus <- function(x, ..., save_path, icc_label = "Intraclass correlation") {

  # Collect condition table
  df_conditions <- do.call(rbind, lapply(x, function(condition) {
    cbind(
      data.frame(
        sample_size = condition$sample_size,
        time_points = condition$time_points,
        ICC = condition$ICC,
        stringsAsFactors = FALSE
      ),
      icondition_reliability_columns(condition)
    )
  }))
  df_conditions <- icollapse_equal_reliability_columns(df_conditions)
  df_conditions <- cbind(condition = 1:length(x), df_conditions)
  reliability_col_names <- character()
  if ("reliability" %in% names(df_conditions)) {
    reliability_col_names <- c(reliability_col_names, "Reliability")
  }
  if ("reliability_A" %in% names(df_conditions)) {
    reliability_col_names <- c(reliability_col_names, "Reliability A")
  }
  if ("reliability_B" %in% names(df_conditions)) {
    reliability_col_names <- c(reliability_col_names, "Reliability B")
  }

  cli::cli_alert_info("Mplus input files for power analysis have been saved to {.path {save_path}}.\n
                      The conditions numbers correspond to the following conditions:")

  # Format condition table
  print(
    knitr::kable(
      df_conditions,
      format = "simple",
      align = rep("r", times = length(colnames(df_conditions))),
      col.names = c("Condition", "Sample size", "Time points", icc_label, reliability_col_names)
    )
  )
}



