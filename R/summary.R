#' Summarize Results from \code{powRICLPM} Object
#'
#' @description
#' S3 method for class \code{powRICLPM}. \code{summary.powRICLPM} summarizes the setup and results of the \code{powRICLPM} analysis. Depending on the arguments that are set, \code{summary.powRICLPM} provides a different summary (see "Details").
#'
#' @param object A \code{powRICLPM} object.
#' @param ... (don't use)
#' @param parameter Character string of length 1 denoting the parameter to visualize the results for.
#' @param sample_size (optional) An \code{integer}, denoting the sample size of the experimental condition of interest.
#' @param time_points (optional) An \code{integer}, denoting the number of time points of the experimental condition of interest.
#' @param intraclass_correlation (optional) A \code{double}, denoting the proportion of variance at the between-unit level of the experimental condition of interest.
#' @param AF_proportion (optional) A \code{double}, denoting the accumulating-factor variance proportion of the DPM experimental condition of interest.
#' @param reliability (optional) A \code{numeric} or \code{character} value denoting the scalar reliability or reliability label of the experimental condition of interest. This is not available for DPM objects.
#' @param ICC Alternative name for \code{intraclass_correlation}.
#'
#' @return No return value, called for side effects.
#'
#' @details
#' \code{summary.powRICLPM} provides a different summary of the \code{powRICLPM} object, depending on the additional arguments that are set:
#' \itemize{
#'   \item When \code{sample_size = ...}, \code{time_points = ...}, and \code{intraclass_correlation = ...} or \code{AF_proportion = ...} are set: Estimation information and results for all parameters in that experimental condition. If multiple RI-CLPM conditions differ only by \code{reliability}, specify \code{reliability = ...} as well.
#'   \item When \code{parameter = "..."} is set: Estimation information and results for a specific parameter across all experimental conditions.
#'   \item No additional arguments: Characteristics of the different experimental conditions are summarized, as well as session info (information that applies to all conditions, such the number of replications, etc.).
#' }
#' \subsection{Interpretation Output}{Depending on the arguments that you set, \code{summary()} prints a table with different analysis outcomes in the columns and where each row refers to a different experimental condition. The following information is available:
#'  \itemize{
#'   \item \code{Sample size}, \code{Time points}, \code{ICC} or \code{AF proportion}, and \code{Reliability} when applicable: The experimental condition that the row refers to.
#'   \item \code{Population}: The true value of the parameter.
#'   \item \code{Avg}: The average (across replications) parameter estimate.
#'   \item \code{Bias}: The difference between the population value and the average parameter estimate.
#'   \item \code{Min}: The lowest (across replications) parameter estimate.
#'   \item \code{EmpSE}: The empirical standard error, computed as the standard deviation of the parameter estimate over replications.
#'   \item \code{SEAvg}: The average (across replications) standard error of the parameter estimate.
#'   \item \code{MSE}: The parameter mean square error, combining a parameter's bias and efficiency.
#'   \item \code{Accuracy}: The average (across replications) width of the confidence interval.
#'   \item \code{Cover}: The coverage rate, representing the proportion of times (across replications) the true parameter estimate fell in the confidence interval.
#'   \item \code{Power}: The proportion of times (across replications) the confidence interval did not contain zero.
#'   \item \code{Error}: The number of replications that failed to run (i.e., \code{lavaan()} produced an error).
#'   \item \code{Not converged}: The number of replications that did not converge to a solution.
#'   \item \code{Inadmissible}: The number of replications that converged to an inadmissible solution (e.g., a variance estimated to be lower than zero).
#'   }
#' }
#'
#' @examples
#' \dontshow{load(system.file("extdata", "out_preliminary.RData", package = "powRICLPM"))}
#' # Get setup of powRICLPM analysis and convergence issues
#' summary(out_preliminary)
#'
#' # Performance measures for "wB2~wA1" parameter across experimental conditions
#' summary(out_preliminary, parameter = "wB2~wA1")
#'
#' # Performance measures for all parameters, for specific experimental condition
#' summary(out_preliminary, sample_size = 700, time_points = 4, intraclass_correlation = .3, reliability = 1)
#'
#' @method summary powRICLPM
#' @export
summary.powRICLPM <- function(
    object,
    ...,
    parameter = NULL,
    sample_size = NULL,
    time_points = NULL,
    intraclass_correlation = NULL,
    AF_proportion = NULL,
    reliability = NULL,
    ICC = NULL
  ) {

  call_summary <- match.call()
  icheck_object_summary(object)
  dpm_object <- iis_DPM_object(object)
  riclpm_icc_supplied <- !is.null(intraclass_correlation) || !is.null(ICC)
  if (!is.null(ICC) && !is.null(intraclass_correlation)) {
    iabort_renamed_argument_conflict("intraclass_correlation", "ICC")
  }
  if (!is.null(AF_proportion) &&
      riclpm_icc_supplied) {
    iabort_renamed_argument_conflict("AF_proportion", "intraclass_correlation/ICC")
  }
  if (dpm_object && riclpm_icc_supplied) {
    cli::cli_abort(
      c(
        "`intraclass_correlation` and `ICC` are not available for DPM summaries:",
        i = "Use `AF_proportion` to select DPM conditions."
      )
    )
  }
  if (!dpm_object && !is.null(AF_proportion)) {
    cli::cli_abort(
      c(
        "`AF_proportion` is only available for DPM summaries:",
        i = "Use `intraclass_correlation` or `ICC` to select RI-CLPM and STARTS conditions."
      )
    )
  }
  if (!is.null(AF_proportion)) {
    intraclass_correlation <- AF_proportion
  }
  if (!is.null(ICC)) {
    intraclass_correlation <- ICC
  }
  ICC <- intraclass_correlation
  icc_table_label <- iicc_table_name(object = object, call = call_summary)

  # Argument validation
  if (!is.null(parameter)) {icheck_parameter_summary(parameter, object)}
  if (!is.null(sample_size)) {icheck_sample_size_summary(sample_size, object)}
  if (!is.null(time_points)) {icheck_time_points_summary(time_points, object)}
  if (!is.null(ICC)) {icheck_ICC_summary(ICC, object)}
  if (iis_DPM_object(object) && !is.null(reliability)) {
    cli::cli_abort(
      c(
        "`reliability` is not available for DPM summaries:",
        i = "The DPM does not separate measurement error, so reliability is not part of the DPM simulation conditions.",
        i = "Select DPM conditions with `sample_size`, `time_points`, and `AF_proportion`."
      )
    )
  }
  if (!is.null(reliability)) {icheck_reliability_summary(reliability, object)}


  # Summarize
  if (!is.null(sample_size) && !is.null(time_points) && !is.null(ICC)) {

    # Collect information for print.summary.powRICLPM.condition()
    condition <- imatch_condition_summary(
      object = object,
      sample_size = sample_size,
      time_points = time_points,
      ICC = ICC,
      reliability = reliability
    )

    ## Simulation results
    results <- condition$estimates[, -1]
    results <- round(results, digits = 3)
    colnames(results) <- c("Population", "Avg", "Bias", "Min", "EmpSE", "SEAvg", "MSE", "Accuracy", "Cover", "Power")
    rownames(results) <- condition$estimates$parameter

    ## Summary of analysis
    replication_values <- c(
      object$session$reps,
      condition$estimation_information$n_completed
    )
    replication_names <- c("Requested:", "Completed:")
    if (condition$estimation_information$n_error > 0) {
      replication_values <- c(replication_values, condition$estimation_information$n_error)
      replication_names <- c(replication_names, "Errors:")
    }
    replication_values <- c(
      replication_values,
      condition$estimation_information$n_nonconvergence,
      condition$estimation_information$n_inadmissible
    )
    replication_names <- c(
      replication_names,
      "Convergence issues:",
      "Inadmissible results:"
    )
    summary_replications <- matrix(replication_values, ncol = 1)
    colnames(summary_replications) <- c("Number of replications")
    rownames(summary_replications) <- replication_names

    ## Summary of condition
    summary_condition <- matrix(c(
      imodel_display_name(object),
      condition$skewness,
      condition$kurtosis,
      format_constraints(object$session$constraints),
      object$session$bounds,
      object$session$estimate_ME,
      condition$significance_criterion
    ), ncol = 1)
    colnames(summary_condition) <- c("Value")
    rownames(summary_condition) <- c("Model:", "Skewness:", "Kurtosis:", "Constraints:", "Bounds:", "Estimated measurement error:", "Significance criterion:")

    summary_list <- list(
      summary_condition = summary_condition,
      results = results,
      replications = summary_replications
    )

    inote_condition_loading_anchor(object, results)
    print.summary.powRICLPM.condition(summary_list, object = object)
    invisible(results)

  } else if (!is.null(parameter)) {

    # Collect information for print.summary.powRICLPM.parameter()
    parameter_df <- give_powRICLPM_results(object, parameter)
    replications_df <- give_powRICLPM_estimation_problems(object)
    parameter_summary <- merge(parameter_df, replications_df, by = c("sample_size","time_points", "ICC", "reliability"))
    parameter_summary <- idrop_DPM_reliability_column(object, parameter_summary)
    parameter_col_names <- c("Sample size", "Time points", icc_table_label)
    if (!iis_DPM_object(object)) {
      parameter_col_names <- c(parameter_col_names, "Reliability")
    }
    colnames(parameter_summary) <- c(parameter_col_names, "Population", "Avg","Bias", "Min", "EmpSE", "SEAvg", "MSE", "Accuracy", "Cover", "Power", "Error", "Not converged", "Inadmissible")
    inote_condition_loading_anchor(object, parameter_summary)
    print.summary.powRICLPM.parameter(parameter_summary, parameter = parameter, object = object)
    invisible(parameter_summary)

  } else {

    # Collect information for print.summary.powRICLPM()
    ## Summary of analysis
    replications_df <- give_powRICLPM_estimation_problems(object)
    replications_df <- idrop_DPM_reliability_column(object, replications_df)
    replications_col_names <- c("Sample size", "Time points", icc_table_label)
    if (!iis_DPM_object(object)) {
      replications_col_names <- c(replications_col_names, "Reliability")
    }
    colnames(replications_df) <- c(replications_col_names, "Error", "Not converged", "Inadmissible")
    inote_condition_loading_anchor(object, replications_df)
    print.summary.powRICLPM(replications_df, object = object)
    iprint_reliability_tables(object$conditions)
  }
}





