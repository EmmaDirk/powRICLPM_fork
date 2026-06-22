library(pkgload)

.libPaths(c("C:/Users/Admin/AppData/Local/R/win-library/4.6", .libPaths()))
pkgload::load_all(".", quiet = TRUE)

results <- data.frame(
  name = character(),
  status = character(),
  detail = character(),
  stringsAsFactors = FALSE
)

add_result <- function(name, status, detail = "") {
  results <<- rbind(
    results,
    data.frame(name = name, status = status, detail = detail, stringsAsFactors = FALSE)
  )
  cat(sprintf("[%s] %s%s\n", status, name, if (nzchar(detail)) paste0(" -- ", detail) else ""))
}

check_ok <- function(name, expr, inspect = NULL) {
  value <- tryCatch(
    eval.parent(substitute(expr)),
    error = function(e) structure(list(message = conditionMessage(e)), class = "sweep_error")
  )
  if (inherits(value, "sweep_error")) {
    add_result(name, "FAIL", value$message)
    return(invisible(NULL))
  }
  if (!is.null(inspect)) {
    inspection <- tryCatch(
      inspect(value),
      error = function(e) structure(list(message = conditionMessage(e)), class = "sweep_error")
    )
    if (inherits(inspection, "sweep_error")) {
      add_result(name, "FAIL", inspection$message)
      return(invisible(NULL))
    }
    if (!isTRUE(inspection)) {
      add_result(name, "FAIL", paste("inspection returned", paste(inspection, collapse = ", ")))
      return(invisible(NULL))
    }
  }
  add_result(name, "PASS")
  invisible(value)
}

check_error <- function(name, expr, pattern) {
  err <- tryCatch(
    {
      eval.parent(substitute(expr))
      NULL
    },
    error = function(e) conditionMessage(e)
  )
  if (is.null(err)) {
    add_result(name, "FAIL", "expected an error but call succeeded")
    return(invisible(NULL))
  }
  if (!grepl(pattern, err, ignore.case = FALSE)) {
    add_result(name, "FAIL", paste0("error did not match /", pattern, "/: ", err))
    return(invisible(NULL))
  }
  add_result(name, "PASS")
  invisible(err)
}

check_error_no_warning <- function(name, expr, pattern) {
  warnings <- character()
  err <- withCallingHandlers(
    tryCatch(
      {
        eval.parent(substitute(expr))
        NULL
      },
      error = function(e) conditionMessage(e)
    ),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  if (is.null(err)) {
    add_result(name, "FAIL", "expected an error but call succeeded")
    return(invisible(NULL))
  }
  if (!grepl(pattern, err, ignore.case = FALSE)) {
    add_result(name, "FAIL", paste0("error did not match /", pattern, "/: ", err))
    return(invisible(NULL))
  }
  if (length(warnings) > 0) {
    add_result(name, "FAIL", paste("unexpected warning(s):", paste(warnings, collapse = " | ")))
    return(invisible(NULL))
  }
  add_result(name, "PASS")
  invisible(err)
}

contains <- function(x, pattern) any(grepl(pattern, x))
not_contains <- function(x, pattern) !any(grepl(pattern, x))

lagged <- matrix(c(0.35, 0.10, 0.12, 0.30), nrow = 2, byrow = TRUE)

riclpm <- check_ok(
  "RICLPM constructor with ICC alias",
  powRICLPM(
    target_power = 0.8,
    sample_size = 600,
    time_points = 3,
    ICC = 0.4,
    RI_cor = 0.25,
    lagged_effects = lagged,
    within_cor = 0.2,
    reps = 1,
    seed = 20260622
  ),
  function(x) inherits(x, "powRICLPM") && identical(x$session$model, "RICLPM")
)

riclpm_long <- check_ok(
  "RICLPM constructor with intraclass_correlation and lagged/residual constraints",
  powRICLPM(
    target_power = 0.8,
    sample_size = 800,
    time_points = 4,
    intraclass_correlation = 0.35,
    RI_cor = 0.20,
    lagged_effects = lagged,
    within_cor = 0.15,
    reps = 1,
    seed = 20260623,
    constraints = c("lagged", "residuals")
  ),
  function(x) inherits(x, "powRICLPM") && identical(x$session$constraints, c("lagged", "residuals"))
)

starts <- check_ok(
  "STARTS constructor with estimated measurement error",
  powRICLPM(
    target_power = 0.8,
    sample_size = 1000,
    time_points = 4,
    ICC = 0.4,
    RI_cor = 0.25,
    lagged_effects = lagged,
    within_cor = 0.2,
    reliability = 0.85,
    estimate_ME = TRUE,
    reps = 1,
    seed = 20260624
  ),
  function(x) inherits(x, "powRICLPM") && isTRUE(x$session$estimate_ME)
)

dpm <- check_ok(
  "DPM constructor with AF names",
  powRICLPM(
    model = "DPM",
    target_power = 0.8,
    sample_size = 600,
    time_points = 3,
    AF_proportion = 0.25,
    AF_cor = 0.2,
    lagged_effects = lagged,
    wave_cor = 0.1,
    reps = 1,
    seed = 20260625
  ),
  function(x) inherits(x, "powRICLPM") && identical(x$session$model, "DPM")
)

dpm_free <- check_ok(
  "DPM constructor with free accumulating-factor loadings",
  powRICLPM(
    model = "DPM",
    target_power = 0.8,
    sample_size = 800,
    time_points = 4,
    AF_proportion = 0.25,
    AF_cor = 0.2,
    lagged_effects = lagged,
    wave_cor = 0.1,
    constraints = "loadings_free",
    loadings = c(NA, 1, 0.85, 0.75),
    reps = 1,
    seed = 20260626
  ),
  function(x) inherits(x, "powRICLPM") && identical(x$session$model, "DPM")
)

dpm_two_wave <- check_ok(
  "DPM constructor allows two waves",
  powRICLPM(
    model = "DPM",
    target_power = 0.8,
    sample_size = 600,
    time_points = 2,
    AF_proportion = 0.25,
    AF_cor = 0.2,
    lagged_effects = lagged,
    wave_cor = 0.1,
    reps = 1,
    seed = 20260627
  ),
  function(x) inherits(x, "powRICLPM") && identical(x$conditions[[1]]$time_points, 2)
)

dpm_constrained <- check_ok(
  "DPM constructor accepts lagged plus residual constraints",
  powRICLPM(
    model = "DPM",
    target_power = 0.8,
    sample_size = 800,
    time_points = 3,
    AF_proportion = 0.25,
    AF_cor = 0.2,
    lagged_effects = lagged,
    wave_cor = 0.1,
    constraints = c("lagged", "residuals"),
    reps = 1,
    seed = 20260628
  ),
  function(x) inherits(x, "powRICLPM") && setequal(x$session$constraints, c("lagged", "residuals"))
)

riclpm_time_varying_rel <- check_ok(
  "RICLPM constructor accepts vector time-varying reliability",
  powRICLPM(
    target_power = 0.8,
    sample_size = 800,
    time_points = 4,
    ICC = 0.4,
    RI_cor = 0.25,
    lagged_effects = lagged,
    within_cor = 0.2,
    reliability = c(0.8, 0.85, 0.9, 1),
    reps = 1,
    seed = 20260629
  ),
  function(x) inherits(x, "powRICLPM") && identical(x$conditions[[1]]$reliability, "time-varying")
)

riclpm_search <- check_ok(
  "RICLPM constructor accepts sample-size search",
  powRICLPM(
    target_power = 0.8,
    search_lower = 600,
    search_upper = 640,
    search_step = 40,
    time_points = 3,
    ICC = 0.4,
    RI_cor = 0.25,
    lagged_effects = lagged,
    within_cor = 0.2,
    reps = 1,
    seed = 20260630
  ),
  function(x) inherits(x, "powRICLPM") && length(x$conditions) == 2L
)

check_error(
  "DPM rejects missing AF_proportion",
  powRICLPM(
    model = "DPM",
    target_power = 0.8,
    sample_size = 600,
    time_points = 3,
    AF_cor = 0.2,
    lagged_effects = lagged,
    wave_cor = 0.1,
    reps = 1
  ),
  "AF_proportion.*must be specified"
)

check_error_no_warning(
  "DPM prioritizes ICC alias before missing AF_proportion",
  powRICLPM(
    model = "DPM",
    target_power = 0.8,
    sample_size = 600,
    time_points = 3,
    ICC = 0.25,
    AF_cor = 0.2,
    lagged_effects = lagged,
    wave_cor = 0.1,
    reps = 1
  ),
  "model = 'RICLPM'"
)

check_error(
  "DPM rejects RI_cor",
  powRICLPM(
    model = "DPM",
    target_power = 0.8,
    sample_size = 600,
    time_points = 3,
    AF_proportion = 0.25,
    RI_cor = 0.2,
    lagged_effects = lagged,
    wave_cor = 0.1,
    reps = 1
  ),
  "RI_cor.*not valid.*AF_cor"
)

check_error(
  "DPM rejects within_cor",
  powRICLPM(
    model = "DPM",
    target_power = 0.8,
    sample_size = 600,
    time_points = 3,
    AF_proportion = 0.25,
    AF_cor = 0.2,
    lagged_effects = lagged,
    within_cor = 0.1,
    reps = 1
  ),
  "within_cor.*not valid.*wave_cor"
)

check_error(
  "DPM rejects reliability below 1",
  powRICLPM(
    model = "DPM",
    target_power = 0.8,
    sample_size = 600,
    time_points = 3,
    AF_proportion = 0.25,
    AF_cor = 0.2,
    lagged_effects = lagged,
    wave_cor = 0.1,
    reliability = 0.8,
    reps = 1
  ),
  "does not separate measurement error"
)

check_error(
  "DPM rejects reliability vector",
  powRICLPM(
    model = "DPM",
    target_power = 0.8,
    sample_size = 600,
    time_points = 3,
    AF_proportion = 0.25,
    AF_cor = 0.2,
    lagged_effects = lagged,
    wave_cor = 0.1,
    reliability = c(1, 1, 1),
    reps = 1
  ),
  "reliability.*DPM"
)

check_error(
  "DPM rejects measurement-error estimation",
  powRICLPM(
    model = "DPM",
    target_power = 0.8,
    sample_size = 600,
    time_points = 3,
    AF_proportion = 0.25,
    AF_cor = 0.2,
    lagged_effects = lagged,
    wave_cor = 0.1,
    estimate_ME = TRUE,
    reps = 1
  ),
  "does not separate measurement error"
)

check_error_no_warning(
  "DPM invalid constraints error does not emit seed warning first",
  powRICLPM(
    model = "DPM",
    target_power = 0.8,
    sample_size = 600,
    time_points = 3,
    AF_proportion = 0.25,
    AF_cor = 0.2,
    lagged_effects = lagged,
    wave_cor = 0.1,
    constraints = "ME",
    reps = 1
  ),
  "does not separate measurement error"
)

check_error(
  "DPM rejects within shorthand constraints",
  powRICLPM(
    model = "DPM",
    target_power = 0.8,
    sample_size = 600,
    time_points = 3,
    AF_proportion = 0.25,
    AF_cor = 0.2,
    lagged_effects = lagged,
    wave_cor = 0.1,
    constraints = "within",
    reps = 1
  ),
  "constraints = 'within'"
)

check_error(
  "DPM rejects ME constraints",
  powRICLPM(
    model = "DPM",
    target_power = 0.8,
    sample_size = 600,
    time_points = 3,
    AF_proportion = 0.25,
    AF_cor = 0.2,
    lagged_effects = lagged,
    wave_cor = 0.1,
    constraints = "ME",
    reps = 1
  ),
  "does not separate measurement error"
)

check_error(
  "DPM rejects Mplus software",
  powRICLPM(
    model = "DPM",
    target_power = 0.8,
    sample_size = 600,
    time_points = 3,
    AF_proportion = 0.25,
    AF_cor = 0.2,
    lagged_effects = lagged,
    wave_cor = 0.1,
    software = "Mplus",
    reps = 1
  ),
  "only available with `software = 'lavaan'`"
)

check_error(
  "DPM rejects bounded estimation",
  powRICLPM(
    model = "DPM",
    target_power = 0.8,
    sample_size = 600,
    time_points = 3,
    AF_proportion = 0.25,
    AF_cor = 0.2,
    lagged_effects = lagged,
    wave_cor = 0.1,
    bounds = TRUE,
    reps = 1
  ),
  "Bounded estimation is not yet available for the DPM"
)

check_error(
  "DPM rejects one wave",
  powRICLPM(
    model = "DPM",
    target_power = 0.8,
    sample_size = 600,
    time_points = 1,
    AF_proportion = 0.25,
    AF_cor = 0.2,
    lagged_effects = lagged,
    wave_cor = 0.1,
    reps = 1
  ),
  "at least 2 time points"
)

check_error(
  "RICLPM rejects two waves",
  powRICLPM(
    target_power = 0.8,
    sample_size = 600,
    time_points = 2,
    ICC = 0.4,
    RI_cor = 0.25,
    lagged_effects = lagged,
    within_cor = 0.2,
    reps = 1
  ),
  "RI-CLPM is not identified"
)

check_error(
  "STARTS rejects three waves",
  powRICLPM(
    target_power = 0.8,
    sample_size = 600,
    time_points = 3,
    ICC = 0.4,
    RI_cor = 0.25,
    lagged_effects = lagged,
    within_cor = 0.2,
    estimate_ME = TRUE,
    reps = 1
  ),
  "larger than 3"
)

check_error(
  "RICLPM rejects DPM names without model = DPM",
  powRICLPM(
    target_power = 0.8,
    sample_size = 600,
    time_points = 3,
    AF_proportion = 0.25,
    AF_cor = 0.2,
    lagged_effects = lagged,
    wave_cor = 0.1,
    reps = 1
  ),
  "model.*defaults to.*RICLPM.*model = 'DPM'"
)

check_error(
  "RICLPM rejects AF_proportion in summary selector",
  summary(riclpm, sample_size = 600, time_points = 3, AF_proportion = 0.4),
  "only available for DPM"
)

check_error(
  "DPM rejects ICC in summary selector",
  summary(dpm, sample_size = 600, time_points = 3, ICC = 0.25),
  "not available for DPM"
)

check_error(
  "DPM rejects reliability in summary selector",
  summary(dpm, reliability = 1),
  "not available for DPM"
)

check_error(
  "RICLPM rejects AF_proportion in give selector",
  give(riclpm, "AF_proportion"),
  "only available for DPM"
)

check_error(
  "DPM rejects ICC in give selector",
  give(dpm, "ICC"),
  "not available for DPM"
)

check_error(
  "DPM rejects reliability in give selector",
  give(dpm, "reliability"),
  "not available for DPM"
)

check_error(
  "give results requires parameter",
  give(riclpm, "results"),
  "No .*parameter.*specified"
)

check_error(
  "give results rejects unknown parameter",
  give(riclpm, "results", parameter = "not_a_parameter"),
  "not found"
)

check_ok(
  "RICLPM give conditions includes reliability",
  give(riclpm, "conditions"),
  function(x) "reliability" %in% names(x) && "ICC" %in% names(x)
)

check_ok(
  "DPM give conditions includes AF_proportion and omits reliability",
  give(dpm, "conditions"),
  function(x) "AF_proportion" %in% names(x) && !("reliability" %in% names(x))
)

check_ok(
  "RICLPM summary parameter labels ICC and Reliability",
  summary(riclpm, parameter = "wB2~wA1"),
  function(x) all(c("ICC", "Reliability") %in% names(x))
)

check_ok(
  "DPM summary parameter labels AF proportion and omits Reliability",
  summary(dpm, parameter = "B2~A1"),
  function(x) "AF proportion" %in% names(x) && !("Reliability" %in% names(x))
)

check_ok(
  "DPM summary condition accepts AF_proportion",
  summary(dpm, sample_size = 600, time_points = 3, AF_proportion = 0.25),
  function(x) all(c("Population", "Avg", "Power") %in% names(x))
)

check_error(
  "summary rejects unknown parameter",
  summary(riclpm, parameter = "not_a_parameter"),
  "parameter.*not found"
)

check_ok(
  "summary infers single time-varying reliability condition",
  summary(riclpm_time_varying_rel, sample_size = 800, time_points = 4, ICC = 0.4),
  function(x) all(c("Population", "Avg", "Power") %in% names(x))
)

check_ok(
  "summary accepts explicit time-varying reliability label",
  summary(riclpm_time_varying_rel, sample_size = 800, time_points = 4, ICC = 0.4, reliability = "time-varying"),
  function(x) all(c("Population", "Avg", "Power") %in% names(x))
)

check_error(
  "RICLPM rejects AF_proportion plot mapping",
  plot(riclpm, parameter = "wB2~wA1", facet_by = "AF_proportion"),
  "only available for DPM"
)

check_error(
  "DPM rejects ICC plot mapping",
  plot(dpm, parameter = "B2~A1", facet_by = "ICC"),
  "not available for DPM"
)

check_error(
  "DPM rejects reliability plot mapping",
  plot(dpm, parameter = "B2~A1", shape_by = "reliability"),
  "not available for DPM"
)

check_error(
  "plot rejects missing parameter",
  plot(riclpm),
  "No .*parameter"
)

check_error(
  "plot rejects unsupported y",
  plot(riclpm, parameter = "wB2~wA1", y = "minimum"),
  "minimum"
)

check_ok(
  "plot accepts SD alias as EmpSE",
  suppressMessages(plot(riclpm, parameter = "wB2~wA1", y = "SD")),
  function(x) inherits(x, "ggplot") && identical(x$labels$y, "EmpSE")
)

check_ok(
  "DPM plot accepts AF_proportion facet",
  suppressMessages(plot(dpm, parameter = "B2~A1", facet_by = "AF_proportion")),
  function(x) inherits(x, "ggplot") && is.null(x$mapping$shape)
)

check_ok(
  "DPM plot default has no shape mapping",
  suppressMessages(plot(dpm, parameter = "B2~A1")),
  function(x) inherits(x, "ggplot") && is.null(x$mapping$shape)
)

check_ok(
  "RICLPM print names model and reliability",
  capture.output(print(riclpm)),
  function(x) contains(x, "RI-CLPM") && contains(x, "Reliability")
)

check_ok(
  "STARTS print names STARTS model",
  capture.output(print(starts)),
  function(x) contains(x, "STARTS model")
)

check_ok(
  "DPM print names model and omits reliability",
  capture.output(print(dpm)),
  function(x) contains(x, "Dynamic Panel Model \\(DPM\\)") && not_contains(x, "Reliability")
)

check_ok(
  "check_loadings interprets RI-CLPM vector",
  capture.output(check_loadings(c(1, 0.9, 1.1), time_points = 3)),
  function(x) contains(x, "random-intercept loadings")
)

check_ok(
  "check_loadings interprets DPM vector",
  capture.output(check_loadings(c(NA, 1, 0.8), time_points = 3, model = "DPM")),
  function(x) contains(x, "accumulating-factor loadings")
)

check_error(
  "check_loadings rejects RI-CLPM NA with DPM guidance",
  check_loadings(c(NA, 1, 0.8), time_points = 3),
  "only valid for DPM loadings"
)

check_error(
  "check_loadings rejects DPM missing first-wave NA",
  check_loadings(c(1, 0.8, 0.7), time_points = 3, model = "DPM"),
  "first DPM loading must be `NA`"
)

check_error(
  "check_loadings rejects DPM second loading not one",
  check_loadings(c(NA, 0.8, 0.7), time_points = 3, model = "DPM"),
  "second DPM loading must be 1"
)

check_ok(
  "check_reliability interprets vector reliability",
  capture.output(check_reliability(c(0.8, 0.85, 1), time_points = 3)),
  function(x) contains(x, "A2 and B2 have reliability 0.85")
)

check_error(
  "check_reliability rejects DPM-style NA",
  check_reliability(c(NA, 1, 0.8), time_points = 3),
  "finite"
)

failed <- results[results$status != "PASS", , drop = FALSE]
cat("\nSweep summary:\n")
print(table(results$status))

if (nrow(failed) > 0) {
  cat("\nFailures:\n")
  print(failed, row.names = FALSE)
  quit(status = 1)
}

cat("\nAll broad interface sweep checks passed.\n")
