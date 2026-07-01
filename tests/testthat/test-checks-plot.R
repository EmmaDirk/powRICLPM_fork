test_that("icheck_plot_parameter() works", {

  out <- powRICLPM(
    target_power = 0.8,
    sample_size = 1000,
    time_points = c(3, 4),
    ICC = 0.5,
    RI_cor = 0.3,
    lagged_effects = matrix(c(0.4, 0.15, 0.2, 0.3), ncol = 2, byrow = TRUE),
    within_cor = 0.3,
    reps = 2,
    seed = 123456
  )

  expect_null(icheck_plot_parameter("wB2~wA1", out))
  expect_error(icheck_plot_parameter(object = out))
  expect_error(icheck_plot_parameter(c("wB2~wA1", "wA2~wB1"), out), "length 2")
  expect_error(icheck_plot_parameter(12, out), "type double")
  expect_error(icheck_plot_parameter("wB4~wA3", out))
  expect_error(icheck_plot_parameter("wY2~wX1", out))

  p <- plot(out, parameter = "wB2~wA1")
  expect_s3_class(p, "ggplot")
  p_visible <- withVisible(plot(out, parameter = "wB2~wA1"))
  expect_false(p_visible$visible)
  expect_s3_class(p_visible$value, "ggplot")
  p_icc <- suppressMessages(plot(out, parameter = "wB2~wA1", facet_by = "ICC"))
  expect_s3_class(p_icc, "ggplot")

  p_average <- plot(out, y = "average", parameter = "wB2~wA1")
  p_empse <- plot(out, y = "EmpSE", parameter = "wB2~wA1")
  p_sd <- plot(out, y = "SD", parameter = "wB2~wA1")
  p_accuracy <- plot(out, y = "accuracy", parameter = "wB2~wA1")
  expect_error(plot(out, parameter = "wB2~wA1", color_by = "none"), "color_by = 'none'")
  expect_error(plot(out, parameter = "wB2~wA1", facet_by = "none"), "facet_by = 'none'")
  expect_error(plot(out, parameter = "wB2~wA1", facet_by = "AF_proportion"), "only available for DPM")
  expect_s3_class(p_average, "ggplot")
  expect_s3_class(p_empse, "ggplot")
  expect_s3_class(p_sd, "ggplot")
  expect_equal(p_sd$labels$y, "EmpSE")
  expect_s3_class(p_accuracy, "ggplot")

  dpm_out <- list(
    conditions = list(list(
      sample_size = 800,
      time_points = 3,
      ICC = 0.2,
      reliability = 1,
      estimates = data.frame(
        parameter = "B2~A1",
        population_value = 0.2,
        average = 0.2,
        bias = 0,
        minimum = 0.1,
        EmpSE = 0.05,
        SEAvg = 0.05,
        MSE = 0.0025,
        accuracy = 0.2,
        coverage = 0.95,
        power = 0.8
      ),
      MCSEs = data.frame(
        MCSE_average = 0.01,
        MCSE_EmpSE = 0.01,
        MCSE_SEAvg = 0.01,
        MCSE_MSE = 0.01,
        MCSE_accuracy = 0.01,
        MCSE_coverage = 0.01,
        MCSE_power = 0.01,
        MCSE_SD = 0.01
      )
    )),
    session = list(
      model = "DPM",
      target_power = 0.8,
      argument_names = list(intraclass_correlation = "AF_proportion")
    )
  )
  class(dpm_out) <- c("powRICLPM", "list")

  p_dpm <- suppressMessages(plot(dpm_out, parameter = "B2~A1"))
  expect_s3_class(p_dpm, "ggplot")
  expect_null(p_dpm$mapping$shape)
  expect_error(
    plot(dpm_out, parameter = "B2~A1", shape_by = "reliability"),
    "not available for DPM"
  )
  expect_error(
    plot(dpm_out, parameter = "B2~A1", facet_by = "ICC"),
    "not available for DPM"
  )
  expect_error(
    plot(dpm_out, parameter = "B2~A1", color_by = "intraclass_correlation"),
    "not available for DPM"
  )
})

test_that("icheck_y() works", {
  expect_null(icheck_y("power"))
  expect_null(icheck_y("bias"))
  expect_null(icheck_y("EmpSE"))
  expect_null(icheck_y("SD"))
  expect_error(icheck_y("minimum"), "minimum")
  expect_error(icheck_y("sample_size"))
  expect_error(icheck_y(3))
})

test_that("plot.powRICLPM() extracts parameter results once", {
  body_text <- paste(deparse(body(plot.powRICLPM)), collapse = "\n")
  matches <- gregexpr(
    "give_powRICLPM_results(x, parameter = parameter)",
    body_text,
    fixed = TRUE
  )[[1]]
  n_matches <- if (identical(matches, -1L)) 0L else length(matches)

  expect_equal(n_matches, 1L)
})

test_that("icheck_plot_options() works", {
  expect_null(icheck_plot_options("time_points"))
  expect_null(icheck_plot_options("intraclass_correlation"))
  expect_null(icheck_plot_options("ICC"))
  expect_null(icheck_plot_options("AF_proportion"))
  expect_null(icheck_plot_options("none"))
  expect_error(icheck_plot_options("sample_size"))
  expect_error(icheck_plot_options(1))
  expect_error(icheck_plot_options(c("ICC", "time_points")))
})
