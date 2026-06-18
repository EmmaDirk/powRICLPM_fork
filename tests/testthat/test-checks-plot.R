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
  expect_error(icheck_plot_parameter(c("wB2~wA1", "wA2~wB1"), out))
  expect_error(icheck_plot_parameter(12, out))
  expect_error(icheck_plot_parameter("wB4~wA3", out))
  expect_error(icheck_plot_parameter("wY2~wX1", out))

  p <- plot(out, parameter = "wB2~wA1")
  expect_s3_class(p, "ggplot")
  p_icc <- suppressMessages(plot(out, parameter = "wB2~wA1", facet_by = "ICC"))
  expect_s3_class(p_icc, "ggplot")

  p_average <- plot(out, y = "average", parameter = "wB2~wA1")
  p_empse <- plot(out, y = "EmpSE", parameter = "wB2~wA1")
  p_sd <- plot(out, y = "SD", parameter = "wB2~wA1")
  p_accuracy <- plot(out, y = "accuracy", parameter = "wB2~wA1")
  expect_s3_class(p_average, "ggplot")
  expect_s3_class(p_empse, "ggplot")
  expect_s3_class(p_sd, "ggplot")
  expect_equal(p_sd$labels$y, "EmpSE")
  expect_s3_class(p_accuracy, "ggplot")
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

test_that("icheck_plot_options() works", {
  expect_null(icheck_plot_options("time_points"))
  expect_null(icheck_plot_options("intraclass_correlation"))
  expect_null(icheck_plot_options("ICC"))
  expect_error(icheck_plot_options("sample_size"))
  expect_error(icheck_plot_options(1))
  expect_error(icheck_plot_options(c("ICC", "time_points")))
})
