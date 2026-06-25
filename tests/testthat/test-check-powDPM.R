test_that("DPM time-point checks work", {
  expect_null(icheck_T_DPM(4))
  expect_error(icheck_T_DPM("4"), "integers")
  expect_error(icheck_T_DPM(NA_real_), "finite")
  expect_error(icheck_T_DPM(2.5), "integers")
  expect_error(icheck_T_DPM(3), "4 time points")
})

test_that("DPM constraint checks reject unsupported options", {
  expect_null(icheck_DPM_constraints("none"))
  expect_null(icheck_DPM_constraints(c("lagged", "residuals")))
  expect_null(icheck_DPM_constraints("ME"))
  expect_error(icheck_DPM_constraints("within"), "within")
})
