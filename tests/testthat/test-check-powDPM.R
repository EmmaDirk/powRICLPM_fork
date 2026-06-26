test_that("DPM time-point checks work", {
  expect_null(icheck_T_DPM(3))
  expect_null(icheck_T_DPM(4))
  expect_error(icheck_T_DPM("4"), "integers")
  expect_error(icheck_T_DPM(NA_real_), "finite")
  expect_error(icheck_T_DPM(2.5), "integers")
  expect_error(icheck_T_DPM(2), "at least 3 time points")
})

test_that("DPM constraint checks reject unsupported options", {
  expect_null(icheck_DPM_constraints("none"))
  expect_null(icheck_DPM_constraints(c("lagged", "residuals")))
  expect_null(icheck_DPM_constraints("AF_loadings_free"))
  expect_null(icheck_DPM_constraints("ME", estimate_ME = TRUE))
  expect_error(icheck_DPM_constraints("ME"), "estimate_ME = TRUE")
  expect_error(icheck_DPM_constraints("within"), "invalid DPM constraints")
  expect_error(icheck_DPM_constraints("loadings_free"), "invalid DPM constraints")
  expect_error(icheck_DPM_constraints("RI_loadings_free"), "invalid DPM constraints")
  expect_error(icheck_DPM_constraints(c("stationarity", "lagged")), "cannot combine")
})

test_that("DPM identification table is enforced", {
  for (i in seq_len(nrow(DPM_identification_table))) {
    row <- DPM_identification_table[i, ]
    constraints <- strsplit(row$constraints_key, "+", fixed = TRUE)[[1]]
    constraints <- constraints[constraints != "none"]
    if (isTRUE(row$AF_loadings_free)) {
      constraints <- c(constraints, "AF_loadings_free")
    }
    if (length(constraints) == 0L) {
      constraints <- "none"
    }
    if (row$min_waves > 3L) {
      expect_error(
        icheck_DPM_identification(
          time_points = row$min_waves - 1L,
          reliability = if (row$estimate_ME) .8 else 1,
          estimate_ME = row$estimate_ME,
          constraints = constraints
        ),
        "not identified"
      )
    }
    expect_null(icheck_DPM_identification(
      time_points = row$min_waves,
      reliability = if (row$estimate_ME) .8 else 1,
      estimate_ME = row$estimate_ME,
      constraints = constraints
    ))
  }
})

test_that("DPM misspecification detection separates restrictive and general cases", {
  restrictive <- detect_DPM_misspecification(
    reliability = .8,
    loadings = NULL,
    time_points = 4,
    constraints = "none",
    estimate_ME = FALSE
  )
  expect_true(restrictive$restrictive)
  expect_match(restrictive$restrictive_reasons, "measurement error")

  general <- detect_DPM_misspecification(
    reliability = 1,
    loadings = NULL,
    time_points = 4,
    constraints = "none",
    estimate_ME = TRUE
  )
  expect_false(general$restrictive)
  expect_true(general$general)
})
