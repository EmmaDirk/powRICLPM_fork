make_lavaan_condition <- function(constraints, estimate_ME = FALSE, loadings = NULL) {
  lagged_effects <- matrix(c(0.4, 0.15, 0.2, 0.3), ncol = 2, byrow = TRUE)

  create_conditions(
    target_power = 0.8,
    sample_size = 300,
    time_points = 3,
    intraclass_correlation = 0.5,
    RI_cor = 0.3,
    lagged_effects = lagged_effects,
    within_cor = 0.3,
    Psi = compute_Psi(lagged_effects, 0.3),
    reliability = 1,
    loadings = loadings,
    skewness = 0,
    kurtosis = 0,
    estimate_ME = estimate_ME,
    significance_criterion = 0.05,
    reps = 2,
    bootstrap_reps = NULL,
    seed = 123456,
    constraints = constraints,
    bounds = FALSE,
    estimator = NA,
    save_path = NULL,
    software = "lavaan"
  )[[1]]
}

expect_no_duplicate_lavaan_elements <- function(condition) {
  model_rows <- condition$est_tab$op %in% c("=~", "~", "~~")
  model_elements <- paste(
    condition$est_tab$lhs[model_rows],
    condition$est_tab$op[model_rows],
    condition$est_tab$rhs[model_rows]
  )

  expect_equal(anyDuplicated(model_elements), 0)
}

test_that("constrained lavaan syntax combines labels and starting values", {
  condition_lagged <- make_lavaan_condition("lagged")
  expect_true(grepl("wA2~a*start(0.4)*wA1", condition_lagged$est_synt, fixed = TRUE))
  expect_true(grepl("wB3~d*start(0.3)*wB2", condition_lagged$est_synt, fixed = TRUE))
  expect_false(grepl("wA2~start(0.4)*wA1", condition_lagged$est_synt, fixed = TRUE))
  expect_no_duplicate_lavaan_elements(condition_lagged)

  condition_residuals <- make_lavaan_condition("residuals")
  expect_true(grepl("wA2~~rvarA*start(0.7815)*wA2", condition_residuals$est_synt, fixed = TRUE))
  expect_true(grepl("wA2~~rcov*start(0.13)*wB2", condition_residuals$est_synt, fixed = TRUE))
  expect_no_duplicate_lavaan_elements(condition_residuals)

  condition_stationarity <- make_lavaan_condition("stationarity")
  expect_true(grepl("wA1=~NA*start(1)*A1", condition_stationarity$est_synt, fixed = TRUE))
  expect_true(grepl("wA1~~cor1*start(0.3)*wB1", condition_stationarity$est_synt, fixed = TRUE))
  expect_true(grepl("wA2~~rvarA2*start(0.7815)*wA2", condition_stationarity$est_synt, fixed = TRUE))
  expect_true(grepl("wA2~~rcov2*start(0.13)*wB2", condition_stationarity$est_synt, fixed = TRUE))
  expect_no_duplicate_lavaan_elements(condition_stationarity)

  condition_ME <- make_lavaan_condition("ME", estimate_ME = TRUE)
  expect_true(grepl("A1~~MEvarA*start(0)*A1", condition_ME$est_synt, fixed = TRUE))
  expect_true(grepl("B3~~MEvarB*start(0)*B3", condition_ME$est_synt, fixed = TRUE))
  expect_no_duplicate_lavaan_elements(condition_ME)
})

test_that("combined RI_loadings_free constraints avoid duplicate lavaan syntax", {
  loadings <- matrix(c(1, 0.5, 1.2, 1, 1.5, 0.8), nrow = 2, byrow = TRUE)
  constraints <- list(
    c("lagged", "RI_loadings_free"),
    c("within", "RI_loadings_free"),
    c("stationarity", "RI_loadings_free")
  )

  for (constraint in constraints) {
    condition <- make_lavaan_condition(constraint, loadings = loadings)
    expect_true(grepl("RI_A=~lx2*start(0.5)*A2", condition$est_synt, fixed = TRUE))
    expect_true(grepl("RI_B=~ly3*start(0.8)*B3", condition$est_synt, fixed = TRUE))
    expect_no_duplicate_lavaan_elements(condition)
  }
})

test_that("constraints in powRICLPM() work", {

  lagged_effects <- matrix(c(0.4, 0.15, 0.2, 0.3), ncol = 2, byrow = TRUE)

  out_lagged <- powRICLPM(
    target_power = 0.8,
    sample_size = 300,
    time_points = 3,
    ICC = 0.5,
    RI_cor = 0.3,
    lagged_effects = lagged_effects,
    within_cor = 0.3,
    reps = 2,
    seed = 123456,
    constraints = "lagged"
  )

  # Cross-lagged effects
  expect_equal(
    out_lagged$conditions[[1]]$estimates$average[which(out_lagged$conditions[[1]]$estimates$parameter == "wB2~wA1")],
    out_lagged$conditions[[1]]$estimates$average[which(out_lagged$conditions[[1]]$estimates$parameter == "wB3~wA2")]
  )

  # Autoregressive effects
  expect_equal(
    out_lagged$conditions[[1]]$estimates$average[which(out_lagged$conditions[[1]]$estimates$parameter == "wA2~wA1")],
    out_lagged$conditions[[1]]$estimates$average[which(out_lagged$conditions[[1]]$estimates$parameter == "wA3~wA2")]
  )

  capture.output(
    lagged_summary <- summary(out_lagged, parameter = "wA2~wB1")
  )
  expect_equal(lagged_summary$Population, 0.15)

  out_residuals <- powRICLPM(
    target_power = 0.8,
    sample_size = 300,
    time_points = 3,
    ICC = 0.5,
    RI_cor = 0.3,
    lagged_effects = lagged_effects,
    within_cor = 0.3,
    reps = 2,
    seed = 123456,
    constraints = "residuals"
  )

  expect_equal(
    out_residuals$conditions[[1]]$estimates$average[which(out_residuals$conditions[[1]]$estimates$parameter == "wA2~~wA2")],
    out_residuals$conditions[[1]]$estimates$average[which(out_residuals$conditions[[1]]$estimates$parameter == "wA3~~wA3")]
  )

  out_within <- powRICLPM(
    target_power = 0.8,
    sample_size = 300,
    time_points = 3,
    ICC = 0.5,
    RI_cor = 0.3,
    lagged_effects = lagged_effects,
    within_cor = 0.3,
    reps = 2,
    seed = 123456,
    constraints = "within"
  )

  # Cross-lagged effects
  expect_equal(
    out_within$conditions[[1]]$estimates$average[which(out_within$conditions[[1]]$estimates$parameter == "wB2~wA1")],
    out_within$conditions[[1]]$estimates$average[which(out_within$conditions[[1]]$estimates$parameter == "wB3~wA2")]
  )

  # Autoregressive effects
  expect_equal(
    out_within$conditions[[1]]$estimates$average[which(out_within$conditions[[1]]$estimates$parameter == "wA2~wA1")],
    out_within$conditions[[1]]$estimates$average[which(out_within$conditions[[1]]$estimates$parameter == "wA3~wA2")]
  )

  # Residual variances
  expect_equal(
    out_within$conditions[[1]]$estimates$average[which(out_within$conditions[[1]]$estimates$parameter == "wA2~~wA2")],
    out_within$conditions[[1]]$estimates$average[which(out_within$conditions[[1]]$estimates$parameter == "wA3~~wA3")]
  )

  out_stat <- powRICLPM(
    target_power = 0.8,
    sample_size = 300,
    time_points = 3,
    ICC = 0.5,
    RI_cor = 0.3,
    lagged_effects = lagged_effects,
    within_cor = 0.3,
    reps = 2,
    seed = 123456,
    constraints = "stationarity"
  )

  # Cross-lagged effects
  expect_equal(
    out_stat$conditions[[1]]$estimates$average[which(out_stat$conditions[[1]]$estimates$parameter == "wB2~wA1")],
    out_stat$conditions[[1]]$estimates$average[which(out_stat$conditions[[1]]$estimates$parameter == "wB3~wA2")],
    tolerance = 1e06
  )

  # Autoregressive effects
  expect_equal(
    out_stat$conditions[[1]]$estimates$average[which(out_stat$conditions[[1]]$estimates$parameter == "wA2~wA1")],
    out_stat$conditions[[1]]$estimates$average[which(out_stat$conditions[[1]]$estimates$parameter == "wA3~wA2")],
    tolerance = 1e06
  )

  capture.output(
    stationarity_loading_summary <- summary(out_stat, parameter = "wA1=~A1")
  )
  expect_equal(stationarity_loading_summary$Population, 1)

  out_ME <- powRICLPM(
    target_power = 0.8,
    sample_size = 1000,
    time_points = 4,
    ICC = 0.5,
    RI_cor = 0.3,
    lagged_effects = lagged_effects,
    within_cor = 0.3,
    reps = 2,
    seed = 123456,
    constraints = "ME",
    estimate_ME = TRUE
  )

  expect_equal(
    out_ME$conditions[[1]]$estimates$average[which(out_ME$conditions[[1]]$estimates$parameter == "A1~~A1")],
    out_ME$conditions[[1]]$estimates$average[which(out_ME$conditions[[1]]$estimates$parameter == "A2~~A2")]
  )

  #expect_false(
  #  out_ME$conditions[[1]]$estimates$average[which(out_ME$conditions[[1]]$estimates$parameter == "A1~~A1")] ==
  #  out_ME$conditions[[1]]$estimates$average[which(out_ME$conditions[[1]]$estimates$parameter == "B2~~B2")]
  #)
})

