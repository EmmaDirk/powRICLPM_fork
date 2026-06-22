#' @title
#' Create `lavaan` Syntax for Data Generation and Model Estimation
#'
#' @description
#' \code{create_lavaan()} creates a \pkg{lavaan} parameter table and model syntax for data generation and estimation of the specified random intercept cross-lagged panel model.
#'
#' @param condition A list, containing values for factors that define an experimental condition (i.e., sample size, number of time points, ICC, and reliability), as well as population values for parameters that depend on these (e.g., the random intercept variances, and covariance between them is a function of `ICC`).
#'
#' @return The `condition` list, extended with a `lavaan` population model table and syntax, and `lavaan` estimation model table and syntax.
#'
#' @details
#' \subsection{Naming conventions}{Details on the naming conventions can be found in the "Details" section of \code{\link{powRICLPM}}.}
#'
#' @noRd
create_lavaan <- function(condition) {
  if (identical(condition[["model"]], "DPM")) {
    return(create_lavaan_DPM(condition))
  }

  # Generate default variable names
  name_var <- LETTERS[1:2]

  # Create matrix of names for observed variable, within, and between components
  name_obs <- sapply(name_var, paste0, 1:condition[["time_points"]])
  name_within <- sapply(name_var, function(x) {
    paste0("w", x, 1:condition[["time_points"]])
  })
  name_RI <- paste0("RI_", name_var)

  # Create population parameter table
  pop_tab <- rbind(
    lav_RI(condition, name_RI, name_obs),
    pop_RI_var(condition, name_RI),
    pop_RI_cor(condition, name_RI),
    pop_within(condition, name_within, name_obs),
    pop_lagged(condition, name_within),
    pop_within_var1(condition, name_within),
    pop_within_cov1(condition, name_within),
    pop_within_var2(condition, name_within),
    pop_within_cov2(condition, name_within),
    pop_ME(condition, name_obs)
  )
  rownames(pop_tab) <- NULL

  # Create population model syntax
  pop_synt <- paste0( # Paste over parameters
    paste0( # Paste over columns
      pop_tab[, 1],
      pop_tab[, 2],
      pop_tab[, 3],
      pop_tab[, 4],
      pop_tab[, 5]
    ),
    collapse = "\n"
  )

  # Create estimation parameter table
  est_tab <- rbind(
    est_RI(condition, name_RI, name_obs),
    est_RI_var(condition, name_RI),
    est_RI_cor(condition, name_RI),
    est_within(condition, name_within, name_obs),
    est_lagged(condition, name_within),
    est_within_var1(condition, name_within),
    est_within_cov1(condition, name_within),
    est_within_var2(condition, name_within),
    est_within_cov2(condition, name_within),
    if (condition[["estimate_ME"]]) {
      create_estimate_ME(condition, name_obs)
    },
    if (has_constraint(condition[["constraints"]], "stationarity")) {
      rbind(
        lav_within_cor(condition),
        lav_stationarity(condition)
      )
    }
  )
  rownames(est_tab) <- NULL

  # Create estimation model syntax
  est_synt <- paste0( # Paste over parameters
    paste0( # Paste over columns
      est_tab[, 1],
      est_tab[, 2],
      est_tab[, 3],
      est_tab[, 4],
      est_tab[, 5]
    ),
    collapse = "\n"
  )

  # Create condition list with extra element space
  list(
    sample_size = condition[["sample_size"]],
    model = condition[["model"]],
    time_points = condition[["time_points"]],
    ICC = condition[["ICC"]],
    reliability = condition[["reliability"]],
    reliability_matrix = condition[["reliability_matrix"]],
    loadings = condition[["loadings"]],
    RI_var = condition[["RI_var"]],
    RI_cov = condition[["RI_cov"]],
    ME_var = condition[["ME_var"]],
    pop_synt = pop_synt,
    pop_tab = pop_tab,
    est_synt = est_synt,
    est_tab = est_tab,
    estimate_ME = condition[["estimate_ME"]],
    skewness = condition[["skewness"]],
    kurtosis = condition[["kurtosis"]],
    significance_criterion = condition[["significance_criterion"]],
    estimates = NA,
    MCSEs = NA,
    estimation_information = NA,
    reps = NA,
    condition_id = condition[["condition_id"]]
  )
}

lav_RI <- function(condition, name_RI, name_obs) {
  lhs <- rep(name_RI, each = condition[["time_points"]])
  op <- rep("=~", times = 2 * condition[["time_points"]])
  pv <- c(t(condition[["loadings"]]))
  con <- rep("*", times = 2 * condition[["time_points"]])
  rhs <- c(unlist(name_obs))
  free <- FALSE
  return(cbind.data.frame(lhs, op, pv, con, rhs, free,
    stringsAsFactors = FALSE
  ))
}

est_RI <- function(condition, name_RI, name_obs) {
  lhs <- rep(name_RI, each = condition[["time_points"]])
  op <- rep("=~", times = 2 * condition[["time_points"]])
  con <- rep("*", times = 2 * condition[["time_points"]])
  rhs <- c(unlist(name_obs))

  if (has_constraint(condition[["constraints"]], "loadings_free")) {
    loading_prefixes <- c("lx", "ly")
    pv <- unlist(lapply(seq_along(loading_prefixes), function(i) {
      loadings <- condition[["loadings"]][i, ]
      c("1", paste0(loading_prefixes[i], 2:condition[["time_points"]], "*start(", loadings[-1], ")"))
    }))
    free <- rep(c(FALSE, rep(TRUE, condition[["time_points"]] - 1)), times = 2)
  } else {
    pv <- rep("1", times = 2 * condition[["time_points"]])
    free <- FALSE
  }

  return(cbind.data.frame(lhs, op, pv, con, rhs, free,
    stringsAsFactors = FALSE
  ))
}

pop_RI_var <- function(condition, name_RI) {
  lhs <- rhs <- name_RI
  op <- "~~"
  con <- "*"
  pv <- condition[["RI_var"]]
  free <- FALSE
  return(cbind.data.frame(lhs, op, pv, con, rhs, free,
    stringsAsFactors = FALSE
  ))
}

est_RI_var <- function(condition, name_RI) {
  lhs <- rhs <- name_RI
  op <- "~~"
  con <- "*"
  pv <- paste0("start(", condition[["RI_var"]], ")")
  free <- TRUE
  return(cbind.data.frame(lhs, op, pv, con, rhs, free,
    stringsAsFactors = FALSE
  ))
}

pop_RI_cor <- function(condition, name_RI) {
  lhs <- name_RI[1]
  rhs <- name_RI[2]
  op <- "~~"
  con <- "*"
  pv <- condition[["RI_cov"]]
  free <- FALSE
  return(cbind.data.frame(lhs, op, pv, con, rhs, free,
    stringsAsFactors = FALSE
  ))
}

est_RI_cor <- function(condition, name_RI) {
  lhs <- name_RI[[1]]
  rhs <- name_RI[[2]]
  op <- "~~"
  con <- "*"
  pv <- paste0("start(", condition[["RI_cov"]], ")")
  free <- TRUE
  return(cbind.data.frame(lhs, op, pv, con, rhs, free,
    stringsAsFactors = FALSE
  ))
}

pop_within <- function(condition, name_within, name_obs) {
  lhs <- c(name_within)
  op <- "=~"
  pv <- "1"
  con <- "*"
  rhs <- c(name_obs)
  free <- FALSE
  return(cbind.data.frame(lhs, op, pv, con, rhs, free,
    stringsAsFactors = FALSE
  ))
}

est_within <- function(condition, name_within, name_obs) {
  op <- "=~"
  if (has_constraint(condition[["constraints"]], "stationarity")) {
    lhs <- c(name_within)
    pv <- rep("NA*start(1)", times = length(name_within))
    free <- TRUE
    rhs <- c(name_obs)
  } else {
    lhs <- c(name_within)
    pv <- "1"
    free <- FALSE
    rhs <- c(name_obs)
  }
  con <- "*"
  return(cbind.data.frame(lhs, op, pv, con, rhs, free,
    stringsAsFactors = FALSE
  ))
}

pop_lagged <- function(condition, name_within) {
  lhs <- rep(c(t(name_within))[-(1:2)], each = 2)
  op <- "~"
  con <- "*"
  pv <- c(t(condition[["lagged_effects"]][[1]]))
  free <- FALSE
  rhs <- c(apply(name_within[-condition[["time_points"]], ], 1, rep, times = 2))
  return(cbind.data.frame(lhs, op, pv, con, rhs, free,
    stringsAsFactors = FALSE
  ))
}

est_lagged <- function(condition, name_within) {
  op <- "~"
  con <- "*"
  if (has_constraint(condition[["constraints"]], "lagged")) {
    lhs <- rep(c(t(name_within))[-(1:2)], each = 2)
    pv <- paste0(
      rep(c("a", "b", "c", "d"), times = condition[["time_points"]] - 1),
      "*",
      rep(
        paste0("start(", c(t(condition[["lagged_effects"]][[1]])), ")"),
        times = condition[["time_points"]] - 1
      )
    )
    rhs <- c(apply(name_within[-condition[["time_points"]], ], 1, rep, times = 2))
  } else {
    lhs <- rep(c(t(name_within))[-(1:2)], each = 2)
    pv <- paste0("start(", c(t(condition[["lagged_effects"]][[1]])), ")")
    rhs <- c(apply(name_within[-condition[["time_points"]], ], 1, rep, times = 2))
  }
  free <- TRUE
  return(cbind.data.frame(lhs, op, pv, con, rhs, free,
    stringsAsFactors = FALSE
  ))
}

pop_within_var1 <- function(condition, name_within) {
  lhs <- rhs <- c(t(name_within[1, ]))
  op <- "~~"
  con <- "*"
  pv <- "1"
  free <- FALSE
  return(cbind.data.frame(lhs, op, pv, con, rhs, free,
    stringsAsFactors = FALSE
  ))
}

est_within_var1 <- function(condition, name_within) {
  lhs <- rhs <- c(t(name_within[1, ]))
  op <- "~~"
  con <- "*"
  if (has_constraint(condition[["constraints"]], "stationarity")) {
    pv <- "1"
    free <- FALSE
  } else {
    pv <- paste0("start(1)")
    free <- TRUE
  }
  return(cbind.data.frame(lhs, op, pv, con, rhs, free,
    stringsAsFactors = FALSE
  ))
}

pop_within_cov1 <- function(condition, name_within) {
  lhs <- name_within[1, "A"]
  rhs <- name_within[1, "B"]
  op <- "~~"
  con <- "*"
  pv <- c(condition[["within_cor"]]) # Get covariances
  free <- FALSE
  return(cbind.data.frame(lhs, op, pv, con, rhs, free,
    stringsAsFactors = FALSE
  ))
}

est_within_cov1 <- function(condition, name_within) {
  op <- "~~"
  con <- "*"
  if (has_constraint(condition[["constraints"]], "stationarity")) { # Label
    lhs <- name_within[1, 1]
    rhs <- name_within[1, 2]
    pv <- paste0("cor1*start(", condition[["within_cor"]], ")")
  } else { # Freely estimate
    lhs <- name_within[1, 1]
    rhs <- name_within[1, 2]
    pv <- paste0("start(", c(condition[["within_cor"]]), ")")
  }
  free <- TRUE
  return(cbind.data.frame(lhs, op, pv, con, rhs, free,
    stringsAsFactors = FALSE
  ))
}

pop_within_var2 <- function(condition, name_within) {
  lhs <- rhs <- c(name_within[-1, ])
  op <- "~~"
  con <- "*"
  pv <- rep(diag(condition[["Psi"]][[1]]), each = (condition[["time_points"]] - 1))
  free <- FALSE
  return(cbind.data.frame(lhs, op, pv, con, rhs, free,
    stringsAsFactors = FALSE
  ))
}

est_within_var2 <- function(condition, name_within) {
  op <- "~~"
  con <- "*"
  if (has_constraint(condition[["constraints"]], "residuals")) { # Constrain over time
    lhs <- rhs <- c(name_within[-1, ])
    pv <- paste0(
      rep(c("rvarA", "rvarB"), each = (condition[["time_points"]] - 1)),
      "*",
      paste0("start(", rep(diag(condition[["Psi"]][[1]]), each = (condition[["time_points"]] - 1)), ")")
    )
  } else if (has_constraint(condition[["constraints"]], "stationarity")) {
    lhs <- rhs <- c(name_within[-1, ])
    pv <- paste0(
      c(
        paste0("rvarA", 2:condition[["time_points"]]),
        paste0("rvarB", 2:condition[["time_points"]])
      ),
      "*",
      paste0("start(", rep(diag(condition[["Psi"]][[1]]), each = (condition[["time_points"]] - 1)), ")")
    )
  } else { # Freely estimate
    lhs <- rhs <- c(name_within[-1, ])
    pv <- paste0("start(", rep(diag(condition[["Psi"]][[1]]), each = (condition[["time_points"]] - 1)), ")")
  }
  free <- TRUE
  return(cbind.data.frame(lhs, op, pv, con, rhs, free,
    stringsAsFactors = FALSE
  ))
}

pop_within_cov2 <- function(condition, name_within) {
  lhs <- name_within[-1, 1]
  rhs <- name_within[-1, 2]
  op <- "~~"
  con <- "*"
  pv <- c(condition[["Psi"]][[1]][lower.tri(condition[["Psi"]][[1]])])
  free <- FALSE
  return(cbind.data.frame(lhs, op, pv, con, rhs, free,
    stringsAsFactors = FALSE
  ))
}

est_within_cov2 <- function(condition, name_within) {
  op <- "~~"
  con <- "*"
  if (has_constraint(condition[["constraints"]], "residuals")) { # Constrain over time
    lhs <- name_within[-1, 1]
    rhs <- name_within[-1, 2]
    pv <- paste0(
      rep("rcov", times = condition[["time_points"]] - 1),
      "*",
      rep(paste0("start(", c(condition[["Psi"]][[1]][lower.tri(condition[["Psi"]][[1]])]), ")"), times = condition[["time_points"]] - 1)
    )
  } else if (has_constraint(condition[["constraints"]], "stationarity")) { # Label
    lhs <- c(name_within[-1, 1])
    rhs <- c(name_within[-1, 2])
    pv <- paste0("rcov", 2:condition[["time_points"]])
  } else { # Freely estimate
    lhs <- name_within[-1, 1]
    rhs <- name_within[-1, 2]
    pv <- paste0("start(", c(condition[["Psi"]][[1]][lower.tri(condition[["Psi"]][[1]])]), ")")
  }
  free <- TRUE
  return(cbind.data.frame(lhs, op, pv, con, rhs, free,
    stringsAsFactors = FALSE
  ))
}

pop_ME <- function(condition, name_obs) {
  lhs <- rhs <- c(name_obs)
  op <- "~~"
  pv <- c(t(condition[["ME_var"]]))
  free <- FALSE
  con <- "*"
  return(cbind.data.frame(lhs, op, pv, con, rhs, free,
    stringsAsFactors = FALSE
  ))
}

create_estimate_ME <- function(condition, name_obs) {
  op <- "~~"
  con <- "*"
  ME_var <- condition[["ME_var"]]
  if (has_constraint(condition[["constraints"]], "stationarity") ||
      has_constraint(condition[["constraints"]], "ME")) {
    lhs <- rhs <- c(name_obs)
    constrained_starts <- rep(rowMeans(ME_var), each = condition[["time_points"]])
    pv <- paste0(
      rep(c("MEvarA", "MEvarB"), each = condition[["time_points"]]),
      "*",
      paste0("start(", constrained_starts, ")")
    )
    free <- TRUE
  } else { # Freely estimate ME
    lhs <- rhs <- c(name_obs)
    pv <- paste0("start(", c(t(ME_var)), ")")
    free <- TRUE
  }
  return(cbind.data.frame(lhs, op, pv, con, rhs, free,
    stringsAsFactors = FALSE
  ))
}

estimate_ME0 <- function(condition, name_obs) {
  lhs <- rhs <- c(name_obs)
  op <- "~~"
  pv <- 0
  free <- FALSE
  con <- "*"
  return(cbind.data.frame(lhs, op, pv, con, rhs, free,
    stringsAsFactors = FALSE
  ))
}

lav_within_cor <- function(condition) {
  lhs <- paste0("cor", 2:(condition[["time_points"]] - 1))
  op <- ":="
  pv <- con <- ""
  rhs1 <- paste0("a*c + b*d + a*d*cor", 1:(condition[["time_points"]] - 2))
  rhs2 <- paste0(" + b*c*cor", 1:(condition[["time_points"]] - 2))
  rhs3 <- paste0(" + rcov", 2:(condition[["time_points"]] - 1))
  rhs <- paste0(rhs1, rhs2, rhs3)
  free <- FALSE
  return(cbind.data.frame(lhs, op, pv, con, rhs, free,
    stringsAsFactors = FALSE
  ))
}

lav_stationarity <- function(condition) {
  lhs <- c(
    paste0("rvarA", 2:condition[["time_points"]]),
    paste0("rvarB", 2:condition[["time_points"]])
  )
  op <- "=="
  rhs <- c(
    paste0("1 - (a^2 + b^2 + 2*a*b*cor", 1:(condition[["time_points"]] - 1), ")"),
    paste0("1 - (c^2 + d^2 + 2*c*d*cor", 1:(condition[["time_points"]] - 1), ")")
  )
  pv <- con <- ""
  free <- FALSE
  return(cbind.data.frame(lhs, op, pv, con, rhs, free,
    stringsAsFactors = FALSE
  ))
}

create_lavaan_DPM <- function(condition) {
  name_var <- LETTERS[1:2]
  name_obs <- sapply(name_var, paste0, 1:condition[["time_points"]])
  name_AF <- paste0("AF_", name_var)

  pop_tab <- rbind(
    pop_DPM_AF(condition, name_AF, name_obs),
    pop_DPM_AF_var(condition, name_AF),
    pop_DPM_AF_cov(condition, name_AF),
    pop_DPM_baseline_var(condition, name_obs),
    pop_DPM_baseline_cov(condition, name_obs),
    pop_DPM_AF_baseline_cov(condition, name_AF, name_obs),
    pop_DPM_lagged(condition, name_obs),
    pop_DPM_residual_var(condition, name_obs),
    pop_DPM_residual_cov(condition, name_obs)
  )
  rownames(pop_tab) <- NULL
  pop_synt <- lavaan_table_to_syntax(pop_tab)

  est_tab <- rbind(
    est_DPM_AF(condition, name_AF, name_obs),
    est_DPM_AF_var(condition, name_AF),
    est_DPM_AF_cov(condition, name_AF),
    est_DPM_baseline_var(condition, name_obs),
    est_DPM_baseline_cov(condition, name_obs),
    est_DPM_AF_baseline_cov(condition, name_AF, name_obs),
    est_DPM_lagged(condition, name_obs),
    est_DPM_residual_var(condition, name_obs),
    est_DPM_residual_cov(condition, name_obs),
    if (has_constraint(condition[["constraints"]], "stationarity")) {
      DPM_stationarity_constraints(condition)
    }
  )
  rownames(est_tab) <- NULL
  est_synt <- lavaan_table_to_syntax(est_tab)

  list(
    sample_size = condition[["sample_size"]],
    model = condition[["model"]],
    time_points = condition[["time_points"]],
    ICC = condition[["ICC"]],
    reliability = condition[["reliability"]],
    reliability_matrix = condition[["reliability_matrix"]],
    loadings = condition[["loadings"]],
    AF_var = condition[["AF_var"]],
    AF_cov = condition[["AF_cov"]],
    DPM_values = condition[["DPM_values"]],
    pop_synt = pop_synt,
    pop_tab = pop_tab,
    est_synt = est_synt,
    est_tab = est_tab,
    estimate_ME = condition[["estimate_ME"]],
    skewness = condition[["skewness"]],
    kurtosis = condition[["kurtosis"]],
    significance_criterion = condition[["significance_criterion"]],
    estimates = NA,
    MCSEs = NA,
    estimation_information = NA,
    reps = NA,
    condition_id = condition[["condition_id"]]
  )
}

lavaan_table_to_syntax <- function(tab) {
  paste0(
    paste0(tab[, 1], tab[, 2], tab[, 3], tab[, 4], tab[, 5]),
    collapse = "\n"
  )
}

pop_DPM_AF <- function(condition, name_AF, name_obs) {
  lhs <- rep(name_AF, each = condition[["time_points"]] - 1L)
  op <- rep("=~", times = 2 * (condition[["time_points"]] - 1L))
  pv <- c(t(condition[["loadings"]]))
  con <- rep("*", times = length(lhs))
  rhs <- c(unlist(name_obs[-1, , drop = FALSE]))
  free <- FALSE
  cbind.data.frame(lhs, op, pv, con, rhs, free, stringsAsFactors = FALSE)
}

est_DPM_AF <- function(condition, name_AF, name_obs) {
  lhs <- rep(name_AF, each = condition[["time_points"]] - 1L)
  op <- rep("=~", times = 2 * (condition[["time_points"]] - 1L))
  con <- rep("*", times = length(lhs))
  rhs <- c(unlist(name_obs[-1, , drop = FALSE]))

  if (has_constraint(condition[["constraints"]], "loadings_free")) {
    loading_prefixes <- c("lx", "ly")
    pv <- unlist(lapply(seq_along(loading_prefixes), function(i) {
      loadings <- condition[["loadings"]][i, ]
      if (condition[["time_points"]] == 2L) {
        return("1")
      }
      c(
        "1",
        paste0(
          loading_prefixes[i], 3:condition[["time_points"]],
          "*start(", loadings[-1], ")"
        )
      )
    }))
    free <- rep(c(FALSE, rep(TRUE, condition[["time_points"]] - 2L)), times = 2)
  } else {
    pv <- rep("1", times = length(lhs))
    free <- FALSE
  }

  cbind.data.frame(lhs, op, pv, con, rhs, free, stringsAsFactors = FALSE)
}

pop_DPM_AF_var <- function(condition, name_AF) {
  lhs <- rhs <- name_AF
  op <- "~~"
  pv <- condition[["AF_var"]]
  con <- "*"
  free <- FALSE
  cbind.data.frame(lhs, op, pv, con, rhs, free, stringsAsFactors = FALSE)
}

est_DPM_AF_var <- function(condition, name_AF) {
  lhs <- rhs <- name_AF
  op <- "~~"
  pv <- paste0(c("vfa", "vfb"), "*start(", condition[["AF_var"]], ")")
  con <- "*"
  free <- TRUE
  cbind.data.frame(lhs, op, pv, con, rhs, free, stringsAsFactors = FALSE)
}

pop_DPM_AF_cov <- function(condition, name_AF) {
  lhs <- name_AF[1]
  rhs <- name_AF[2]
  op <- "~~"
  pv <- condition[["AF_cov"]]
  con <- "*"
  free <- FALSE
  cbind.data.frame(lhs, op, pv, con, rhs, free, stringsAsFactors = FALSE)
}

est_DPM_AF_cov <- function(condition, name_AF) {
  lhs <- name_AF[1]
  rhs <- name_AF[2]
  op <- "~~"
  pv <- paste0("cff*start(", condition[["AF_cov"]], ")")
  con <- "*"
  free <- TRUE
  cbind.data.frame(lhs, op, pv, con, rhs, free, stringsAsFactors = FALSE)
}

pop_DPM_baseline_var <- function(condition, name_obs) {
  lhs <- rhs <- c(t(name_obs[1, ]))
  op <- "~~"
  pv <- "1"
  con <- "*"
  free <- FALSE
  cbind.data.frame(lhs, op, pv, con, rhs, free, stringsAsFactors = FALSE)
}

est_DPM_baseline_var <- function(condition, name_obs) {
  lhs <- rhs <- c(t(name_obs[1, ]))
  op <- "~~"
  con <- "*"
  if (has_constraint(condition[["constraints"]], "stationarity")) {
    pv <- "1"
    free <- FALSE
  } else {
    pv <- "start(1)"
    free <- TRUE
  }
  cbind.data.frame(lhs, op, pv, con, rhs, free, stringsAsFactors = FALSE)
}

pop_DPM_baseline_cov <- function(condition, name_obs) {
  lhs <- name_obs[1, 1]
  rhs <- name_obs[1, 2]
  op <- "~~"
  pv <- condition[["within_cor"]]
  con <- "*"
  free <- FALSE
  cbind.data.frame(lhs, op, pv, con, rhs, free, stringsAsFactors = FALSE)
}

est_DPM_baseline_cov <- function(condition, name_obs) {
  lhs <- name_obs[1, 1]
  rhs <- name_obs[1, 2]
  op <- "~~"
  con <- "*"
  if (has_constraint(condition[["constraints"]], "stationarity")) {
    pv <- paste0("wc*start(", condition[["within_cor"]], ")")
  } else {
    pv <- paste0("start(", condition[["within_cor"]], ")")
  }
  free <- TRUE
  cbind.data.frame(lhs, op, pv, con, rhs, free, stringsAsFactors = FALSE)
}

pop_DPM_AF_baseline_cov <- function(condition, name_AF, name_obs) {
  lhs <- rep(name_AF, each = 2)
  rhs <- rep(c(t(name_obs[1, ])), times = 2)
  op <- "~~"
  pv <- c(t(condition[["DPM_values"]]$gamma_start))
  con <- "*"
  free <- FALSE
  cbind.data.frame(lhs, op, pv, con, rhs, free, stringsAsFactors = FALSE)
}

est_DPM_AF_baseline_cov <- function(condition, name_AF, name_obs) {
  lhs <- rep(name_AF, each = 2)
  rhs <- rep(c(t(name_obs[1, ])), times = 2)
  op <- "~~"
  labels <- c("gAA1", "gAB1", "gBA1", "gBB1")
  pv <- paste0(labels, "*start(", c(t(condition[["DPM_values"]]$gamma_start)), ")")
  con <- "*"
  free <- TRUE
  cbind.data.frame(lhs, op, pv, con, rhs, free, stringsAsFactors = FALSE)
}

pop_DPM_lagged <- function(condition, name_obs) {
  lhs <- rep(c(t(name_obs))[-(1:2)], each = 2)
  op <- "~"
  con <- "*"
  pv <- c(t(condition[["lagged_effects"]][[1]]))
  free <- FALSE
  rhs <- c(apply(name_obs[-condition[["time_points"]], , drop = FALSE], 1, rep, times = 2))
  cbind.data.frame(lhs, op, pv, con, rhs, free, stringsAsFactors = FALSE)
}

est_DPM_lagged <- function(condition, name_obs) {
  lhs <- rep(c(t(name_obs))[-(1:2)], each = 2)
  op <- "~"
  con <- "*"
  starts <- paste0("start(", c(t(condition[["lagged_effects"]][[1]])), ")")
  if (has_constraint(condition[["constraints"]], "lagged")) {
    pv <- paste0(
      rep(c("a", "b", "c", "d"), times = condition[["time_points"]] - 1L),
      "*",
      rep(starts, times = condition[["time_points"]] - 1L)
    )
  } else {
    pv <- starts
  }
  free <- TRUE
  rhs <- c(apply(name_obs[-condition[["time_points"]], , drop = FALSE], 1, rep, times = 2))
  cbind.data.frame(lhs, op, pv, con, rhs, free, stringsAsFactors = FALSE)
}

pop_DPM_residual_var <- function(condition, name_obs) {
  lhs <- rhs <- c(name_obs[-1, , drop = FALSE])
  op <- "~~"
  con <- "*"
  pv <- c(
    condition[["DPM_values"]]$Psi[1, 1, ],
    condition[["DPM_values"]]$Psi[2, 2, ]
  )
  free <- FALSE
  cbind.data.frame(lhs, op, pv, con, rhs, free, stringsAsFactors = FALSE)
}

est_DPM_residual_var <- function(condition, name_obs) {
  lhs <- rhs <- c(name_obs[-1, , drop = FALSE])
  op <- "~~"
  con <- "*"
  starts <- c(
    condition[["DPM_values"]]$Psi[1, 1, ],
    condition[["DPM_values"]]$Psi[2, 2, ]
  )
  if (has_constraint(condition[["constraints"]], "stationarity")) {
    labels <- c(
      paste0("rvarA", 2:condition[["time_points"]]),
      paste0("rvarB", 2:condition[["time_points"]])
    )
    pv <- paste0(labels, "*start(", starts, ")")
  } else if (has_constraint(condition[["constraints"]], "residuals")) {
    pv <- paste0(
      rep(c("rvarA", "rvarB"), each = condition[["time_points"]] - 1L),
      "*start(", starts, ")"
    )
  } else {
    pv <- paste0("start(", starts, ")")
  }
  free <- TRUE
  cbind.data.frame(lhs, op, pv, con, rhs, free, stringsAsFactors = FALSE)
}

pop_DPM_residual_cov <- function(condition, name_obs) {
  lhs <- name_obs[-1, 1]
  rhs <- name_obs[-1, 2]
  op <- "~~"
  con <- "*"
  pv <- condition[["DPM_values"]]$Psi[1, 2, ]
  free <- FALSE
  cbind.data.frame(lhs, op, pv, con, rhs, free, stringsAsFactors = FALSE)
}

est_DPM_residual_cov <- function(condition, name_obs) {
  lhs <- name_obs[-1, 1]
  rhs <- name_obs[-1, 2]
  op <- "~~"
  con <- "*"
  starts <- condition[["DPM_values"]]$Psi[1, 2, ]
  if (has_constraint(condition[["constraints"]], "stationarity")) {
    pv <- paste0("rcov", 2:condition[["time_points"]], "*start(", starts, ")")
  } else if (has_constraint(condition[["constraints"]], "residuals")) {
    pv <- paste0("rcov*start(", starts, ")")
  } else {
    pv <- paste0("start(", starts, ")")
  }
  free <- TRUE
  cbind.data.frame(lhs, op, pv, con, rhs, free, stringsAsFactors = FALSE)
}

DPM_stationarity_constraints <- function(condition) {
  if (has_constraint(condition[["constraints"]], "loadings_free")) {
    return(DPM_stationarity_constraints_free_loadings(condition))
  }
  DPM_stationarity_constraints_fixed_loadings(condition)
}

DPM_stationarity_constraints_fixed_loadings <- function(condition) {
  waves <- 2:condition[["time_points"]]
  lhs <- c(
    "gAA1", "gAB1", "gBA1", "gBB1",
    paste0("rvarA", waves),
    paste0("rvarB", waves),
    paste0("rcov", waves)
  )
  op <- "=="
  pv <- con <- ""
  rhs <- c(
    "vfa + gAA1*a + gAB1*b",
    "cff + gAA1*c + gAB1*d",
    "cff + gBA1*a + gBB1*b",
    "vfb + gBA1*c + gBB1*d",
    rep("1 - vfa - (a^2 + b^2 + 2*a*b*wc) - 2*(gAA1*a + gAB1*b)", length(waves)),
    rep("1 - vfb - (c^2 + d^2 + 2*c*d*wc) - 2*(gBA1*c + gBB1*d)", length(waves)),
    rep(
      paste0(
        "wc - cff - (a*c + b*d + (a*d + b*c)*wc) - ",
        "(gAA1*c + gAB1*d) - (gBA1*a + gBB1*b)"
      ),
      length(waves)
    )
  )
  free <- FALSE
  cbind.data.frame(lhs, op, pv, con, rhs, free, stringsAsFactors = FALSE)
}

DPM_stationarity_constraints_free_loadings <- function(condition) {
  waves <- 2:condition[["time_points"]]
  gamma_tab <- DPM_gamma_definitions_free_loadings(condition)
  residual_tab <- DPM_residual_constraints_free_loadings(condition, waves)
  rbind(gamma_tab, residual_tab)
}

DPM_gamma_definitions_free_loadings <- function(condition) {
  waves <- 2:condition[["time_points"]]
  lhs <- character()
  rhs <- character()

  for (wave in waves) {
    prev <- wave - 1L
    lx <- if (wave == 2L) "1" else paste0("lx", wave)
    ly <- if (wave == 2L) "1" else paste0("ly", wave)
    lhs <- c(lhs, paste0(c("gAA", "gAB", "gBA", "gBB"), wave))
    rhs <- c(
      rhs,
      paste0("vfa*", lx, " + gAA", prev, "*a + gAB", prev, "*b"),
      paste0("cff*", ly, " + gAA", prev, "*c + gAB", prev, "*d"),
      paste0("cff*", lx, " + gBA", prev, "*a + gBB", prev, "*b"),
      paste0("vfb*", ly, " + gBA", prev, "*c + gBB", prev, "*d")
    )
  }

  op <- ":="
  pv <- con <- ""
  free <- FALSE
  cbind.data.frame(lhs, op, pv, con, rhs, free, stringsAsFactors = FALSE)
}

DPM_residual_constraints_free_loadings <- function(condition, waves) {
  lhs <- c(
    paste0("rvarA", waves),
    paste0("rvarB", waves),
    paste0("rcov", waves)
  )
  rhs_A <- rhs_B <- rhs_cov <- character()
  for (wave in waves) {
    prev <- wave - 1L
    lx <- if (wave == 2L) "1" else paste0("lx", wave)
    ly <- if (wave == 2L) "1" else paste0("ly", wave)
    rhs_A <- c(
      rhs_A,
      paste0(
        "1 - ", lx, "^2*vfa - (a^2 + b^2 + 2*a*b*wc) - ",
        "2*", lx, "*(gAA", prev, "*a + gAB", prev, "*b)"
      )
    )
    rhs_B <- c(
      rhs_B,
      paste0(
        "1 - ", ly, "^2*vfb - (c^2 + d^2 + 2*c*d*wc) - ",
        "2*", ly, "*(gBA", prev, "*c + gBB", prev, "*d)"
      )
    )
    rhs_cov <- c(
      rhs_cov,
      paste0(
        "wc - ", lx, "*", ly, "*cff - (a*c + b*d + (a*d + b*c)*wc) - ",
        lx, "*(gAA", prev, "*c + gAB", prev, "*d) - ",
        ly, "*(gBA", prev, "*a + gBB", prev, "*b)"
      )
    )
  }

  op <- "=="
  pv <- con <- ""
  rhs <- c(rhs_A, rhs_B, rhs_cov)
  free <- FALSE
  cbind.data.frame(lhs, op, pv, con, rhs, free, stringsAsFactors = FALSE)
}
