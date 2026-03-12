# transport.R - Transport weighting, covariate shift detection

#' Decide whether to apply transport weighting
#' @param X_rct RCT covariate matrix
#' @param X_obs Observational covariate matrix
#' @param transport Transport method
#' @param auto_method Method for auto detection
#' @param auto_alpha Alpha for shift tests
#' @param auto_energy_R Permutations for energy test
#' @return List with transport decision and diagnostics
decide_transport <- function(
  X_rct,
  X_obs,
  transport,
  auto_method,
  auto_alpha,
  auto_energy_R
) {
  if (transport == "none") {
    return(list(
      apply_transport = FALSE,
      method = "none",
      auto_diagnostics = NULL
    ))
  }

  if (transport == "rct_to_obs") {
    return(list(
      apply_transport = TRUE,
      method = "rct_to_obs",
      auto_diagnostics = NULL
    ))
  }

  # Auto detection
  auto_diag <- detect_covariate_shift(
    X_rct,
    X_obs,
    auto_method,
    auto_alpha,
    auto_energy_R
  )

  list(
    apply_transport = auto_diag$shift_detected,
    method = if (auto_diag$shift_detected) "rct_to_obs" else "none",
    auto_diagnostics = auto_diag
  )
}

#' Detect covariate shift between datasets
#' @param X_rct RCT covariate matrix
#' @param X_obs Observational covariate matrix
#' @param method Detection method
#' @param alpha Significance level
#' @param energy_R Permutations for energy test
#' @return List with shift detection results
detect_covariate_shift <- function(X_rct, X_obs, method, alpha, energy_R) {
  results <- list()
  detected <- FALSE

  # Remove intercept for tests
  X_rct_test <- X_rct[, !colnames(X_rct) %in% "(Intercept)", drop = FALSE]
  X_obs_test <- X_obs[, !colnames(X_obs) %in% "(Intercept)", drop = FALSE]

  # KS test
  if (method %in% c("both", "ks") && ncol(X_rct_test) > 0) {
    ks_result <- perform_ks_tests(X_rct_test, X_obs_test, alpha)
    results$ks <- ks_result
    if (ks_result$any_significant) detected <- TRUE
  }

  # Energy test
  if (method %in% c("both", "energy")) {
    energy_result <- perform_energy_test(
      X_rct_test,
      X_obs_test,
      alpha,
      energy_R
    )
    results$energy <- energy_result
    if (energy_result$available && energy_result$significant) detected <- TRUE
  }

  # Combined standardized mean differences
  smd_result <- compute_standardized_mean_differences(X_rct_test, X_obs_test)
  results$smd <- smd_result

  list(
    shift_detected = detected,
    method = method,
    alpha = alpha,
    tests = results
  )
}

#' Perform KS tests for each covariate
#' @param X_rct RCT covariates
#' @param X_obs Observational covariates
#' @param alpha Significance level
#' @return KS test results
perform_ks_tests <- function(X_rct, X_obs, alpha) {
  if (ncol(X_rct) == 0) {
    return(list(
      available = TRUE,
      p_values = numeric(0),
      p_adj = numeric(0),
      any_significant = FALSE,
      significant_vars = character(0)
    ))
  }

  var_names <- colnames(X_rct)
  p_values <- numeric(length(var_names))
  names(p_values) <- var_names

  for (i in seq_along(var_names)) {
    var_name <- var_names[i]
    x_rct <- X_rct[, var_name]
    x_obs <- X_obs[, var_name]

    # Handle discrete variables
    if (length(unique(c(x_rct, x_obs))) <= 10) {
      # Chi-square test for discrete variables
      tab_rct <- table(x_rct)
      tab_obs <- table(x_obs)
      all_levels <- sort(unique(c(names(tab_rct), names(tab_obs))))

      counts_rct <- rep(0, length(all_levels))
      counts_obs <- rep(0, length(all_levels))
      names(counts_rct) <- names(counts_obs) <- all_levels

      counts_rct[names(tab_rct)] <- tab_rct
      counts_obs[names(tab_obs)] <- tab_obs

      test_result <- chisq.test(rbind(counts_rct, counts_obs))
      p_values[i] <- test_result$p.value
    } else {
      # KS test for continuous variables
      test_result <- ks.test(x_rct, x_obs)
      p_values[i] <- test_result$p.value
    }
  }

  # Multiple testing correction
  p_adj <- p.adjust(p_values, method = "BH")
  significant <- p_adj < alpha

  list(
    available = TRUE,
    p_values = p_values,
    p_adj = p_adj,
    any_significant = any(significant),
    significant_vars = names(p_adj)[significant]
  )
}

#' Perform energy test for multivariate shift
#' @param X_rct RCT covariates
#' @param X_obs Observational covariates
#' @param alpha Significance level
#' @param R Number of permutations
#' @return Energy test results
perform_energy_test <- function(X_rct, X_obs, alpha, R) {
  if (!requireNamespace("energy", quietly = TRUE)) {
    return(list(
      available = FALSE,
      p_value = NA_real_,
      significant = FALSE,
      R = R
    ))
  }

  if (ncol(X_rct) == 0) {
    return(list(
      available = TRUE,
      p_value = 1.0,
      significant = FALSE,
      R = R
    ))
  }

  # Combine data
  Z <- rbind(X_rct, X_obs)
  sizes <- c(nrow(X_rct), nrow(X_obs))

  # Perform test
  test_result <- energy::eqdist.etest(Z, sizes = sizes, R = as.integer(R))

  list(
    available = TRUE,
    p_value = test_result$p.value,
    significant = test_result$p.value < alpha,
    R = R,
    statistic = test_result$statistic
  )
}

#' Compute standardized mean differences
#' @param X_rct RCT covariates
#' @param X_obs Observational covariates
#' @return SMD results
compute_standardized_mean_differences <- function(X_rct, X_obs) {
  if (ncol(X_rct) == 0) {
    return(list(
      smd = numeric(0),
      max_smd = 0
    ))
  }

  var_names <- colnames(X_rct)
  smd <- numeric(length(var_names))
  names(smd) <- var_names

  for (i in seq_along(var_names)) {
    var_name <- var_names[i]
    x_rct <- X_rct[, var_name]
    x_obs <- X_obs[, var_name]

    mean_diff <- mean(x_obs) - mean(x_rct)
    pooled_sd <- sqrt((var(x_rct) + var(x_obs)) / 2)

    smd[i] <- if (pooled_sd > 0) mean_diff / pooled_sd else 0
  }

  list(
    smd = smd,
    max_smd = max(abs(smd))
  )
}

#' Compute transport weights
#' @param X_rct RCT covariate matrix
#' @param X_obs Observational covariate matrix
#' @param apply_transport Whether to apply transport weighting
#' @return List with weights
compute_transport_weights <- function(X_rct, X_obs, apply_transport) {
  n_rct <- nrow(X_rct)
  n_obs <- nrow(X_obs)

  if (!apply_transport) {
    return(list(
      w_rct = rep(1, n_rct),
      w_obs = rep(1, n_obs),
      transport_applied = FALSE,
      diagnostics = NULL
    ))
  }

  # Fit density ratio model
  Z <- rbind(X_rct, X_obs)
  indicator <- c(rep(0, n_rct), rep(1, n_obs)) # 0 = RCT, 1 = OBS

  # Remove intercept from predictors if present
  if ("(Intercept)" %in% colnames(Z)) {
    Z <- Z[, !colnames(Z) %in% "(Intercept)", drop = FALSE]
  }

  # Fit logistic regression
  if (ncol(Z) == 0) {
    # No covariates - equal weights
    w_rct <- rep(1, n_rct)
  } else {
    df_fit <- data.frame(S = indicator, Z)
    fit <- glm(S ~ ., data = df_fit, family = binomial())

    # Predict probability of being in observational study
    newdata_rct <- data.frame(Z[1:n_rct, , drop = FALSE])
    names(newdata_rct) <- names(df_fit)[-1]

    p_obs_given_rct <- predict(fit, newdata = newdata_rct, type = "response")
    p_obs_given_rct <- pmax(pmin(p_obs_given_rct, 0.999), 0.001) # Stabilize

    # Compute weights: P(S=1|X) / P(S=0|X)
    w_rct <- p_obs_given_rct / (1 - p_obs_given_rct)
    w_rct <- w_rct / mean(w_rct) # Normalize
  }

  w_obs <- rep(1, n_obs)

  # Compute diagnostics
  diagnostics <- list(
    weight_range = range(w_rct),
    weight_quantiles = quantile(
      w_rct,
      c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99)
    ),
    effective_sample_size = sum(w_rct)^2 / sum(w_rct^2),
    coefficient_of_variation = sd(w_rct) / mean(w_rct)
  )

  list(
    w_rct = w_rct,
    w_obs = w_obs,
    transport_applied = TRUE,
    diagnostics = diagnostics
  )
}
