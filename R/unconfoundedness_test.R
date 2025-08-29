# unconfoundedness_test.R
# R >= 4.1

#' @title Test unconfoundedness with comprehensive diagnostics and multiple estimators
#' @description
#' Compare a marginal treatment effect ω in an RCT-like dataset with the same
#' estimand in an observational dataset. Enhanced version supports:
#' - Estimators: IPW, AIPW, TMLE, G-computation, matching
#' - Effect measures: risk difference, risk ratio, odds ratio
#' - Inference: bootstrap CI, analytical SE, robust SE
#' - Transport weighting with auto-detection
#' - Comprehensive diagnostics and visualization
#' - Sensitivity analysis tools
#'
#' @param data_rct,data_obs data.frame with same Y/A definitions and covariates
#' @param formula model formula `Y ~ A + X1 + X2`. A must be binary and the first RHS term
#' @param estimator character: "aipw" (default), "ipw", "tmle", "gcomp", "matching"
#' @param effect_measure character: "rd" (risk difference), "rr" (risk ratio), "or" (odds ratio)
#' @param stabilize logical; stabilized IPW (default TRUE)
#' @param trim length-2 numeric in [0,1] for weight trimming; NULL disables
#' @param family_y "gaussian" or "binomial"
#' @param transport "none", "rct_to_obs", or "auto" (default "none")
#' @param auto_method "both", "ks", "energy" for shift detection
#' @param auto_alpha significance level for shift tests (default 0.01)
#' @param auto_energy_R permutations for energy test (default 199)
#' @param inference_method "bootstrap", "analytical", "robust"
#' @param B bootstrap replicates (default 1000)
#' @param alpha CI level (default 0.05 for 95% CI)
#' @param parallel logical; use parallel processing for bootstrap
#' @param n_cores integer; cores for parallel processing (NULL = auto-detect)
#' @param seed RNG seed or NULL
#' @param validate logical; perform extensive input validation
#' @return enhanced unconf_test object with comprehensive diagnostics
#' @examples
#' \dontrun{
#' # Enhanced example with multiple estimators
#' set.seed(42)
#' data_rct <- generate_rct_data(n = 800)
#' data_obs <- generate_obs_data(n = 2000, confounding = TRUE)
#'
#' # Comprehensive test
#' result <- unconfoundedness_test(
#'   data_rct, data_obs,
#'   Y ~ A + X1 + X2,
#'   estimator = "aipw",
#'   effect_measure = "rd",
#'   transport = "auto",
#'   inference_method = "bootstrap",
#'   B = 1000,
#'   parallel = TRUE,
#'   validate = TRUE
#' )
#'
#' print(result)
#' plot(result)
#' summary(result)
#' }
#' @export
unconfoundedness_test <- function(
  data_rct,
  data_obs,
  formula,
  estimator = c("aipw", "ipw", "tmle", "gcomp", "matching"),
  effect_measure = c("rd", "rr", "or"),
  stabilize = TRUE,
  trim = c(0.01, 0.99),
  family_y = c("gaussian", "binomial"),
  transport = c("none", "rct_to_obs", "auto"),
  auto_method = c("both", "ks", "energy"),
  auto_alpha = 0.01,
  auto_energy_R = 199L,
  inference_method = c("bootstrap", "analytical", "robust"),
  B = 1000L,
  alpha = 0.05,
  parallel = FALSE,
  n_cores = NULL,
  seed = NULL,
  validate = TRUE
) {
  # Match arguments
  estimator <- match.arg(estimator)
  effect_measure <- match.arg(effect_measure)
  family_y <- match.arg(family_y)
  transport <- match.arg(transport)
  auto_method <- match.arg(auto_method)
  inference_method <- match.arg(inference_method)

  # Set seed early if provided
  if (!is.null(seed)) {
    set.seed(as.integer(seed))
  }

  # Enhanced input validation
  if (validate) {
    validation_result <- validate_inputs(
      data_rct,
      data_obs,
      formula,
      estimator,
      effect_measure,
      family_y,
      alpha,
      trim
    )
    if (!validation_result$valid) {
      stop("Input validation failed: ", validation_result$message)
    }
    # Display warnings
    if (length(validation_result$warnings) > 0) {
      for (w in validation_result$warnings) {
        warning(w, call. = FALSE)
      }
    }
  }

  # Parse formula and prepare data
  parsed_data <- parse_formula_and_data(formula, data_rct, data_obs)

  # Auto shift detection and transport decision
  transport_result <- decide_transport(
    parsed_data$X_rct,
    parsed_data$X_obs,
    transport,
    auto_method,
    auto_alpha,
    auto_energy_R
  )

  # Apply transport weights
  weights <- compute_transport_weights(
    parsed_data$X_rct,
    parsed_data$X_obs,
    transport_result$apply_transport
  )

  # Compute effects with chosen estimator
  effects <- compute_effects(
    parsed_data,
    estimator,
    effect_measure,
    family_y,
    weights$w_rct,
    weights$w_obs,
    stabilize,
    trim
  )

  # Compute difference
  diff_estimate <- compute_effect_difference(
    effects$rct_effect,
    effects$obs_effect,
    effect_measure
  )

  # Inference
  inference_result <- compute_inference(
    parsed_data,
    estimator,
    effect_measure,
    family_y,
    weights,
    effects,
    diff_estimate,
    inference_method,
    B,
    alpha,
    parallel,
    n_cores
  )

  # Comprehensive diagnostics
  diagnostics <- compute_diagnostics(
    parsed_data,
    effects,
    weights,
    transport_result,
    estimator,
    family_y
  )

  # Sensitivity analysis
  sensitivity <- compute_sensitivity_analysis(
    effects,
    diagnostics,
    effect_measure
  )

  # Create result object
  result <- create_result_object(
    effects,
    diff_estimate,
    inference_result,
    diagnostics,
    sensitivity,
    transport_result,
    estimator,
    effect_measure,
    family_y,
    inference_method,
    match.call()
  )

  class(result) <- c("unconf_test", "list")
  result
}

#' Input validation function
#' @param data_rct RCT data
#' @param data_obs Observational data
#' @param formula Model formula
#' @param estimator Chosen estimator
#' @param effect_measure Effect measure
#' @param family_y Outcome family
#' @param alpha Alpha level
#' @param trim Trimming bounds
#' @return List with validation results
validate_inputs <- function(
  data_rct,
  data_obs,
  formula,
  estimator,
  effect_measure,
  family_y,
  alpha,
  trim
) {
  warnings <- character(0)

  # Basic data checks
  if (!is.data.frame(data_rct) || !is.data.frame(data_obs)) {
    return(list(
      valid = FALSE,
      message = "data_rct and data_obs must be data.frames"
    ))
  }

  if (nrow(data_rct) < 50) {
    warnings <- c(
      warnings,
      "Small RCT sample size (n < 50) may lead to unstable estimates"
    )
  }

  if (nrow(data_obs) < 100) {
    warnings <- c(
      warnings,
      "Small observational sample size (n < 100) may lead to unstable estimates"
    )
  }

  # Formula validation
  tryCatch(
    {
      terms(formula)
    },
    error = function(e) {
      return(list(valid = FALSE, message = "Invalid formula specification"))
    }
  )

  # Check required variables
  vars <- all.vars(formula)
  if (length(vars) < 2) {
    return(list(
      valid = FALSE,
      message = "Formula must include outcome and treatment"
    ))
  }

  missing_rct <- setdiff(vars, names(data_rct))
  missing_obs <- setdiff(vars, names(data_obs))

  if (length(missing_rct) > 0) {
    return(list(
      valid = FALSE,
      message = paste(
        "Missing variables in RCT data:",
        paste(missing_rct, collapse = ", ")
      )
    ))
  }

  if (length(missing_obs) > 0) {
    return(list(
      valid = FALSE,
      message = paste(
        "Missing variables in observational data:",
        paste(missing_obs, collapse = ", ")
      )
    ))
  }

  # Check for missing values
  rct_missing <- sum(is.na(data_rct[vars]))
  obs_missing <- sum(is.na(data_obs[vars]))

  if (rct_missing > 0) {
    warnings <- c(
      warnings,
      paste("RCT data has", rct_missing, "missing values")
    )
  }

  if (obs_missing > 0) {
    warnings <- c(
      warnings,
      paste("Observational data has", obs_missing, "missing values")
    )
  }

  # Parameter validation
  if (alpha <= 0 || alpha >= 1) {
    return(list(valid = FALSE, message = "alpha must be between 0 and 1"))
  }

  if (!is.null(trim)) {
    if (
      length(trim) != 2 || trim[1] >= trim[2] || any(trim < 0) || any(trim > 1)
    ) {
      return(list(
        valid = FALSE,
        message = "trim must be c(lower, upper) with 0 <= lower < upper <= 1"
      ))
    }
  }

  # Estimator-specific warnings
  if (estimator == "tmle") {
    if (!requireNamespace("tmle", quietly = TRUE)) {
      return(list(
        valid = FALSE,
        message = "Package 'tmle' required for TMLE estimator"
      ))
    }
  }

  if (estimator == "matching") {
    if (!requireNamespace("MatchIt", quietly = TRUE)) {
      return(list(
        valid = FALSE,
        message = "Package 'MatchIt' required for matching estimator"
      ))
    }
  }

  # Effect measure compatibility
  if (effect_measure %in% c("rr", "or") && family_y != "binomial") {
    return(list(
      valid = FALSE,
      message = "Risk ratio and odds ratio only available for binary outcomes"
    ))
  }

  list(valid = TRUE, warnings = warnings)
}

#' Parse formula and prepare data matrices
#' @param formula Model formula
#' @param data_rct RCT data
#' @param data_obs Observational data
#' @return List with parsed components
parse_formula_and_data <- function(formula, data_rct, data_obs) {
  # Parse formula components
  tt <- terms(formula)
  vars <- all.vars(formula)
  y_name <- vars[1]
  rhs_terms <- attr(tt, "term.labels")
  a_name <- strsplit(rhs_terms[1], ":|\\*")[[1]][1]

  # Create model frames
  mf_rct <- model.frame(formula, data_rct)
  mf_obs <- model.frame(formula, data_obs)

  # Extract variables
  Y_rct <- mf_rct[[y_name]]
  A_rct <- mf_rct[[a_name]]
  Y_obs <- mf_obs[[y_name]]
  A_obs <- mf_obs[[a_name]]

  # Coerce treatment to 0/1
  A_rct <- coerce_binary(A_rct, "Treatment in RCT data")
  A_obs <- coerce_binary(A_obs, "Treatment in observational data")

  # Create covariate matrices
  X_terms <- setdiff(rhs_terms, a_name)
  if (length(X_terms) > 0) {
    form_X <- reformulate(X_terms)
    X_rct <- model.matrix(form_X, data = data_rct)
    X_obs <- model.matrix(form_X, data = data_obs)
  } else {
    X_rct <- matrix(
      1,
      nrow = nrow(data_rct),
      dimnames = list(NULL, "(Intercept)")
    )
    X_obs <- matrix(
      1,
      nrow = nrow(data_obs),
      dimnames = list(NULL, "(Intercept)")
    )
  }

  list(
    Y_rct = Y_rct,
    A_rct = A_rct,
    X_rct = X_rct,
    Y_obs = Y_obs,
    A_obs = A_obs,
    X_obs = X_obs,
    y_name = y_name,
    a_name = a_name,
    data_rct = data_rct,
    data_obs = data_obs,
    formula = formula
  )
}

#' Coerce variable to binary 0/1
#' @param x Variable to coerce
#' @param var_name Variable name for error messages
#' @return Binary 0/1 vector
coerce_binary <- function(x, var_name = "Variable") {
  if (is.logical(x)) {
    return(as.integer(x))
  }

  if (is.factor(x)) {
    x <- droplevels(x)
    if (nlevels(x) != 2) {
      stop(paste(var_name, "must have exactly 2 levels"))
    }
    return(as.integer(x) - 1L)
  }

  if (is.character(x)) {
    unique_vals <- unique(x)
    if (length(unique_vals) != 2) {
      stop(paste(var_name, "must have exactly 2 unique values"))
    }
    return(as.integer(x == sort(unique_vals)[2]))
  }

  if (is.numeric(x)) {
    unique_vals <- sort(unique(x))
    if (length(unique_vals) != 2) {
      stop(paste(var_name, "must have exactly 2 unique values"))
    }
    if (all(unique_vals == c(0, 1))) {
      return(as.integer(x))
    } else {
      return(as.integer(x == unique_vals[2]))
    }
  }

  stop(paste(
    var_name,
    "must be logical, factor, character, or numeric with 2 levels"
  ))
}

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

#' Compute effects using specified estimator
#' @param parsed_data Parsed data object
#' @param estimator Estimation method
#' @param effect_measure Effect measure
#' @param family_y Outcome family
#' @param w_rct RCT weights
#' @param w_obs Observational weights
#' @param stabilize Stabilize weights
#' @param trim Trimming bounds
#' @return Effects object
compute_effects <- function(
  parsed_data,
  estimator,
  effect_measure,
  family_y,
  w_rct,
  w_obs,
  stabilize,
  trim
) {
  # Compute RCT effect
  rct_effect <- estimate_single_effect(
    Y = parsed_data$Y_rct,
    A = parsed_data$A_rct,
    X = parsed_data$X_rct,
    weights = w_rct,
    estimator = estimator,
    effect_measure = effect_measure,
    family_y = family_y,
    stabilize = stabilize,
    trim = trim
  )

  # Compute observational effect
  obs_effect <- estimate_single_effect(
    Y = parsed_data$Y_obs,
    A = parsed_data$A_obs,
    X = parsed_data$X_obs,
    weights = w_obs,
    estimator = estimator,
    effect_measure = effect_measure,
    family_y = family_y,
    stabilize = stabilize,
    trim = trim
  )

  list(
    rct_effect = rct_effect,
    obs_effect = obs_effect
  )
}

#' Estimate effect for single dataset
#' @param Y Outcome vector
#' @param A Treatment vector
#' @param X Covariate matrix
#' @param weights Sample weights
#' @param estimator Estimation method
#' @param effect_measure Effect measure
#' @param family_y Outcome family
#' @param stabilize Stabilize weights
#' @param trim Trimming bounds
#' @return Single effect estimate
estimate_single_effect <- function(
  Y,
  A,
  X,
  weights,
  estimator,
  effect_measure,
  family_y,
  stabilize,
  trim
) {
  switch(
    estimator,
    "ipw" = estimate_ipw(
      Y,
      A,
      X,
      weights,
      effect_measure,
      family_y,
      stabilize,
      trim
    ),
    "aipw" = estimate_aipw(
      Y,
      A,
      X,
      weights,
      effect_measure,
      family_y,
      stabilize,
      trim
    ),
    "tmle" = estimate_tmle(Y, A, X, weights, effect_measure, family_y),
    "gcomp" = estimate_gcomp(Y, A, X, weights, effect_measure, family_y),
    "matching" = estimate_matching(Y, A, X, weights, effect_measure, family_y),
    stop("Unknown estimator: ", estimator)
  )
}

# estimators.R - Complete estimator implementations

#' IPW Estimator
#' @param Y Outcome vector
#' @param A Treatment vector
#' @param X Covariate matrix
#' @param weights Sample weights
#' @param effect_measure Effect measure
#' @param family_y Outcome family
#' @param stabilize Stabilize weights
#' @param trim Trimming bounds
#' @return IPW estimate
estimate_ipw <- function(
  Y,
  A,
  X,
  weights,
  effect_measure,
  family_y,
  stabilize,
  trim
) {
  # Fit propensity score model
  ps_fit <- fit_propensity_model(A, X, weights)
  e <- ps_fit$fitted_values

  # Compute IPW weights
  ipw_weights <- compute_ipw_weights(A, e, weights, stabilize)

  # Apply trimming if specified
  if (!is.null(trim)) {
    keep_idx <- apply_weight_trimming(ipw_weights, trim)
    Y <- Y[keep_idx]
    A <- A[keep_idx]
    ipw_weights <- ipw_weights[keep_idx]
    e <- e[keep_idx]
  }

  # Estimate potential outcomes
  mu1 <- weighted.mean(Y[A == 1], ipw_weights[A == 1])
  mu0 <- weighted.mean(Y[A == 0], ipw_weights[A == 0])

  # Compute effect measure
  effect <- compute_effect_measure(mu1, mu0, effect_measure)

  # Compute influence function for variance
  IF <- compute_ipw_influence_function(
    Y,
    A,
    e,
    ipw_weights,
    mu1,
    mu0,
    effect_measure
  )

  list(
    estimate = effect,
    mu1 = mu1,
    mu0 = mu0,
    variance = var(IF, na.rm = TRUE),
    influence_function = IF,
    propensity_scores = e,
    weights = ipw_weights,
    details = list(
      ps_model = ps_fit,
      method = "ipw"
    )
  )
}

#' AIPW (Doubly Robust) Estimator
#' @param Y Outcome vector
#' @param A Treatment vector
#' @param X Covariate matrix
#' @param weights Sample weights
#' @param effect_measure Effect measure
#' @param family_y Outcome family
#' @param stabilize Stabilize weights
#' @param trim Trimming bounds
#' @return AIPW estimate
estimate_aipw <- function(
  Y,
  A,
  X,
  weights,
  effect_measure,
  family_y,
  stabilize,
  trim
) {
  # Fit propensity score model
  ps_fit <- fit_propensity_model(A, X, weights)
  e <- ps_fit$fitted_values

  # Fit outcome regression models
  outcome_fits <- fit_outcome_models(Y, A, X, weights, family_y)
  mu1_hat <- outcome_fits$mu1_pred
  mu0_hat <- outcome_fits$mu0_pred

  # Compute AIPW pseudo-outcomes
  pseudo_outcomes <- compute_aipw_pseudo_outcomes(
    Y,
    A,
    e,
    mu1_hat,
    mu0_hat,
    weights,
    stabilize
  )

  # Apply trimming if specified
  if (!is.null(trim)) {
    base_weights <- weights * (A / e + (1 - A) / (1 - e))
    keep_idx <- apply_weight_trimming(base_weights, trim)
    Y <- Y[keep_idx]
    A <- A[keep_idx]
    e <- e[keep_idx]
    mu1_hat <- mu1_hat[keep_idx]
    mu0_hat <- mu0_hat[keep_idx]
    weights <- weights[keep_idx]
    pseudo_outcomes <- compute_aipw_pseudo_outcomes(
      Y,
      A,
      e,
      mu1_hat,
      mu0_hat,
      weights,
      stabilize
    )
  }

  # Estimate potential outcomes
  mu1 <- weighted.mean(pseudo_outcomes$mu1, weights)
  mu0 <- weighted.mean(pseudo_outcomes$mu0, weights)

  # Compute effect measure
  effect <- compute_effect_measure(mu1, mu0, effect_measure)

  # Compute influence function
  IF <- compute_aipw_influence_function(
    pseudo_outcomes,
    weights,
    mu1,
    mu0,
    effect_measure
  )

  list(
    estimate = effect,
    mu1 = mu1,
    mu0 = mu0,
    variance = var(IF, na.rm = TRUE),
    influence_function = IF,
    propensity_scores = e,
    outcome_predictions = list(mu1 = mu1_hat, mu0 = mu0_hat),
    details = list(
      ps_model = ps_fit,
      outcome_models = outcome_fits,
      method = "aipw"
    )
  )
}

#' TMLE Estimator
#' @param Y Outcome vector
#' @param A Treatment vector
#' @param X Covariate matrix
#' @param weights Sample weights
#' @param effect_measure Effect measure
#' @param family_y Outcome family
#' @return TMLE estimate
estimate_tmle <- function(Y, A, X, weights, effect_measure, family_y) {
  if (!requireNamespace("tmle", quietly = TRUE)) {
    stop("Package 'tmle' required for TMLE estimator")
  }

  # Convert matrix to data frame for tmle
  if (ncol(X) > 1 || !all(X[, 1] == 1)) {
    W <- X[, !colnames(X) %in% "(Intercept)", drop = FALSE]
    if (ncol(W) == 0) W <- NULL
  } else {
    W <- NULL
  }

  # Fit TMLE
  tmle_fit <- tryCatch(
    {
      tmle::tmle(
        Y = Y,
        A = A,
        W = W,
        family = if (family_y == "binomial") "binomial" else "gaussian",
        V = 5 # 5-fold cross-validation
      )
    },
    error = function(e) {
      # Fallback to simpler TMLE if cross-validation fails
      tmle::tmle(
        Y = Y,
        A = A,
        W = W,
        family = if (family_y == "binomial") "binomial" else "gaussian"
      )
    }
  )

  # Extract estimates
  mu1 <- tmle_fit$estimates$EY1$psi
  mu0 <- tmle_fit$estimates$EY0$psi
  effect_rd <- tmle_fit$estimates$ATE$psi

  # Convert to requested effect measure
  if (effect_measure == "rd") {
    effect <- effect_rd
  } else if (effect_measure == "rr") {
    effect <- if (mu0 != 0) mu1 / mu0 else NA_real_
  } else if (effect_measure == "or") {
    if (mu0 != 0 && mu0 != 1 && mu1 != 0 && mu1 != 1) {
      effect <- (mu1 / (1 - mu1)) / (mu0 / (1 - mu0))
    } else {
      effect <- NA_real_
    }
  }

  # Variance (for risk difference)
  variance <- tmle_fit$estimates$ATE$var.psi

  # Adjust variance for other measures using delta method
  if (effect_measure == "rr" && mu0 != 0 && !is.na(effect)) {
    # Delta method for log risk ratio
    grad <- c(1 / mu0, -mu1 / (mu0^2))
    var_mu <- matrix(
      c(tmle_fit$estimates$EY1$var.psi, 0, 0, tmle_fit$estimates$EY0$var.psi),
      2,
      2
    )
    variance <- as.numeric(t(grad) %*% var_mu %*% grad * (effect^2))
  } else if (effect_measure == "or" && !is.na(effect)) {
    # Delta method for log odds ratio (approximate)
    if (mu1 != 0 && mu1 != 1 && mu0 != 0 && mu0 != 1) {
      variance <- tmle_fit$estimates$ATE$var.psi *
        (effect^2) /
        ((mu1 * (1 - mu1) + mu0 * (1 - mu0))^2)
    }
  }

  list(
    estimate = effect,
    mu1 = mu1,
    mu0 = mu0,
    variance = variance,
    influence_function = tmle_fit$estimates$ATE$IC,
    details = list(
      tmle_fit = tmle_fit,
      method = "tmle"
    )
  )
}

#' G-computation Estimator
#' @param Y Outcome vector
#' @param A Treatment vector
#' @param X Covariate matrix
#' @param weights Sample weights
#' @param effect_measure Effect measure
#' @param family_y Outcome family
#' @return G-computation estimate
estimate_gcomp <- function(Y, A, X, weights, effect_measure, family_y) {
  # Fit outcome models
  outcome_fits <- fit_outcome_models(Y, A, X, weights, family_y)

  # Predict potential outcomes for all subjects
  mu1_hat <- outcome_fits$mu1_pred
  mu0_hat <- outcome_fits$mu0_pred

  # G-computation estimates
  mu1 <- weighted.mean(mu1_hat, weights)
  mu0 <- weighted.mean(mu0_hat, weights)

  # Compute effect measure
  effect <- compute_effect_measure(mu1, mu0, effect_measure)

  # Bootstrap-based variance estimation
  n <- length(Y)
  W_scaled <- weights / sum(weights)

  # Influence function approximation for G-computation
  IF <- W_scaled * (mu1_hat - mu0_hat - effect)
  variance <- var(IF, na.rm = TRUE)

  list(
    estimate = effect,
    mu1 = mu1,
    mu0 = mu0,
    variance = variance,
    influence_function = IF,
    outcome_predictions = list(mu1 = mu1_hat, mu0 = mu0_hat),
    details = list(
      outcome_models = outcome_fits,
      method = "gcomp"
    )
  )
}

#' Matching Estimator
#' @param Y Outcome vector
#' @param A Treatment vector
#' @param X Covariate matrix
#' @param weights Sample weights
#' @param effect_measure Effect measure
#' @param family_y Outcome family
#' @return Matching estimate
estimate_matching <- function(Y, A, X, weights, effect_measure, family_y) {
  if (!requireNamespace("MatchIt", quietly = TRUE)) {
    stop("Package 'MatchIt' required for matching estimator")
  }

  # Prepare data
  W <- X[, !colnames(X) %in% "(Intercept)", drop = FALSE]
  if (ncol(W) == 0) {
    stop("Matching requires covariates")
  }

  data_match <- data.frame(Y = Y, A = A, W, weights = weights)

  # Perform matching
  match_fit <- tryCatch(
    {
      MatchIt::matchit(
        A ~ .,
        data = data_match[, !names(data_match) %in% c("Y", "weights")],
        method = "nearest",
        distance = "glm",
        replace = FALSE,
        ratio = 1
      )
    },
    error = function(e) {
      # Fallback to exact matching if available
      if (ncol(W) <= 3) {
        MatchIt::matchit(
          A ~ .,
          data = data_match[, !names(data_match) %in% c("Y", "weights")],
          method = "exact"
        )
      } else {
        stop("Matching failed: ", e$message)
      }
    }
  )

  # Extract matched data
  matched_data <- MatchIt::match.data(match_fit)

  # Compute effect on matched sample
  Y_matched <- matched_data$Y
  A_matched <- matched_data$A
  match_weights <- matched_data$weights

  # Simple difference in means on matched sample
  mu1 <- weighted.mean(Y_matched[A_matched == 1], match_weights[A_matched == 1])
  mu0 <- weighted.mean(Y_matched[A_matched == 0], match_weights[A_matched == 0])

  # Compute effect measure
  effect <- compute_effect_measure(mu1, mu0, effect_measure)

  # Approximate variance using matched sample
  n1 <- sum(A_matched == 1)
  n0 <- sum(A_matched == 0)

  if (n1 > 1 && n0 > 1) {
    s1_sq <- wtd.var(Y_matched[A_matched == 1], match_weights[A_matched == 1])
    s0_sq <- wtd.var(Y_matched[A_matched == 0], match_weights[A_matched == 0])

    if (effect_measure == "rd") {
      variance <- s1_sq / n1 + s0_sq / n0
    } else {
      # Delta method approximation for other measures
      if (effect_measure == "rr" && mu0 != 0) {
        variance <- (s1_sq / n1) / (mu0^2) + (s0_sq / n0) * (mu1^2) / (mu0^4)
      } else if (
        effect_measure == "or" && mu0 != 0 && mu0 != 1 && mu1 != 0 && mu1 != 1
      ) {
        variance <- (s1_sq / n1) /
          (mu1^2 * (1 - mu1)^2) +
          (s0_sq / n0) / (mu0^2 * (1 - mu0)^2)
      } else {
        variance <- NA_real_
      }
    }
  } else {
    variance <- NA_real_
  }

  list(
    estimate = effect,
    mu1 = mu1,
    mu0 = mu0,
    variance = variance,
    influence_function = NULL, # Not easily available for matching
    details = list(
      match_fit = match_fit,
      matched_data = matched_data,
      n_matched = nrow(matched_data),
      method = "matching"
    )
  )
}

# Helper functions for estimators

#' Extract covariate names from a model formula
#'
#' Given a formula like `Y ~ A + X1 + X2`, return the covariate names
#' to be used in the propensity model: all symbols on the RHS except `treat`.
#' Works with interactions (X1:X2), functions (log(X1)), and backticked names.
#' If a variable isn't a column in `data`, it's dropped.
#'
#' @param formula stats::formula like `Y ~ A + X1 + X2`
#' @param data data.frame used to validate which names actually exist
#' @param treat character(1) name of the treatment variable (default "A")
#' @return character() vector of covariate names present in `data`
#' @examples
#' df <- data.frame(Y=0, A=0, X1=1, X2=2, X3=3)
#' extract_covariates(Y ~ A + X1 + X2, df)          # c("X1","X2")
#' extract_covariates(Y ~ A + log(X1) + X2:X3, df)  # c("X1","X2","X3")
#'
#' @keywords internal
extract_covariates <- function(formula, data, treat = "A") {
  # --- type checks ------------------------------------------------------------
  if (!inherits(formula, "formula")) {
    stop("`formula` must be a stats::formula.", call. = FALSE)
  }
  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame.", call. = FALSE)
  }
  if (!is.character(treat) || length(treat) != 1L || is.na(treat) || treat == "") {
    stop("`treat` must be a non-empty character(1).", call. = FALSE)
  }

  # --- pull symbols from RHS only --------------------------------------------
  # all.vars() returns symbol names appearing in the expression,
  # expanding interactions/functions etc. We avoid evaluating anything.
  rhs <- formula[[3L]]
  rhs_vars <- all.vars(rhs)

  # --- drop treatment and keep only columns that exist ------------------------
  covars <- setdiff(rhs_vars, treat)
  covars <- intersect(covars, colnames(data))

  covars
}

#' Build a propensity-score formula safely
#'
#' Constructs a formula of the form \code{A ~ x1 + x2 + ...} from a set of
#' predictor names. If \code{predictors} is empty, it falls back to the
#' intercept-only model \code{A ~ 1} and emits a warning.
#'
#' @param predictors character()
#'   Vector of predictor (covariate) column names. May be length 0.
#'
#' @return formula
#'   A valid formula object for use in \code{glm()}.
#'
#' @examples
#' build_ps_formula(c("X1","X2"))  # A ~ X1 + X2
#' build_ps_formula(character(0))  # A ~ 1 (with warning)
#'
#' @keywords internal
build_ps_formula <- function(predictors) {
  # ---- Type checks -----------------------------------------------------------
  if (!is.character(predictors)) {
    stop("`predictors` must be a character vector of column names.", call. = FALSE)
  }
  if (anyNA(predictors)) {
    stop("`predictors` contains NA values.", call. = FALSE)
  }

  # ---- Empty set => intercept-only ------------------------------------------
  if (length(predictors) == 0L) {
    warning(
      "No covariates available for the propensity model; using intercept-only model (A ~ 1).",
      call. = FALSE
    )
    return(A ~ 1)
  }

  # ---- Construct A ~ x1 + x2 + ... ------------------------------------------
  stats::as.formula(paste("A ~", paste(predictors, collapse = " + ")))
}

#' Fit propensity score model
#'
#' Builds a logistic regression for \code{A} given covariates \code{X}.
#' Handles the cases where \code{X} has only an intercept or no columns
#' by falling back to \code{A ~ 1}. Also stabilizes fitted probabilities to
#' avoid exact 0/1 weights downstream.
#'
#' @param A numeric|integer|logical
#'   Treatment indicator (0/1). Coerced to integer 0/1 if logical.
#' @param X matrix|data.frame
#'   Covariate design (may include an \code{"(Intercept)"} column). Row count
#'   must match \code{length(A)}. Column names are required (added if missing).
#' @param weights numeric|NULL
#'   Optional sampling weights of length \code{length(A)}. Defaults to 1.
#'
#' @return list
#'   \itemize{
#'     \item \code{model}: the \code{glm} fit (or a placeholder list on fallback)
#'     \item \code{fitted_values}: numeric vector of stabilized propensities
#'     \item \code{formula}: the formula used
#'   }
#'
#' @examples
#' set.seed(1)
#' A <- rbinom(100, 1, 0.4)
#' X <- cbind(`(Intercept)` = 1, X1 = rnorm(100), X2 = rnorm(100))
#' fit <- fit_propensity_model(A, X, weights = NULL)
#' head(fit$fitted_values)
#'
#' @keywords internal
fit_propensity_model <- function(A, X, weights = NULL) {
  # ---- Type checks -----------------------------------------------------------
  n <- length(A)
  if (is.logical(A)) A <- as.integer(A)
  if (!is.numeric(A) || any(!A %in% c(0, 1))) {
    stop("`A` must be a binary vector (0/1 or logical).", call. = FALSE)
  }

  if (!is.matrix(X) && !is.data.frame(X)) {
    stop("`X` must be a matrix or data.frame.", call. = FALSE)
  }
  if (nrow(X) != n) {
    stop("`nrow(X)` must equal length(A).", call. = FALSE)
  }

  # Normalize colnames (required for formula terms)
  if (is.null(colnames(X))) {
    colnames(X) <- paste0("X", seq_len(ncol(X)))
  }

  # Weights default to 1
  if (is.null(weights)) {
    weights <- rep(1, n)
  } else {
    if (!is.numeric(weights) || length(weights) != n) {
      stop("`weights` must be numeric and the same length as `A`.", call. = FALSE)
    }
  }

  # ---- Build predictors data.frame (drop intercept-like columns) ------------
  # Intercept columns are commonly named "(Intercept)" (from model.matrix),
  # or sometimes a constant-1 unnamed column. We drop any column that is all 1s
  # or is named "(Intercept)".
  X_df <- as.data.frame(X, check.names = TRUE)
  is_intercept_col <- function(col, nm) {
    (nm == "(Intercept)") || (is.numeric(col) && all(col == 1))
  }
  keep <- !mapply(is_intercept_col, X_df, names(X_df))
  predictors_df <- X_df[, keep, drop = FALSE]

  # ---- Build safe formula ----------------------------------------------------
  predictor_names <- colnames(predictors_df)
  f_ps <- build_ps_formula(predictor_names)

  # ---- Assemble modeling data.frame -----------------------------------------
  df <- if (length(predictor_names)) {
    data.frame(A = A, predictors_df, check.names = TRUE)
  } else {
    data.frame(A = A, check.names = TRUE)
  }

  # ---- Fit model (robust fallback) ------------------------------------------
  fit <- tryCatch(
    stats::glm(f_ps, data = df, family = stats::binomial(), weights = weights),
    error = function(e) {
      # Separation / singularities / other failures: fallback to marginal p
      warning(
        "Propensity score model fitting failed; using marginal treated proportion as PS.",
        call. = FALSE
      )
      list(fitted.values = rep(mean(A), n))
    }
  )

  # ---- Stabilize fitted probabilities ---------------------------------------
  if (inherits(fit, "glm")) {
    p <- stats::fitted(fit)
  } else {
    p <- fit$fitted.values
  }
  # clamp to (1e-3, 1-1e-3) to avoid extreme weights
  p <- pmin(pmax(p, 1e-3), 1 - 1e-3)

  list(
    model = fit,
    fitted_values = p,
    formula = f_ps
  )
}

#' Fit outcome regression models
#' @param Y Outcome vector
#' @param A Treatment vector
#' @param X Covariate matrix
#' @param weights Sample weights
#' @param family_y Outcome family
#' @return Outcome model fits
fit_outcome_models <- function(Y, A, X, weights, family_y) {
  # Prepare predictors
  if (ncol(X) > 1 || !all(X[, 1] == 1)) {
    predictors <- X[, !colnames(X) %in% "(Intercept)", drop = FALSE]
    if (ncol(predictors) == 0) {
      use_predictors <- FALSE
    } else {
      use_predictors <- TRUE
    }
  } else {
    use_predictors <- FALSE
  }

  # Fit model for treated (A=1)
  idx_1 <- A == 1
  if (sum(idx_1) < 2) {
    stop("Insufficient treated observations for outcome modeling")
  }

  if (use_predictors) {
    df_1 <- data.frame(Y = Y[idx_1], predictors[idx_1, , drop = FALSE])
    formula_1 <- as.formula(paste(
      "Y ~",
      paste(colnames(predictors), collapse = " + ")
    ))
  } else {
    df_1 <- data.frame(Y = Y[idx_1])
    formula_1 <- Y ~ 1
  }

  family_obj <- if (family_y == "binomial") binomial() else gaussian()

  fit_1 <- tryCatch(
    {
      glm(formula_1, data = df_1, family = family_obj, weights = weights[idx_1])
    },
    error = function(e) {
      # Fallback to intercept-only model
      glm(
        Y ~ 1,
        data = data.frame(Y = Y[idx_1]),
        family = family_obj,
        weights = weights[idx_1]
      )
    }
  )

  # Fit model for control (A=0)
  idx_0 <- A == 0
  if (sum(idx_0) < 2) {
    stop("Insufficient control observations for outcome modeling")
  }

  if (use_predictors) {
    df_0 <- data.frame(Y = Y[idx_0], predictors[idx_0, , drop = FALSE])
    formula_0 <- as.formula(paste(
      "Y ~",
      paste(colnames(predictors), collapse = " + ")
    ))
  } else {
    df_0 <- data.frame(Y = Y[idx_0])
    formula_0 <- Y ~ 1
  }

  fit_0 <- tryCatch(
    {
      glm(formula_0, data = df_0, family = family_obj, weights = weights[idx_0])
    },
    error = function(e) {
      # Fallback to intercept-only model
      glm(
        Y ~ 1,
        data = data.frame(Y = Y[idx_0]),
        family = family_obj,
        weights = weights[idx_0]
      )
    }
  )

  # Predict for all observations
  if (use_predictors && inherits(fit_1, "glm") && inherits(fit_0, "glm")) {
    newdata <- data.frame(predictors)
    names(newdata) <- colnames(predictors)
  } else {
    newdata <- data.frame(row.names = seq_along(Y))
  }

  mu1_pred <- tryCatch(
    {
      predict(fit_1, newdata = newdata, type = "response")
    },
    error = function(e) {
      rep(mean(Y[idx_1]), length(Y))
    }
  )

  mu0_pred <- tryCatch(
    {
      predict(fit_0, newdata = newdata, type = "response")
    },
    error = function(e) {
      rep(mean(Y[idx_0]), length(Y))
    }
  )

  list(
    fit_1 = fit_1,
    fit_0 = fit_0,
    mu1_pred = as.numeric(mu1_pred),
    mu0_pred = as.numeric(mu0_pred),
    formulas = list(treated = formula_1, control = formula_0)
  )
}

#' Compute IPW weights
#' @param A Treatment vector
#' @param e Propensity scores
#' @param weights Sample weights
#' @param stabilize Whether to stabilize
#' @return IPW weights
compute_ipw_weights <- function(A, e, weights, stabilize) {
  if (stabilize) {
    # Stabilized weights
    p_a <- mean(A) # Marginal treatment probability
    numerator <- A * p_a + (1 - A) * (1 - p_a)
  } else {
    numerator <- 1
  }

  ipw_weights <- weights * numerator / (A * e + (1 - A) * (1 - e))

  # Check for extreme weights
  if (any(is.infinite(ipw_weights)) || any(ipw_weights < 0)) {
    warning("Extreme IPW weights detected, check propensity score model")
    ipw_weights[is.infinite(ipw_weights)] <- max(
      ipw_weights[is.finite(ipw_weights)],
      na.rm = TRUE
    )
    ipw_weights[ipw_weights < 0] <- 0
  }

  ipw_weights
}

#' Apply weight trimming
#' @param weights Weight vector
#' @param trim Trimming bounds
#' @return Indices to keep
apply_weight_trimming <- function(weights, trim) {
  if (is.null(trim)) {
    return(seq_along(weights))
  }

  lower_bound <- quantile(weights, trim[1], na.rm = TRUE)
  upper_bound <- quantile(weights, trim[2], na.rm = TRUE)

  which(weights >= lower_bound & weights <= upper_bound & !is.na(weights))
}

#' Compute effect measure from potential outcome means
#' @param mu1 Mean under treatment
#' @param mu0 Mean under control
#' @param effect_measure Effect measure
#' @return Effect estimate
compute_effect_measure <- function(mu1, mu0, effect_measure) {
  switch(
    effect_measure,
    "rd" = mu1 - mu0, # Risk difference
    "rr" = {
      # Risk ratio
      if (mu0 == 0) {
        warning("Control mean is 0, risk ratio undefined")
        NA_real_
      } else {
        mu1 / mu0
      }
    },
    "or" = {
      # Odds ratio
      if (mu0 == 0 || mu0 == 1 || mu1 == 0 || mu1 == 1) {
        warning("Probabilities at boundary, odds ratio undefined")
        NA_real_
      } else {
        (mu1 / (1 - mu1)) / (mu0 / (1 - mu0))
      }
    },
    stop("Unknown effect measure: ", effect_measure)
  )
}

#' Compute AIPW pseudo-outcomes
#' @param Y Outcome vector
#' @param A Treatment vector
#' @param e Propensity scores
#' @param mu1_hat Outcome predictions under treatment
#' @param mu0_hat Outcome predictions under control
#' @param weights Sample weights
#' @param stabilize Whether to stabilize
#' @return Pseudo-outcomes
compute_aipw_pseudo_outcomes <- function(
  Y,
  A,
  e,
  mu1_hat,
  mu0_hat,
  weights,
  stabilize
) {
  if (stabilize) {
    p_a <- mean(A)
    numerator_1 <- p_a
    numerator_0 <- 1 - p_a
  } else {
    numerator_1 <- numerator_0 <- 1
  }

  # AIPW pseudo-outcomes
  mu1_pseudo <- mu1_hat + numerator_1 * A * (Y - mu1_hat) / e
  mu0_pseudo <- mu0_hat + numerator_0 * (1 - A) * (Y - mu0_hat) / (1 - e)

  list(mu1 = mu1_pseudo, mu0 = mu0_pseudo)
}

#' Compute IPW influence function
#' @param Y Outcome vector
#' @param A Treatment vector
#' @param e Propensity scores
#' @param weights IPW weights
#' @param mu1 Mean under treatment
#' @param mu0 Mean under control
#' @param effect_measure Effect measure
#' @return Influence function
compute_ipw_influence_function <- function(
  Y,
  A,
  e,
  weights,
  mu1,
  mu0,
  effect_measure
) {
  n <- length(Y)

  # IPW components
  w1 <- A / e
  w0 <- (1 - A) / (1 - e)

  # Normalize by sum of weights
  sum_w1 <- sum(w1)
  sum_w0 <- sum(w0)

  if1 <- w1 * (Y - mu1) / sum_w1
  if0 <- w0 * (Y - mu0) / sum_w0

  # Influence function for effect measure
  switch(
    effect_measure,
    "rd" = if1 - if0,
    "rr" = {
      if (mu0 != 0) {
        if1 / mu0 - mu1 * if0 / (mu0^2)
      } else {
        rep(NA_real_, n)
      }
    },
    "or" = {
      # Delta method for log odds ratio
      if (mu1 != 0 && mu1 != 1 && mu0 != 0 && mu0 != 1) {
        d_mu1 <- 1 / (mu1 * (1 - mu1))
        d_mu0 <- -1 / (mu0 * (1 - mu0))
        or_val <- (mu1 / (1 - mu1)) / (mu0 / (1 - mu0))
        or_val * (d_mu1 * if1 + d_mu0 * if0)
      } else {
        rep(NA_real_, n)
      }
    },
    stop("Unknown effect measure: ", effect_measure)
  )
}

#' Compute AIPW influence function
#' @param pseudo_outcomes Pseudo-outcomes
#' @param weights Sample weights
#' @param mu1 Mean under treatment
#' @param mu0 Mean under control
#' @param effect_measure Effect measure
#' @return Influence function
compute_aipw_influence_function <- function(
  pseudo_outcomes,
  weights,
  mu1,
  mu0,
  effect_measure
) {
  # Weighted mean influence functions
  W <- weights / sum(weights)

  if1 <- W * (pseudo_outcomes$mu1 - mu1)
  if0 <- W * (pseudo_outcomes$mu0 - mu0)

  # Influence function for effect measure
  switch(
    effect_measure,
    "rd" = if1 - if0,
    "rr" = {
      if (mu0 != 0) {
        if1 / mu0 - mu1 * if0 / (mu0^2)
      } else {
        rep(NA_real_, length(if1))
      }
    },
    "or" = {
      # Delta method for log odds ratio
      if (mu1 != 0 && mu1 != 1 && mu0 != 0 && mu0 != 1) {
        d_mu1 <- 1 / (mu1 * (1 - mu1))
        d_mu0 <- -1 / (mu0 * (1 - mu0))
        or_val <- (mu1 / (1 - mu1)) / (mu0 / (1 - mu0))
        or_val * (d_mu1 * if1 + d_mu0 * if0)
      } else {
        rep(NA_real_, length(if1))
      }
    },
    stop("Unknown effect measure: ", effect_measure)
  )
}

#' Weighted variance
#' @param x Vector
#' @param weights Weights
#' @return Weighted variance
wtd.var <- function(x, weights = NULL) {
  if (is.null(weights)) {
    return(var(x, na.rm = TRUE))
  }

  # Remove missing values
  valid_idx <- !is.na(x) & !is.na(weights)
  x <- x[valid_idx]
  weights <- weights[valid_idx]

  if (length(x) < 2) {
    return(NA_real_)
  }

  if (all(weights == weights[1])) {
    return(var(x))
  }

  w_mean <- weighted.mean(x, weights)
  sum(weights * (x - w_mean)^2) / (sum(weights) - 1)
}
