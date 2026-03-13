#' Analytical inference using influence functions
#' @param effects Effect estimates
#' @param diff_estimate Difference estimate
#' @param effect_measure Effect measure
#' @param alpha Alpha level
#' @return Analytical inference results
analytical_inference <- function(
  effects,
  diff_estimate,
  effect_measure,
  alpha
) {
  # Check if influence functions are available
  if (
    is.null(effects$rct_effect$influence_function) ||
      is.null(effects$obs_effect$influence_function)
  ) {
    stop("Influence functions not available for analytical inference")
  }

  # Compute variance of difference
  var_rct <- effects$rct_effect$variance
  var_obs <- effects$obs_effect$variance

  # Assume independence between RCT and observational studies
  var_diff <- var_rct + var_obs
  se <- sqrt(var_diff)

  # Test statistic and p-value
  test_stat <- diff_estimate / se
  p_value <- 2 * (1 - pnorm(abs(test_stat)))

  # Confidence interval
  critical_value <- qnorm(1 - alpha / 2)
  ci <- c(
    lower = diff_estimate - critical_value * se,
    upper = diff_estimate + critical_value * se
  )

  list(
    ci = ci,
    se = se,
    p_value = p_value,
    test_stat = test_stat
  )
}

#' Robust inference with sandwich estimator
#' @param parsed_data Parsed data
#' @param effects Effect estimates
#' @param diff_estimate Difference estimate
#' @param effect_measure Effect measure
#' @param alpha Alpha level
#' @return Robust inference results
robust_inference <- function(
  parsed_data,
  effects,
  diff_estimate,
  effect_measure,
  alpha
) {
  # For now, use analytical inference with robust standard errors
  # In practice, would implement proper sandwich estimator
  if (
    !is.null(effects$rct_effect$influence_function) &&
      !is.null(effects$obs_effect$influence_function)
  ) {
    return(analytical_inference(
      effects,
      diff_estimate,
      effect_measure,
      alpha
    ))
  } else {
    # Fallback to bootstrap
    warning("Robust inference not fully implemented, using bootstrap")
    return(bootstrap_inference(
      parsed_data,
      "aipw",
      effect_measure,
      "gaussian",
      FALSE,
      diff_estimate,
      500,
      alpha,
      FALSE,
      NULL
    ))
  }
}

#' Compute comprehensive diagnostics
#' @param parsed_data Parsed data
#' @param effects Effect estimates
#' @param weights Weight object
#' @param transport_result Transport results
#' @param estimator Estimator used
#' @param family_y Outcome family
#' @return Comprehensive diagnostics
compute_diagnostics <- function(
  parsed_data,
  effects,
  weights,
  transport_result,
  estimator,
  family_y
) {
  # Overlap/positivity diagnostics
  overlap_diag <- compute_overlap_diagnostics(
    parsed_data$A_rct,
    parsed_data$A_obs,
    effects$rct_effect$propensity_scores,
    effects$obs_effect$propensity_scores
  )

  # Balance diagnostics
  balance_diag <- compute_balance_diagnostics(
    parsed_data$X_rct,
    parsed_data$X_obs,
    parsed_data$A_rct,
    parsed_data$A_obs,
    weights$w_rct,
    weights$w_obs
  )

  # Model diagnostics
  model_diag <- compute_model_diagnostics(effects, family_y)

  # Sample size diagnostics
  sample_diag <- compute_sample_diagnostics(
    parsed_data,
    weights$w_rct,
    weights$w_obs
  )

  # Transport diagnostics (if applicable)
  transport_diag <- if (weights$transport_applied) {
    weights$diagnostics
  } else {
    NULL
  }

  # Auto-detection diagnostics
  auto_diag <- transport_result$auto_diagnostics

  list(
    overlap = overlap_diag,
    balance = balance_diag,
    model = model_diag,
    sample_size = sample_diag,
    transport = transport_diag,
    auto_detection = auto_diag
  )
}

#' Compute overlap/positivity diagnostics
#' @param A_rct RCT treatment
#' @param A_obs Observational treatment
#' @param e_rct RCT propensity scores
#' @param e_obs Observational propensity scores
#' @return Overlap diagnostics
compute_overlap_diagnostics <- function(A_rct, A_obs, e_rct, e_obs) {
  # Handle estimators that don't produce propensity scores (e.g., gcomp)
  if (is.null(e_rct) || is.null(e_obs) || length(e_rct) == 0 || length(e_obs) == 0) {
    prev_rct <- mean(A_rct)
    prev_obs <- mean(A_obs)
    return(list(
      propensity_score_range = list(rct = c(NA, NA), obs = c(NA, NA)),
      extreme_ps_proportion = list(rct = 0, obs = 0),
      treatment_prevalence = list(rct = prev_rct, obs = prev_obs),
      overlap_statistic = NA_real_,
      overlap_quality = "Unknown (no propensity scores)"
    ))
  }

  # Propensity score ranges
  ps_range_rct <- range(e_rct)
  ps_range_obs <- range(e_obs)

  # Extreme propensity scores
  extreme_threshold <- c(0.01, 0.99)
  extreme_rct <- mean(
    e_rct < extreme_threshold[1] | e_rct > extreme_threshold[2]
  )
  extreme_obs <- mean(
    e_obs < extreme_threshold[1] | e_obs > extreme_threshold[2]
  )

  # Overlap between RCT and observational PS distributions
  overlap_stat <- compute_overlap_statistic(e_rct, e_obs)

  # Treatment prevalence
  prev_rct <- mean(A_rct)
  prev_obs <- mean(A_obs)

  list(
    propensity_score_range = list(rct = ps_range_rct, obs = ps_range_obs),
    extreme_ps_proportion = list(rct = extreme_rct, obs = extreme_obs),
    treatment_prevalence = list(rct = prev_rct, obs = prev_obs),
    overlap_statistic = overlap_stat,
    overlap_quality = classify_overlap_quality(
      extreme_rct,
      extreme_obs,
      overlap_stat
    )
  )
}

#' Compute overlap statistic between two distributions
#' @param x1 First distribution
#' @param x2 Second distribution
#' @return Overlap statistic
compute_overlap_statistic <- function(x1, x2) {
  # Compute overlap using kernel density estimation
  range_x <- range(c(x1, x2))
  x_grid <- seq(range_x[1], range_x[2], length.out = 100)

  # Kernel density estimates
  d1 <- density(x1, from = range_x[1], to = range_x[2], n = 100)
  d2 <- density(x2, from = range_x[1], to = range_x[2], n = 100)

  # Overlap coefficient
  overlap <- sum(pmin(d1$y, d2$y)) / sum(pmax(d1$y, d2$y))

  overlap
}

#' Classify overlap quality
#' @param extreme_rct Proportion of extreme PS in RCT
#' @param extreme_obs Proportion of extreme PS in observational
#' @param overlap_stat Overlap statistic
#' @return Overlap quality classification
classify_overlap_quality <- function(extreme_rct, extreme_obs, overlap_stat) {
  if (extreme_rct > 0.1 || extreme_obs > 0.1 || overlap_stat < 0.1) {
    "Poor"
  } else if (extreme_rct > 0.05 || extreme_obs > 0.05 || overlap_stat < 0.3) {
    "Fair"
  } else if (overlap_stat > 0.7) {
    "Excellent"
  } else {
    "Good"
  }
}

#' Compute balance diagnostics
#' @param X_rct RCT covariates
#' @param X_obs Observational covariates
#' @param A_rct RCT treatment
#' @param A_obs Observational treatment
#' @param w_rct RCT weights
#' @param w_obs Observational weights
#' @return Balance diagnostics
compute_balance_diagnostics <- function(
  X_rct,
  X_obs,
  A_rct,
  A_obs,
  w_rct,
  w_obs
) {
  # Remove intercept column
  X_rct_vars <- X_rct[, !colnames(X_rct) %in% "(Intercept)", drop = FALSE]
  X_obs_vars <- X_obs[, !colnames(X_obs) %in% "(Intercept)", drop = FALSE]

  if (ncol(X_rct_vars) == 0) {
    return(list(
      covariate_balance = NULL,
      treatment_balance = list(
        rct = compute_treatment_balance(X_rct, A_rct, w_rct),
        obs = compute_treatment_balance(X_obs, A_obs, w_obs)
      )
    ))
  }

  # Standardized mean differences between datasets
  smd_between <- compute_standardized_mean_differences(X_rct_vars, X_obs_vars)

  # Balance within each dataset (treated vs control)
  balance_rct <- compute_treatment_balance(X_rct_vars, A_rct, w_rct)
  balance_obs <- compute_treatment_balance(X_obs_vars, A_obs, w_obs)

  # Overall balance assessment
  balance_quality <- assess_balance_quality(
    smd_between$smd,
    balance_rct$smd,
    balance_obs$smd
  )

  list(
    between_datasets = smd_between,
    within_rct = balance_rct,
    within_obs = balance_obs,
    overall_quality = balance_quality
  )
}

#' Compute treatment balance within dataset
#' @param X Covariate matrix
#' @param A Treatment vector
#' @param weights Sample weights
#' @return Treatment balance
compute_treatment_balance <- function(X, A, weights) {
  if (ncol(X) == 0 || all(colnames(X) == "(Intercept)")) {
    return(list(smd = numeric(0), max_smd = 0))
  }

  X_vars <- X[, !colnames(X) %in% "(Intercept)", drop = FALSE]
  var_names <- colnames(X_vars)
  smd <- numeric(length(var_names))
  names(smd) <- var_names

  for (i in seq_along(var_names)) {
    var_name <- var_names[i]
    x <- X_vars[, var_name]

    # Weighted means by treatment group
    mean_1 <- weighted.mean(x[A == 1], weights[A == 1])
    mean_0 <- weighted.mean(x[A == 0], weights[A == 0])

    # Weighted variances
    var_1 <- wtd.var(x[A == 1], weights[A == 1])
    var_0 <- wtd.var(x[A == 0], weights[A == 0])

    # Pooled standard deviation
    pooled_sd <- sqrt((var_1 + var_0) / 2)

    # Standardized mean difference
    smd[i] <- if (pooled_sd > 0) (mean_1 - mean_0) / pooled_sd else 0
  }

  list(smd = smd, max_smd = max(abs(smd)))
}

#' Assess overall balance quality
#' @param smd_between SMD between datasets
#' @param smd_rct SMD within RCT
#' @param smd_obs SMD within observational
#' @return Balance quality assessment
assess_balance_quality <- function(smd_between, smd_rct, smd_obs) {
  max_between <- if (length(smd_between) > 0) max(abs(smd_between)) else 0
  max_within <- max(c(
    if (length(smd_rct) > 0) max(abs(smd_rct)) else 0,
    if (length(smd_obs) > 0) max(abs(smd_obs)) else 0
  ))

  if (max_between > 0.25 || max_within > 0.25) {
    "Poor"
  } else if (max_between > 0.1 || max_within > 0.1) {
    "Fair"
  } else {
    "Good"
  }
}

#' Compute model diagnostics
#' @param effects Effect estimates
#' @param family_y Outcome family
#' @return Model diagnostics
compute_model_diagnostics <- function(effects, family_y) {
  diagnostics <- list()

  # Propensity score model diagnostics
  if (!is.null(effects$rct_effect$details$ps_model)) {
    diagnostics$rct_ps_model <- extract_model_diagnostics(
      effects$rct_effect$details$ps_model$model
    )
  }

  if (!is.null(effects$obs_effect$details$ps_model)) {
    diagnostics$obs_ps_model <- extract_model_diagnostics(
      effects$obs_effect$details$ps_model$model
    )
  }

  # Outcome model diagnostics (for AIPW/G-comp)
  if (!is.null(effects$rct_effect$details$outcome_models)) {
    diagnostics$rct_outcome_models <- list(
      treated = extract_model_diagnostics(
        effects$rct_effect$details$outcome_models$fit_1
      ),
      control = extract_model_diagnostics(
        effects$rct_effect$details$outcome_models$fit_0
      )
    )
  }

  if (!is.null(effects$obs_effect$details$outcome_models)) {
    diagnostics$obs_outcome_models <- list(
      treated = extract_model_diagnostics(
        effects$obs_effect$details$outcome_models$fit_1
      ),
      control = extract_model_diagnostics(
        effects$obs_effect$details$outcome_models$fit_0
      )
    )
  }

  diagnostics
}

#' Extract diagnostics from GLM model
#' @param model GLM model object
#' @return Model diagnostics
extract_model_diagnostics <- function(model) {
  if (is.null(model) || !inherits(model, "glm")) {
    return(NULL)
  }

  list(
    aic = AIC(model),
    bic = BIC(model),
    deviance = deviance(model),
    null_deviance = model$null.deviance,
    df_residual = df.residual(model),
    converged = model$converged,
    pseudo_r2 = 1 - deviance(model) / model$null.deviance
  )
}

#' Compute sample size diagnostics
#' @param parsed_data Parsed data
#' @param w_rct RCT weights
#' @param w_obs Observational weights
#' @return Sample size diagnostics
compute_sample_diagnostics <- function(parsed_data, w_rct, w_obs) {
  # Raw sample sizes
  n_rct <- length(parsed_data$Y_rct)
  n_obs <- length(parsed_data$Y_obs)

  # Treatment group sizes
  n_rct_treated <- sum(parsed_data$A_rct == 1)
  n_rct_control <- sum(parsed_data$A_rct == 0)
  n_obs_treated <- sum(parsed_data$A_obs == 1)
  n_obs_control <- sum(parsed_data$A_obs == 0)

  # Effective sample sizes (accounting for weights)
  ess_rct <- sum(w_rct)^2 / sum(w_rct^2)
  ess_obs <- sum(w_obs)^2 / sum(w_obs^2)

  # Weight efficiency
  weight_efficiency_rct <- ess_rct / n_rct
  weight_efficiency_obs <- ess_obs / n_obs

  list(
    raw_n = list(rct = n_rct, obs = n_obs),
    treatment_n = list(
      rct_treated = n_rct_treated,
      rct_control = n_rct_control,
      obs_treated = n_obs_treated,
      obs_control = n_obs_control
    ),
    effective_n = list(rct = ess_rct, obs = ess_obs),
    weight_efficiency = list(
      rct = weight_efficiency_rct,
      obs = weight_efficiency_obs
    ),
    minimum_group_size = min(
      n_rct_treated,
      n_rct_control,
      n_obs_treated,
      n_obs_control
    )
  )
}

#' Compute sensitivity analysis
#' @param effects Effect estimates
#' @param diagnostics Diagnostic results
#' @param effect_measure Effect measure
#' @return Sensitivity analysis results
compute_sensitivity_analysis <- function(effects, diagnostics, effect_measure) {
  # E-value calculation for unmeasured confounding
  e_value <- compute_e_value(
    effects$obs_effect$estimate,
    effects$rct_effect$estimate,
    effect_measure
  )

  # Fragility assessment
  fragility <- assess_fragility(diagnostics)

  # Robustness checks
  robustness <- assess_robustness(effects, diagnostics)

  list(
    e_value = e_value,
    fragility = fragility,
    robustness = robustness
  )
}

#' Compute E-value for sensitivity to unmeasured confounding
#' @param obs_effect Observational effect
#' @param rct_effect RCT effect
#' @param effect_measure Effect measure
#' @return E-value
compute_e_value <- function(obs_effect, rct_effect, effect_measure) {
  if (effect_measure == "rd") {
    # For risk difference, convert to risk ratio approximation
    # This is approximate and assumes baseline risk
    baseline_risk <- 0.1 # Assumed baseline risk
    rr_obs <- (baseline_risk + obs_effect) / baseline_risk
    rr_rct <- (baseline_risk + rct_effect) / baseline_risk
    ratio_of_rr <- rr_obs / rr_rct
  } else if (effect_measure == "rr") {
    ratio_of_rr <- obs_effect / rct_effect
  } else if (effect_measure == "or") {
    # Convert OR to RR approximation (requires baseline risk assumption)
    baseline_risk <- 0.1
    rr_obs <- obs_effect / (1 - baseline_risk + baseline_risk * obs_effect)
    rr_rct <- rct_effect / (1 - baseline_risk + baseline_risk * rct_effect)
    ratio_of_rr <- rr_obs / rr_rct
  }

  # Guard against NA/NaN/Inf
  if (is.na(ratio_of_rr) || !is.finite(ratio_of_rr)) {
    return(list(
      e_value = NA_real_,
      interpretation = interpret_e_value(NA_real_)
    ))
  }

  # E-value formula
  if (abs(ratio_of_rr - 1) < 1e-6) {
    e_value <- 1.0
  } else {
    e_value <- ratio_of_rr + sqrt(ratio_of_rr * (ratio_of_rr - 1))
    if (ratio_of_rr < 1) {
      e_value <- 1 / e_value
    }
  }

  list(
    e_value = e_value,
    interpretation = interpret_e_value(e_value)
  )
}

#' Interpret E-value
#' @param e_value E-value
#' @return Interpretation
interpret_e_value <- function(e_value) {
  if (is.na(e_value) || !is.finite(e_value)) {
    return("Cannot compute (undefined)")
  }
  if (e_value < 1.5) {
    "Very fragile to unmeasured confounding"
  } else if (e_value < 2.0) {
    "Fragile to unmeasured confounding"
  } else if (e_value < 3.0) {
    "Moderately robust to unmeasured confounding"
  } else {
    "Robust to unmeasured confounding"
  }
}

#' Assess fragility of results
#' @param diagnostics Diagnostic results
#' @return Fragility assessment
assess_fragility <- function(diagnostics) {
  fragility_factors <- character(0)

  # Check overlap quality
  if (diagnostics$overlap$overlap_quality %in% c("Poor", "Fair")) {
    fragility_factors <- c(fragility_factors, "Poor overlap")
  }

  # Check balance quality
  if (
    !is.null(diagnostics$balance$overall_quality) &&
      diagnostics$balance$overall_quality == "Poor"
  ) {
    fragility_factors <- c(fragility_factors, "Poor covariate balance")
  }

  # Check sample sizes
  if (diagnostics$sample_size$minimum_group_size < 50) {
    fragility_factors <- c(fragility_factors, "Small sample size")
  }

  # Check weight efficiency
  if (
    !is.null(diagnostics$transport) &&
      diagnostics$transport$effective_sample_size <
        0.5 * diagnostics$sample_size$raw_n$rct
  ) {
    fragility_factors <- c(
      fragility_factors,
      "Inefficient transport weights"
    )
  }

  list(
    factors = fragility_factors,
    overall = if (length(fragility_factors) == 0) {
      "Robust"
    } else if (length(fragility_factors) <= 2) {
      "Moderate"
    } else {
      "High"
    }
  )
}

#' Assess robustness of results
#' @param effects Effect estimates
#' @param diagnostics Diagnostic results
#' @return Robustness assessment
assess_robustness <- function(effects, diagnostics) {
  # Model specification robustness
  model_robust <- assess_model_robustness(diagnostics$model)

  # Estimation method robustness
  method_robust <- assess_method_robustness(effects)

  # Data quality robustness
  data_robust <- assess_data_robustness(diagnostics)

  list(
    model_specification = model_robust,
    estimation_method = method_robust,
    data_quality = data_robust
  )
}

#' Assess model robustness
#' @param model_diagnostics Model diagnostics
#' @return Model robustness assessment
assess_model_robustness <- function(model_diagnostics) {
  if (is.null(model_diagnostics)) {
    return("Cannot assess - no model diagnostics")
  }

  issues <- character(0)

  # Check for convergence issues
  for (model_type in names(model_diagnostics)) {
    model_info <- model_diagnostics[[model_type]]
    if (is.list(model_info)) {
      # Nested models (e.g., outcome models)
      for (sub_model in names(model_info)) {
        if (
          is.list(model_info[[sub_model]]) &&
            !is.null(model_info[[sub_model]]$converged) &&
            !model_info[[sub_model]]$converged
        ) {
          issues <- c(
            issues,
            paste("Convergence failure in", model_type, sub_model)
          )
        }
      }
    } else if (!is.null(model_info$converged) && !model_info$converged) {
      issues <- c(issues, paste("Convergence failure in", model_type))
    }
  }

  if (length(issues) == 0) {
    "Good model fit"
  } else {
    paste("Issues:", paste(issues, collapse = "; "))
  }
}

#' Assess estimation method robustness
#' @param effects Effect estimates
#' @return Method robustness assessment
assess_method_robustness <- function(effects) {
  rct_method <- effects$rct_effect$details$method
  obs_method <- effects$obs_effect$details$method

  if (rct_method == "aipw" && obs_method == "aipw") {
    "High (doubly robust estimator)"
  } else if (
    rct_method %in% c("ipw", "gcomp") && obs_method %in% c("ipw", "gcomp")
  ) {
    "Moderate (single robust estimator)"
  } else {
    "Variable (mixed methods)"
  }
}

#' Assess data quality robustness
#' @param diagnostics Diagnostic results
#' @return Data robustness assessment
assess_data_robustness <- function(diagnostics) {
  issues <- character(0)

  if (diagnostics$overlap$overlap_quality == "Poor") {
    issues <- c(issues, "Poor overlap")
  }

  if (
    !is.null(diagnostics$balance$overall_quality) &&
      diagnostics$balance$overall_quality == "Poor"
  ) {
    issues <- c(issues, "Poor balance")
  }

  if (diagnostics$sample_size$minimum_group_size < 30) {
    issues <- c(issues, "Very small sample")
  }

  if (length(issues) == 0) {
    "Good data quality"
  } else {
    paste("Issues:", paste(issues, collapse = "; "))
  }
}

#' Create comprehensive result object
#' @param effects Effect estimates
#' @param diff_estimate Difference estimate
#' @param inference_result Inference results
#' @param diagnostics Diagnostic results
#' @param sensitivity Sensitivity analysis
#' @param transport_result Transport results
#' @param estimator Estimator used
#' @param effect_measure Effect measure
#' @param family_y Outcome family
#' @param inference_method Inference method
#' @param call Function call
#' @return Complete result object
create_result_object <- function(
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
  call
) {
  # Create summary table
  summary_df <- data.frame(
    Dataset = c("RCT", "Observational", "Difference"),
    Estimate = c(
      effects$rct_effect$estimate,
      effects$obs_effect$estimate,
      diff_estimate
    ),
    SE = c(
      sqrt(effects$rct_effect$variance %||% NA_real_),
      sqrt(effects$obs_effect$variance %||% NA_real_),
      inference_result$standard_error
    ),
    CI_Lower = c(
      NA_real_,
      NA_real_,
      inference_result$confidence_interval[1]
    ),
    CI_Upper = c(
      NA_real_,
      NA_real_,
      inference_result$confidence_interval[2]
    ),
    P_Value = c(NA_real_, NA_real_, inference_result$p_value),
    stringsAsFactors = FALSE
  )

  # Main result object
  result <- list(
    # Main results
    estimates = list(
      rct = effects$rct_effect$estimate,
      observational = effects$obs_effect$estimate,
      difference = diff_estimate
    ),

    # Inference
    inference = inference_result,

    # Diagnostics
    diagnostics = diagnostics,

    # Sensitivity analysis
    sensitivity = sensitivity,

    # Method information
    methods = list(
      estimator = estimator,
      effect_measure = effect_measure,
      outcome_family = family_y,
      inference_method = inference_method,
      transport_method = transport_result$method,
      transport_applied = transport_result$apply_transport
    ),

    # Raw effects (for advanced users)
    raw_effects = effects,

    # Summary
    summary = summary_df,

    # Call
    call = call
  )

  result
}

#' Compute effect difference
#' @param rct_effect RCT effect estimate
#' @param obs_effect Observational effect estimate
#' @param effect_measure Effect measure
#' @return Difference estimate
compute_effect_difference <- function(rct_effect, obs_effect, effect_measure) {
  switch(
    effect_measure,
    "rd" = obs_effect$estimate - rct_effect$estimate,
    "rr" = obs_effect$estimate / rct_effect$estimate, # Ratio of ratios
    "or" = obs_effect$estimate / rct_effect$estimate, # Ratio of odds ratios
    stop("Unknown effect measure: ", effect_measure)
  )
}

#' Compute inference (confidence intervals and p-values)
#' @param parsed_data Parsed data
#' @param estimator Estimator
#' @param effect_measure Effect measure
#' @param family_y Outcome family
#' @param weights Weight object
#' @param effects Effect estimates
#' @param diff_estimate Difference estimate
#' @param inference_method Inference method
#' @param B Bootstrap replicates
#' @param alpha Alpha level
#' @param parallel Use parallel processing
#' @param n_cores Number of cores
#' @return Inference results
compute_inference <- function(
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
) {
  if (inference_method == "bootstrap") {
    # Bootstrap inference
    boot_results <- bootstrap_inference(
      parsed_data,
      estimator,
      effect_measure,
      family_y,
      weights$transport_applied,
      diff_estimate,
      B,
      alpha,
      parallel,
      n_cores
    )

    return(list(
      method = "bootstrap",
      confidence_interval = boot_results$ci,
      standard_error = boot_results$se,
      p_value = boot_results$p_value,
      bootstrap_distribution = boot_results$boot_diffs,
      test_statistic = boot_results$test_stat
    ))
  } else if (inference_method == "analytical") {
    # Analytical inference using influence functions
    analytical_results <- analytical_inference(
      effects,
      diff_estimate,
      effect_measure,
      alpha
    )

    return(list(
      method = "analytical",
      confidence_interval = analytical_results$ci,
      standard_error = analytical_results$se,
      p_value = analytical_results$p_value,
      test_statistic = analytical_results$test_stat
    ))
  } else if (inference_method == "robust") {
    # Robust inference with sandwich estimator
    robust_results <- robust_inference(
      parsed_data,
      effects,
      diff_estimate,
      effect_measure,
      alpha
    )

    return(list(
      method = "robust",
      confidence_interval = robust_results$ci,
      standard_error = robust_results$se,
      p_value = robust_results$p_value,
      test_statistic = robust_results$test_stat
    ))
  }
}

#' Bootstrap inference
#' @param parsed_data Parsed data
#' @param estimator Estimator
#' @param effect_measure Effect measure
#' @param family_y Outcome family
#' @param transport_applied Whether transport was applied
#' @param diff_estimate Original difference estimate
#' @param B Bootstrap replicates
#' @param alpha Alpha level
#' @param parallel Use parallel processing
#' @param n_cores Number of cores
#' @return Bootstrap results
bootstrap_inference <- function(
  parsed_data,
  estimator,
  effect_measure,
  family_y,
  transport_applied,
  diff_estimate,
  B,
  alpha,
  parallel,
  n_cores
) {
  if (parallel && requireNamespace("parallel", quietly = TRUE)) {
    if (is.null(n_cores)) {
      n_cores <- max(1, parallel::detectCores() - 1)
    }

    # Set up cluster
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl))

    # Export package functions and data to workers
    pkg_ns <- asNamespace("unconfoundedr")
    parallel::clusterExport(
      cl,
      c("parsed_data", "estimator", "effect_measure",
        "family_y", "transport_applied"),
      envir = environment()
    )
    parallel::clusterExport(
      cl,
      c("single_bootstrap_replicate", "estimate_single_effect",
        "compute_effect_difference", "compute_transport_weights",
        "estimate_ipw", "estimate_aipw", "estimate_gcomp",
        "estimate_matching", "estimate_tmle",
        "fit_propensity_model", "fit_outcome_models",
        "build_ps_formula", "compute_ipw_weights",
        "apply_weight_trimming", "compute_effect_measure",
        "compute_aipw_pseudo_outcomes",
        "compute_ipw_influence_function",
        "compute_aipw_influence_function",
        "wtd.var"),
      envir = pkg_ns
    )

    # Run bootstrap in parallel
    boot_diffs <- parallel::parLapply(cl, 1:B, function(b) {
      single_bootstrap_replicate(
        parsed_data,
        estimator,
        effect_measure,
        family_y,
        transport_applied
      )
    })
    boot_diffs <- unlist(boot_diffs)
  } else {
    # Sequential bootstrap
    boot_diffs <- numeric(B)
    for (b in 1:B) {
      boot_diffs[b] <- single_bootstrap_replicate(
        parsed_data,
        estimator,
        effect_measure,
        family_y,
        transport_applied
      )
    }
  }

  # Remove failed replicates
  boot_diffs <- boot_diffs[!is.na(boot_diffs)]

  # Compute statistics
  se <- sd(boot_diffs)
  ci <- quantile(boot_diffs, c(alpha / 2, 1 - alpha / 2))
  test_stat <- diff_estimate / se
  p_value <- 2 * (1 - pnorm(abs(test_stat)))

  list(
    ci = ci,
    se = se,
    p_value = p_value,
    test_stat = test_stat,
    boot_diffs = boot_diffs
  )
}

#' Single bootstrap replicate
#' @param parsed_data Parsed data
#' @param estimator Estimator
#' @param effect_measure Effect measure
#' @param family_y Outcome family
#' @param transport_applied Whether transport was applied
#' @return Bootstrap difference estimate
single_bootstrap_replicate <- function(
  parsed_data,
  estimator,
  effect_measure,
  family_y,
  transport_applied
) {
  # Bootstrap sample indices
  n_rct <- nrow(parsed_data$data_rct)
  n_obs <- nrow(parsed_data$data_obs)

  idx_rct <- sample(1:n_rct, n_rct, replace = TRUE)
  idx_obs <- sample(1:n_obs, n_obs, replace = TRUE)

  # Bootstrap samples
  Y_rct_b <- parsed_data$Y_rct[idx_rct]
  A_rct_b <- parsed_data$A_rct[idx_rct]
  X_rct_b <- parsed_data$X_rct[idx_rct, , drop = FALSE]

  Y_obs_b <- parsed_data$Y_obs[idx_obs]
  A_obs_b <- parsed_data$A_obs[idx_obs]
  X_obs_b <- parsed_data$X_obs[idx_obs, , drop = FALSE]

  # Recompute transport weights for bootstrap sample
  if (transport_applied) {
    weights_b <- compute_transport_weights(X_rct_b, X_obs_b, TRUE)
  } else {
    weights_b <- list(
      w_rct = rep(1, length(idx_rct)),
      w_obs = rep(1, length(idx_obs))
    )
  }

  # Estimate effects on bootstrap samples (skip degenerate resamples)
  rct_effect_b <- tryCatch(
    estimate_single_effect(
      Y_rct_b,
      A_rct_b,
      X_rct_b,
      weights_b$w_rct,
      estimator,
      effect_measure,
      family_y,
      stabilize = TRUE,
      trim = NULL
    ),
    error = function(e) return(NULL)
  )
  if (is.null(rct_effect_b)) return(NA_real_)

  obs_effect_b <- tryCatch(
    estimate_single_effect(
      Y_obs_b,
      A_obs_b,
      X_obs_b,
      weights_b$w_obs,
      estimator,
      effect_measure,
      family_y,
      stabilize = TRUE,
      trim = NULL
    ),
    error = function(e) return(NULL)
  )
  if (is.null(obs_effect_b)) return(NA_real_)

  # Return difference
  compute_effect_difference(rct_effect_b, obs_effect_b, effect_measure)
}
