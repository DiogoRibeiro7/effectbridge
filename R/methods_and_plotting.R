# methods_and_plotting.R - Print, plot, summary, and other S3 methods

#' Enhanced print method for unconf_test objects
#' @param x An object of class `unconf_test`
#' @param digits Number of digits to display
#' @param ... Additional arguments (unused)
#' @return Invisibly returns `x`
#' @export
print.unconf_test <- function(x, digits = 4, ...) {
  
  cat("# Unconfoundedness Test Results\n")
  cat("=" * 50, "\n\n")
  
  # Method summary
  cat("## Methods\n")
  cat(sprintf("- Estimator: %s\n", x$methods$estimator))
  cat(sprintf("- Effect measure: %s\n", x$methods$effect_measure))
  cat(sprintf("- Outcome family: %s\n", x$methods$outcome_family))
  cat(sprintf("- Inference: %s\n", x$methods$inference_method))
  cat(sprintf("- Transport: %s (applied: %s)\n", 
              x$methods$transport_method, x$methods$transport_applied))
  cat("\n")
  
  # Main results
  cat("## Results\n")
  cat(sprintf("- RCT effect:          %.*f\n", digits, x$estimates$rct))
  cat(sprintf("- Observational effect: %.*f\n", digits, x$estimates$observational))
  cat(sprintf("- Difference:          %.*f\n", digits, x$estimates$difference))
  cat(sprintf("- 95%% CI for difference: [%.*f, %.*f]\n", 
              digits, x$inference$confidence_interval[1], 
              digits, x$inference$confidence_interval[2]))
  cat(sprintf("- Test statistic:      %.*f\n", digits, x$inference$test_statistic))
  cat(sprintf("- P-value:             %.*g\n", digits, x$inference$p_value))
  cat("\n")
  
  # Interpretation
  cat("## Interpretation\n")
  interpretation <- interpret_test_result(x$inference$p_value, x$estimates$difference)
  cat(sprintf("- %s\n", interpretation))
  cat("\n")
  
  # Key diagnostics
  cat("## Key Diagnostics\n")
  
  # Overlap quality
  cat(sprintf("- Overlap quality: %s\n", x$diagnostics$overlap$overlap_quality))
  
  # Balance quality  
  if (!is.null(x$diagnostics$balance$overall_quality)) {
    cat(sprintf("- Balance quality: %s\n", x$diagnostics$balance$overall_quality))
  }
  
  # Sample sizes
  cat(sprintf("- Sample sizes: RCT=%d, Obs=%d\n", 
              x$diagnostics$sample_size$raw_n$rct,
              x$diagnostics$sample_size$raw_n$obs))
  
  # E-value
  if (!is.null(x$sensitivity$e_value)) {
    cat(sprintf("- E-value: %.*f (%s)\n", 
                digits, x$sensitivity$e_value$e_value,
                x$sensitivity$e_value$interpretation))
  }
  
  cat("\n")
  
  # Warnings/recommendations
  warnings <- generate_warnings(x)
  if (length(warnings) > 0) {
    cat("## Warnings & Recommendations\n")
    for (w in warnings) {
      cat(sprintf("⚠ %s\n", w))
    }
    cat("\n")
  }
  
  cat("Use plot(x) for diagnostic plots and summary(x) for detailed output.\n")
  
  invisible(x)
}

#' Enhanced summary method for unconf_test objects
#' @param object An object of class `unconf_test`
#' @param ... Additional arguments (unused)
#' @return Summary object
#' @export
summary.unconf_test <- function(object, ...) {
  
  # Create comprehensive summary
  summary_obj <- list(
    call = object$call,
    methods = object$methods,
    estimates = object$estimates,
    inference = object$inference,
    diagnostics = create_diagnostic_summary(object$diagnostics),
    sensitivity = object$sensitivity,
    recommendations = generate_recommendations(object)
  )
  
  class(summary_obj) <- "summary.unconf_test"
  summary_obj
}

#' Print method for summary.unconf_test
#' @param x Summary object
#' @param ... Additional arguments (unused)
#' @export
print.summary.unconf_test <- function(x, ...) {
  
  cat("# Comprehensive Unconfoundedness Test Summary\n")
  cat("=" * 60, "\n\n")
  
  # Call
  cat("## Call\n")
  print(x$call)
  cat("\n")
  
  # Detailed results table
  cat("## Detailed Results\n")
  print(format_summary_table(x), row.names = FALSE)
  cat("\n")
  
  # Diagnostic summary
  cat("## Diagnostic Summary\n")
  print_diagnostic_summary(x$diagnostics)
  cat("\n")
  
  # Sensitivity analysis
  if (!is.null(x$sensitivity)) {
    cat("## Sensitivity Analysis\n")
    print_sensitivity_summary(x$sensitivity)
    cat("\n")
  }
  
  # Recommendations
  if (length(x$recommendations) > 0) {
    cat("## Recommendations\n")
    for (i in seq_along(x$recommendations)) {
      cat(sprintf("%d. %s\n", i, x$recommendations[i]))
    }
  }
  
  invisible(x)
}

#' Plot method for unconf_test objects
#' @param x An object of class `unconf_test`
#' @param type Type of plot: "effects", "overlap", "balance", "weights", "bootstrap"
#' @param ... Additional plotting arguments
#' @return Invisibly returns the plot object (if applicable)
#' @export
plot.unconf_test <- function(x, type = c("effects", "overlap", "balance", "weights", "bootstrap"), ...) {
  
  type <- match.arg(type)
  
  switch(type,
    "effects" = plot_effects(x, ...),
    "overlap" = plot_overlap(x, ...),
    "balance" = plot_balance(x, ...),
    "weights" = plot_weights(x, ...),
    "bootstrap" = plot_bootstrap(x, ...)
  )
}

#' Plot effect estimates with confidence intervals
#' @param x unconf_test object
#' @param ... Additional arguments
plot_effects <- function(x, ...) {
  
  # Prepare data for plotting
  estimates <- c(x$estimates$rct, x$estimates$observational, x$estimates$difference)
  labels <- c("RCT", "Observational", "Difference")
  
  # Get confidence intervals (approximate for RCT and Obs)
  ci_lower <- c(NA, NA, x$inference$confidence_interval[1])
  ci_upper <- c(NA, NA, x$inference$confidence_interval[2])
  
  # Create plot
  par(mfrow = c(1, 1), mar = c(5, 8, 4, 2))
  
  # Plot points
  plot(estimates, 1:3, xlim = range(c(estimates, ci_lower, ci_upper), na.rm = TRUE),
       ylim = c(0.5, 3.5), yaxt = "n", ylab = "", xlab = "Effect Estimate",
       main = "Effect Estimates Comparison", pch = 19, cex = 1.5)
  
  # Add confidence interval for difference
  if (!is.na(ci_lower[3])) {
    arrows(ci_lower[3], 3, ci_upper[3], 3, angle = 90, code = 3, length = 0.1, lwd = 2)
  }
  
  # Add reference line at 0 for difference
  abline(v = 0, lty = 2, col = "gray")
  
  # Add labels
  axis(2, at = 1:3, labels = labels, las = 2)
  
  # Add effect values as text
  text(estimates, 1:3, sprintf("%.3f", estimates), pos = 4, offset = 0.5)
  
  invisible(NULL)
}

#' Plot propensity score overlap
#' @param x unconf_test object
#' @param ... Additional arguments
plot_overlap <- function(x, ...) {
  
  # Extract propensity scores
  ps_rct <- x$raw_effects$rct_effect$propensity_scores
  ps_obs <- x$raw_effects$obs_effect$propensity_scores
  
  if (is.null(ps_rct) || is.null(ps_obs)) {
    cat("Propensity scores not available for overlap plot.\n")
    return(invisible(NULL))
  }
  
  par(mfrow = c(2, 1), mar = c(4, 4, 3, 2))
  
  # Histogram comparison
  hist_range <- range(c(ps_rct, ps_obs))
  breaks <- seq(hist_range[1], hist_range[2], length.out = 20)
  
  hist(ps_rct, breaks = breaks, col = rgb(1, 0, 0, 0.3), 
       main = "Propensity Score Distributions", 
       xlab = "Propensity Score", ylab = "Frequency")
  hist(ps_obs, breaks = breaks, col = rgb(0, 0, 1, 0.3), add = TRUE)
  legend("topright", c("RCT", "Observational"), 
         fill = c(rgb(1, 0, 0, 0.3), rgb(0, 0, 1, 0.3)))
  
  # Density comparison
  plot(density(ps_rct), col = "red", lwd = 2, 
       main = "Propensity Score Density Comparison",
       xlab = "Propensity Score", ylab = "Density")
  lines(density(ps_obs), col = "blue", lwd = 2)
  legend("topright", c("RCT", "Observational"), col = c("red", "blue"), lwd = 2)
  
  invisible(NULL)
}

#' Plot covariate balance
#' @param x unconf_test object
#' @param ... Additional arguments
plot_balance <- function(x, ...) {
  
  balance_info <- x$diagnostics$balance
  
  if (is.null(balance_info) || is.null(balance_info$between_datasets)) {
    cat("Balance information not available.\n")
    return(invisible(NULL))
  }
  
  # Extract standardized mean differences
  smd_between <- balance_info$between_datasets$smd
  smd_rct <- balance_info$within_rct$smd
  smd_obs <- balance_info$within_obs$smd
  
  if (length(smd_between) == 0) {
    cat("No covariates available for balance plot.\n")
    return(invisible(NULL))
  }
  
  # Create balance plot (Love plot style)
  par(mfrow = c(1, 1), mar = c(5, 8, 4, 2))
  
  var_names <- names(smd_between)
  n_vars <- length(var_names)
  
  plot(range(c(smd_between, smd_rct, smd_obs)), c(0.5, n_vars + 0.5),
       type = "n", yaxt = "n", ylab = "", xlab = "Standardized Mean Difference",
       main = "Covariate Balance Assessment")
  
  # Add reference lines
  abline(v = c(-0.1, 0, 0.1), lty = c(2, 1, 2), col = c("red", "black", "red"))
  
  # Plot points
  points(smd_between, 1:n_vars, pch = 19, col = "blue", cex = 1.2)
  
  if (length(smd_rct) > 0) {
    points(smd_rct, (1:n_vars) - 0.1, pch = 17, col = "red", cex = 1)
  }
  
  if (length(smd_obs) > 0) {
    points(smd_obs, (1:n_vars) + 0.1, pch = 15, col = "green", cex = 1)
  }
  
  # Add variable names
  axis(2, at = 1:n_vars, labels = var_names, las = 2)
  
  # Add legend
  legend("topright", 
         c("Between datasets", "Within RCT", "Within Obs"),
         pch = c(19, 17, 15), 
         col = c("blue", "red", "green"))
  
  invisible(NULL)
}

#' Plot transport weights distribution
#' @param x unconf_test object  
#' @param ... Additional arguments
plot_weights <- function(x, ...) {
  
  if (!x$methods$transport_applied) {
    cat("No transport weights applied.\n")
    return(invisible(NULL))
  }
  
  # Extract weights from raw effects
  weights <- x$raw_effects$rct_effect$weights %||% 
            attr(x$raw_effects$rct_effect, "weights")
  
  if (is.null(weights)) {
    cat("Weight information not available.\n")
    return(invisible(NULL))
  }
  
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 2))
  
  # Histogram of weights
  hist(weights, breaks = 30, col = "lightblue", border = "black",
       main = "Distribution of Transport Weights", 
       xlab = "Weight", ylab = "Frequency")
  
  # Box plot
  boxplot(weights, col = "lightgreen", 
          main = "Transport Weights Box Plot",
          ylab = "Weight")
  
  # Q-Q plot against normal
  qqnorm(log(weights), main = "Q-Q Plot (Log Weights vs Normal)")
  qqline(log(weights))
  
  # Weight vs index (to check for patterns)
  plot(weights, main = "Weights by Observation Index",
       xlab = "Observation Index", ylab = "Weight", 
       pch = 19, cex = 0.5, col = "darkblue")
  
  invisible(NULL)
}

#' Plot bootstrap distribution
#' @param x unconf_test object
#' @param ... Additional arguments  
plot_bootstrap <- function(x, ...) {
  
  if (x$methods$inference_method != "bootstrap" || 
      is.null(x$inference$bootstrap_distribution)) {
    cat("Bootstrap distribution not available.\n")
    return(invisible(NULL))
  }
  
  boot_dist <- x$inference$bootstrap_distribution
  
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 2))
  
  # Histogram
  hist(boot_dist, breaks = 30, col = "lightcoral", border = "black",
       main = "Bootstrap Distribution", 
       xlab = "Difference Estimate", ylab = "Frequency")
  abline(v = x$estimates$difference, col = "red", lwd = 2, lty = 2)
  abline(v = 0, col = "black", lwd = 1, lty = 2)
  
  # Q-Q plot
  qqnorm(boot_dist, main = "Q-Q Plot vs Normal")
  qqline(boot_dist)
  
  # Time series plot
  plot(boot_dist, type = "l", main = "Bootstrap Sequence",
       xlab = "Bootstrap Replicate", ylab = "Difference Estimate")
  abline(h = x$estimates$difference, col = "red", lwd = 2, lty = 2)
  
  # Cumulative average
  cumavg <- cumsum(boot_dist) / seq_along(boot_dist)
  plot(cumavg, type = "l", main = "Cumulative Average",
       xlab = "Bootstrap Replicate", ylab = "Cumulative Average")
  abline(h = x$estimates$difference, col = "red", lwd = 2, lty = 2)
  
  invisible(NULL)
}

# Helper functions for methods

#' Interpret test result
#' @param p_value P-value from test
#' @param difference Estimated difference
#' @return Interpretation string
interpret_test_result <- function(p_value, difference) {
  
  significance <- if (p_value < 0.001) {
    "highly significant"
  } else if (p_value < 0.01) {
    "significant"
  } else if (p_value < 0.05) {
    "marginally significant"
  } else {
    "non-significant"
  }
  
  direction <- if (abs(difference) < 1e-6) {
    "no evidence of"
  } else if (difference > 0) {
    "evidence of positive"
  } else {
    "evidence of negative"
  }
  
  paste0("Results show ", direction, " confounding (", significance, 
         " at p = ", sprintf("%.3f", p_value), ")")
}

#' Generate warnings based on diagnostics
#' @param x unconf_test object
#' @return Character vector of warnings
generate_warnings <- function(x) {
  
  warnings <- character(0)
  
  # Check overlap
  if (x$diagnostics$overlap$overlap_quality %in% c("Poor", "Fair")) {
    warnings <- c(warnings, 
                  paste("Poor propensity score overlap detected.",
                        "Consider restricting analysis to region of common support."))
  }
  
  # Check balance
  if (!is.null(x$diagnostics$balance$overall_quality) && 
      x$diagnostics$balance$overall_quality == "Poor") {
    warnings <- c(warnings, 
                  "Poor covariate balance. Results may be sensitive to model specification.")
  }
  
  # Check sample size
  if (x$diagnostics$sample_size$minimum_group_size < 30) {
    warnings <- c(warnings, 
                  "Very small sample size in some treatment groups. Results may be unstable.")
  }
  
  # Check effective sample size from transport
  if (!is.null(x$diagnostics$transport) && 
      x$diagnostics$transport$effective_sample_size < 0.3 * x$diagnostics$sample_size$raw_n$rct) {
    warnings <- c(warnings, 
                  "Transport weighting substantially reduces effective sample size.")
  }
  
  # Check E-value
  if (!is.null(x$sensitivity$e_value) && x$sensitivity$e_value$e_value < 2.0) {
    warnings <- c(warnings, 
                  "Low E-value suggests results may be fragile to unmeasured confounding.")
  }
  
  warnings
}

#' Create diagnostic summary for summary method
#' @param diagnostics Diagnostics object
#' @return Summarized diagnostics
create_diagnostic_summary <- function(diagnostics) {
  
  summary_diag <- list()
  
  # Overlap summary
  summary_diag$overlap <- list(
    quality = diagnostics$overlap$overlap_quality,
    rct_extreme = diagnostics$overlap$extreme_ps_proportion$rct,
    obs_extreme = diagnostics$overlap$extreme_ps_proportion$obs,
    treatment_prevalence_diff = abs(diff(unlist(diagnostics$overlap$treatment_prevalence)))
  )
  
  # Balance summary
  if (!is.null(diagnostics$balance)) {
    summary_diag$balance <- list(
      overall_quality = diagnostics$balance$overall_quality,
      max_smd_between = diagnostics$balance$between_datasets$max_smd,
      max_smd_within_rct = diagnostics$balance$within_rct$max_smd,
      max_smd_within_obs = diagnostics$balance$within_obs$max_smd
    )
  }
  
  # Sample size summary
  summary_diag$sample_size <- list(
    total_n = diagnostics$sample_size$raw_n$rct + diagnostics$sample_size$raw_n$obs,
    min_group_size = diagnostics$sample_size$minimum_group_size,
    rct_efficiency = diagnostics$sample_size$weight_efficiency$rct,
    obs_efficiency = diagnostics$sample_size$weight_efficiency$obs
  )
  
  # Transport summary
  if (!is.null(diagnostics$transport)) {
    summary_diag$transport <- list(
      effective_sample_size = diagnostics$transport$effective_sample_size,
      weight_cv = diagnostics$transport$coefficient_of_variation,
      extreme_weights = sum(diagnostics$transport$weight_range)
    )
  }
  
  summary_diag
}

#' Format summary table
#' @param summary_obj Summary object
#' @return Formatted data frame
format_summary_table <- function(summary_obj) {
  
  data.frame(
    Estimate = c(
      sprintf("%.4f", summary_obj$estimates$rct),
      sprintf("%.4f", summary_obj$estimates$observational),
      sprintf("%.4f", summary_obj$estimates$difference)
    ),
    SE = c(
      "—", "—", 
      sprintf("%.4f", summary_obj$inference$standard_error)
    ),
    CI_95 = c(
      "—", "—",
      sprintf("[%.4f, %.4f]", 
              summary_obj$inference$confidence_interval[1],
              summary_obj$inference$confidence_interval[2])
    ),
    P_value = c(
      "—", "—",
      format_p_value(summary_obj$inference$p_value)
    ),
    row.names = c("RCT Effect", "Observational Effect", "Difference")
  )
}

#' Format p-value for display
#' @param p P-value
#' @return Formatted p-value string
format_p_value <- function(p) {
  if (p < 0.001) {
    "< 0.001"
  } else {
    sprintf("%.3f", p)
  }
}

#' Print diagnostic summary
#' @param diag_summary Diagnostic summary object
print_diagnostic_summary <- function(diag_summary) {
  
  cat(sprintf("- Overlap quality: %s\n", diag_summary$overlap$quality))
  cat(sprintf("- Extreme PS proportion: RCT=%.1f%%, Obs=%.1f%%\n",
              100 * diag_summary$overlap$rct_extreme,
              100 * diag_summary$overlap$obs_extreme))
  
  if (!is.null(diag_summary$balance)) {
    cat(sprintf("- Balance quality: %s (max SMD: %.3f)\n",
                diag_summary$balance$overall_quality,
                diag_summary$balance$max_smd_between))
  }
  
  cat(sprintf("- Sample size: total=%d, min group=%d\n",
              diag_summary$sample_size$total_n,
              diag_summary$sample_size$min_group_size))
  
  if (!is.null(diag_summary$transport)) {
    cat(sprintf("- Transport efficiency: ESS=%.1f, Weight CV=%.2f\n",
                diag_summary$transport$effective_sample_size,
                diag_summary$transport$weight_cv))
  }
}

#' Print sensitivity summary
#' @param sensitivity Sensitivity analysis object
print_sensitivity_summary <- function(sensitivity) {
  
  if (!is.null(sensitivity$e_value)) {
    cat(sprintf("- E-value: %.2f (%s)\n",
                sensitivity$e_value$e_value,
                sensitivity$e_value$interpretation))
  }
  
  if (!is.null(sensitivity$fragility)) {
    cat(sprintf("- Fragility: %s\n", sensitivity$fragility$overall))
    if (length(sensitivity$fragility$factors) > 0) {
      cat("  Factors:", paste(sensitivity$fragility$factors, collapse = ", "), "\n")
    }
  }
  
  if (!is.null(sensitivity$robustness)) {
    cat(sprintf("- Model robustness: %s\n", sensitivity$robustness$model_specification))
    cat(sprintf("- Method robustness: %s\n", sensitivity$robustness$estimation_method))
    cat(sprintf("- Data robustness: %s\n", sensitivity$robustness$data_quality))
  }
}

#' Generate recommendations
#' @param x unconf_test object
#' @return Character vector of recommendations
generate_recommendations <- function(x) {
  
  recommendations <- character(0)
  
  # Based on overlap quality
  if (x$diagnostics$overlap$overlap_quality == "Poor") {
    recommendations <- c(recommendations,
                        "Consider trimming extreme propensity scores or using matching methods")
  }
  
  # Based on balance
  if (!is.null(x$diagnostics$balance$overall_quality) && 
      x$diagnostics$balance$overall_quality == "Poor") {
    recommendations <- c(recommendations,
                        "Include additional covariates or use more flexible models")
  }
  
  # Based on sample size
  if (x$diagnostics$sample_size$minimum_group_size < 50) {
    recommendations <- c(recommendations,
                        "Increase sample size or combine analysis with similar studies")
  }
  
  # Based on method
  if (x$methods$estimator == "ipw") {
    recommendations <- c(recommendations,
                        "Consider using AIPW (doubly robust) estimator for better protection against model misspecification")
  }
  
  # Based on significance
  if (x$inference$p_value < 0.05) {
    recommendations <- c(recommendations,
                        "Investigate potential sources of confounding and consider sensitivity analyses")
  } else {
    recommendations <- c(recommendations,
                        "Non-significant result supports unconfoundedness assumption, but consider power analysis")
  }
  
  # Based on E-value
  if (!is.null(x$sensitivity$e_value) && x$sensitivity$e_value$e_value < 2.0) {
    recommendations <- c(recommendations,
                        "Low E-value suggests vulnerability to unmeasured confounding. Consider instrumental variables or other approaches")
  }
  
  recommendations
}

# Utility operator (defined again for safety)
`%||%` <- function(a, b) if (is.null(a)) b else a
