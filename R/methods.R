# methods.R - Print, summary, and other S3 methods

#' Enhanced print method for unconf_test objects
#' @param x An object of class `unconf_test`
#' @param digits Number of digits to display
#' @param ... Additional arguments (unused)
#' @return Invisibly returns `x`
#' @export
print.unconf_test <- function(x, digits = 4, ...) {

  cat("# Unconfoundedness Test Results\n")
  cat(strrep("=", 50), "\n\n")

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
      cat(sprintf("! %s\n", w))
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
  cat(strrep("=", 60), "\n\n")

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
  if (!is.null(x$diagnostics$overlap$overlap_quality) &&
      x$diagnostics$overlap$overlap_quality %in% c("Poor", "Fair")) {
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
  if (!is.null(x$sensitivity$e_value) &&
      !is.na(x$sensitivity$e_value$e_value) &&
      x$sensitivity$e_value$e_value < 2.0) {
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
      "---", "---",
      sprintf("%.4f", summary_obj$inference$standard_error)
    ),
    CI_95 = c(
      "---", "---",
      sprintf("[%.4f, %.4f]",
              summary_obj$inference$confidence_interval[1],
              summary_obj$inference$confidence_interval[2])
    ),
    P_value = c(
      "---", "---",
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
  if (!is.null(x$sensitivity$e_value) &&
      !is.na(x$sensitivity$e_value$e_value) &&
      x$sensitivity$e_value$e_value < 2.0) {
    recommendations <- c(recommendations,
                        "Low E-value suggests vulnerability to unmeasured confounding. Consider instrumental variables or other approaches")
  }

  recommendations
}
