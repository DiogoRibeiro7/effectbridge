# plotting.R - Plot methods and visualization functions

#' Plot method for unconf_test objects
#' @param x An object of class `unconf_test`
#' @param type Type of plot: "effects", "overlap", "balance", "weights", "bootstrap"
#' @param ... Additional plotting arguments
#' @return Invisibly returns the plot object (if applicable)
#' @seealso [print.unconf_test()], [summary.unconf_test()],
#'   [unconfoundedness_test()] for creating the object.
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
  plot(estimates, 1:3,
    xlim = range(c(estimates, ci_lower, ci_upper), na.rm = TRUE),
    ylim = c(0.5, 3.5), yaxt = "n", ylab = "", xlab = "Effect Estimate",
    main = "Effect Estimates Comparison", pch = 19, cex = 1.5
  )

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

  hist(ps_rct,
    breaks = breaks, col = rgb(1, 0, 0, 0.3),
    main = "Propensity Score Distributions",
    xlab = "Propensity Score", ylab = "Frequency"
  )
  hist(ps_obs, breaks = breaks, col = rgb(0, 0, 1, 0.3), add = TRUE)
  legend("topright", c("RCT", "Observational"),
    fill = c(rgb(1, 0, 0, 0.3), rgb(0, 0, 1, 0.3))
  )

  # Density comparison
  plot(density(ps_rct),
    col = "red", lwd = 2,
    main = "Propensity Score Density Comparison",
    xlab = "Propensity Score", ylab = "Density"
  )
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
    main = "Covariate Balance Assessment"
  )

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
    col = c("blue", "red", "green")
  )

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
  hist(weights,
    breaks = 30, col = "lightblue", border = "black",
    main = "Distribution of Transport Weights",
    xlab = "Weight", ylab = "Frequency"
  )

  # Box plot
  boxplot(weights,
    col = "lightgreen",
    main = "Transport Weights Box Plot",
    ylab = "Weight"
  )

  # Q-Q plot against normal
  qqnorm(log(weights), main = "Q-Q Plot (Log Weights vs Normal)")
  qqline(log(weights))

  # Weight vs index (to check for patterns)
  plot(weights,
    main = "Weights by Observation Index",
    xlab = "Observation Index", ylab = "Weight",
    pch = 19, cex = 0.5, col = "darkblue"
  )

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
  hist(boot_dist,
    breaks = 30, col = "lightcoral", border = "black",
    main = "Bootstrap Distribution",
    xlab = "Difference Estimate", ylab = "Frequency"
  )
  abline(v = x$estimates$difference, col = "red", lwd = 2, lty = 2)
  abline(v = 0, col = "black", lwd = 1, lty = 2)

  # Q-Q plot
  qqnorm(boot_dist, main = "Q-Q Plot vs Normal")
  qqline(boot_dist)

  # Time series plot
  plot(boot_dist,
    type = "l", main = "Bootstrap Sequence",
    xlab = "Bootstrap Replicate", ylab = "Difference Estimate"
  )
  abline(h = x$estimates$difference, col = "red", lwd = 2, lty = 2)

  # Cumulative average
  cumavg <- cumsum(boot_dist) / seq_along(boot_dist)
  plot(cumavg,
    type = "l", main = "Cumulative Average",
    xlab = "Bootstrap Replicate", ylab = "Cumulative Average"
  )
  abline(h = x$estimates$difference, col = "red", lwd = 2, lty = 2)

  invisible(NULL)
}
