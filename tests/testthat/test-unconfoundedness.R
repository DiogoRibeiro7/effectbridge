# Enhanced test suite for effectbridge package

# Test basic functionality
test_that("basic unconfoundedness test works", {
  set.seed(123)

  # Generate simple test data
  rct_data <- generate_rct_data(n = 200, treatment_effect = 1.0, seed = 123)
  obs_data <- generate_obs_data(
    n = 300,
    treatment_effect = 1.0,
    confounding_strength = 0,
    seed = 456
  )

  # Basic test
  result <- unconfoundedness_test(
    rct_data,
    obs_data,
    Y ~ A + X1 + X2,
    estimator = "aipw",
    B = 100, # Small B for testing
    validate = FALSE, # Skip validation for speed
    seed = 789
  )

  # Check structure
  expect_s3_class(result, "unconf_test")
  expect_true(is.list(result))
  expect_true(all(
    c("estimates", "inference", "diagnostics", "methods") %in% names(result)
  ))

  # Check estimates exist and are numeric
  expect_true(is.numeric(result$estimates$rct))
  expect_true(is.numeric(result$estimates$observational))
  expect_true(is.numeric(result$estimates$difference))

  # Check inference
  expect_true(is.numeric(result$inference$p_value))
  expect_true(length(result$inference$confidence_interval) == 2)
  expect_true(
    result$inference$confidence_interval[1] <
      result$inference$confidence_interval[2]
  )
})

# Test different estimators
test_that("different estimators work", {
  set.seed(456)

  rct_data <- generate_rct_data(n = 150, seed = 123)
  obs_data <- generate_obs_data(n = 200, confounding_strength = 0.1, seed = 456)

  estimators <- c("ipw", "aipw")

  for (est in estimators) {
    result <- unconfoundedness_test(
      rct_data,
      obs_data,
      Y ~ A + X1 + X2,
      estimator = est,
      B = 50,
      validate = FALSE,
      seed = 789
    )

    expect_s3_class(result, "unconf_test")
    expect_equal(result$methods$estimator, est)
    expect_true(is.numeric(result$estimates$rct))
  }
})

# Test effect measures
test_that("different effect measures work for binary outcomes", {
  set.seed(789)

  rct_data <- generate_rct_data(n = 200, outcome_type = "binary", seed = 123)
  obs_data <- generate_obs_data(
    n = 250,
    outcome_type = "binary",
    confounding_strength = 0.1,
    seed = 456
  )

  effect_measures <- c("rd", "rr", "or")

  for (measure in effect_measures) {
    result <- unconfoundedness_test(
      rct_data,
      obs_data,
      Y ~ A + X1 + X2,
      effect_measure = measure,
      family_y = "binomial",
      B = 50,
      validate = FALSE,
      seed = 789
    )

    expect_s3_class(result, "unconf_test")
    expect_equal(result$methods$effect_measure, measure)
    expect_true(is.numeric(result$estimates$rct))
    expect_true(is.finite(result$estimates$rct))
  }
})

# Test transport methods
test_that("transport methods work", {
  set.seed(111)

  # Create data with covariate shift
  rct_data <- generate_rct_data(n = 200, seed = 123)
  obs_data <- generate_obs_data(n = 300, covariate_shift = TRUE, seed = 456)

  transport_methods <- c("none", "rct_to_obs")

  for (transport in transport_methods) {
    result <- unconfoundedness_test(
      rct_data,
      obs_data,
      Y ~ A + X1 + X2,
      transport = transport,
      B = 50,
      validate = FALSE,
      seed = 789
    )

    expect_s3_class(result, "unconf_test")
    expect_equal(result$methods$transport_method, transport)

    if (transport == "rct_to_obs") {
      expect_true(result$methods$transport_applied)
    } else {
      expect_false(result$methods$transport_applied)
    }
  }
})

# Test auto transport detection
test_that("auto transport stays off with no covariate shift", {
  set.seed(100)

  # Generate data with same distribution
  gen <- function(n, seed_offset = 0) {
    set.seed(100 + seed_offset)
    X1 <- rnorm(n)
    X2 <- rbinom(n, 1, 0.4)
    A <- rbinom(n, 1, plogis(-0.2 + 0.8 * X1 + 0.6 * X2))
    Y0 <- 0.5 + 0.5 * X1 + 0.3 * X2 + rnorm(n)
    Y1 <- Y0 + 1.0
    data.frame(Y = ifelse(A == 1, Y1, Y0), A, X1, X2)
  }

  d_rct <- gen(300, 0)
  d_obs <- gen(400, 1000)

  result <- unconfoundedness_test(
    d_rct,
    d_obs,
    Y ~ A + X1 + X2,
    estimator = "aipw",
    family_y = "gaussian",
    transport = "auto",
    auto_method = "both",
    auto_alpha = 0.01,
    B = 100,
    validate = FALSE,
    seed = 100
  )

  expect_false(result$methods$transport_applied)
  expect_false(result$diagnostics$auto_detection$shift_detected)
})

# Test auto transport turns on with covariate shift
test_that("auto transport turns on with covariate shift", {
  set.seed(101)

  gen_rct <- function(n) {
    X1 <- rnorm(n)
    X2 <- rbinom(n, 1, 0.3)
    A <- rbinom(n, 1, 0.5)
    Y0 <- 0.5 + 0.5 * X1 + 0.3 * X2 + rnorm(n)
    Y1 <- Y0 + 1.0
    data.frame(Y = ifelse(A == 1, Y1, Y0), A, X1, X2)
  }

  gen_obs <- function(n) {
    X1 <- rnorm(n, 0.7)
    X2 <- rbinom(n, 1, 0.7) # Clear shift
    A <- rbinom(n, 1, plogis(-0.2 + 0.8 * X1 + 0.6 * X2))
    Y0 <- 0.5 + 0.5 * X1 + 0.3 * X2 + rnorm(n)
    Y1 <- Y0 + 1.0
    data.frame(Y = ifelse(A == 1, Y1, Y0), A, X1, X2)
  }

  d_rct <- gen_rct(300)
  d_obs <- gen_obs(400)

  result <- unconfoundedness_test(
    d_rct,
    d_obs,
    Y ~ A + X1 + X2,
    estimator = "aipw",
    family_y = "gaussian",
    transport = "auto",
    auto_method = "ks",
    auto_alpha = 0.01,
    B = 100,
    validate = FALSE,
    seed = 101
  )

  expect_true(result$methods$transport_applied)
  expect_true(result$diagnostics$auto_detection$shift_detected)
})

# Test inference methods
test_that("different inference methods work", {
  set.seed(222)

  rct_data <- generate_rct_data(n = 200, seed = 123)
  obs_data <- generate_obs_data(n = 250, confounding_strength = 0.1, seed = 456)

  # Bootstrap inference
  result_boot <- unconfoundedness_test(
    rct_data,
    obs_data,
    Y ~ A + X1 + X2,
    estimator = "aipw",
    inference_method = "bootstrap",
    B = 100,
    validate = FALSE,
    seed = 789
  )

  expect_equal(result_boot$methods$inference_method, "bootstrap")
  expect_true("bootstrap_distribution" %in% names(result_boot$inference))

  # Analytical inference
  result_analytical <- unconfoundedness_test(
    rct_data,
    obs_data,
    Y ~ A + X1 + X2,
    estimator = "aipw",
    inference_method = "analytical",
    validate = FALSE,
    seed = 789
  )

  expect_equal(result_analytical$methods$inference_method, "analytical")
  expect_true(is.numeric(result_analytical$inference$standard_error))
})

# Test input validation
test_that("input validation works", {
  set.seed(333)

  rct_data <- generate_rct_data(n = 100, seed = 123)
  obs_data <- generate_obs_data(n = 150, seed = 456)

  # Missing variables
  expect_error(
    unconfoundedness_test(rct_data[, -1], obs_data, Y ~ A + X1 + X2),
    "Missing variables"
  )

  # Invalid alpha
  expect_error(
    unconfoundedness_test(rct_data, obs_data, Y ~ A + X1 + X2, alpha = 1.5),
    "alpha must be between"
  )

  # Invalid trim
  expect_error(
    unconfoundedness_test(
      rct_data,
      obs_data,
      Y ~ A + X1 + X2,
      trim = c(0.9, 0.1)
    ),
    "trim must be"
  )
})

# Test print and summary methods
test_that("print and summary methods work", {
  set.seed(444)

  rct_data <- generate_rct_data(n = 150, seed = 123)
  obs_data <- generate_obs_data(n = 200, seed = 456)

  result <- unconfoundedness_test(
    rct_data,
    obs_data,
    Y ~ A + X1 + X2,
    B = 50,
    validate = FALSE,
    seed = 789
  )

  # Test print method doesn't error
  expect_output(print(result))

  # Test summary method
  summ <- summary(result)
  expect_s3_class(summ, "summary.unconf_test")
  expect_output(print(summ))
})

# Test data generation functions
test_that("data generation functions work", {
  # RCT data generation
  rct_data <- generate_rct_data(n = 100, treatment_effect = 0.5, seed = 123)
  expect_s3_class(rct_data, "data.frame")
  expect_equal(nrow(rct_data), 100)
  expect_true(all(c("Y", "A", "X1", "X2") %in% names(rct_data)))
  expect_true(all(rct_data$A %in% c(0, 1)))

  # Observational data generation
  obs_data <- generate_obs_data(n = 200, confounding_strength = 0.3, seed = 456)
  expect_s3_class(obs_data, "data.frame")
  expect_equal(nrow(obs_data), 200)
  expect_true(all(c("Y", "A", "X1", "X2") %in% names(obs_data)))
  expect_true(all(obs_data$A %in% c(0, 1)))

  # Binary outcomes
  binary_data <- generate_rct_data(n = 100, outcome_type = "binary", seed = 789)
  expect_true(all(binary_data$Y %in% c(0, 1)))
})

# Test utility functions
test_that("utility functions work", {
  # P-value adjustment
  p_vals <- c(0.01, 0.03, 0.05, 0.10)
  adj_p <- adjust_p_values(p_vals, method = "BH")
  expect_equal(length(adj_p), length(p_vals))
  expect_true(all(adj_p >= p_vals)) # Adjusted p-values should be >= original

  # Sample size calculation
  ss <- compute_required_sample_size(power = 0.8, effect_size = 0.5)
  expect_true(is.list(ss))
  expect_true(all(c("n_rct", "n_obs") %in% names(ss)))
  expect_true(all(sapply(ss[c("n_rct", "n_obs")], is.numeric)))
})

# Test error handling and edge cases
test_that("edge cases are handled properly", {
  set.seed(555)

  # Small sample sizes trigger validation warnings
  small_rct <- generate_rct_data(n = 30, seed = 123)
  small_obs <- generate_obs_data(n = 50, seed = 456)

  expect_warning(
    result <- unconfoundedness_test(
      small_rct,
      small_obs,
      Y ~ A + X1 + X2,
      B = 50,
      seed = 789
    ),
    "Small.*sample size"
  )

  # Perfect separation (all treated or all control)
  perfect_sep <- generate_rct_data(n = 30, seed = 123)
  perfect_sep$A <- rep(1, nrow(perfect_sep)) # All treated

  # Should handle gracefully or give informative error
  expect_error(
    unconfoundedness_test(
      perfect_sep,
      small_obs,
      Y ~ A + X1 + X2,
      validate = FALSE
    ),
    class = "error" # Some kind of error expected
  )
})

# Test batch analysis
test_that("batch analysis works", {
  set.seed(666)

  # Create multiple datasets
  rct_list <- list(
    study1 = generate_rct_data(n = 100, seed = 123),
    study2 = generate_rct_data(n = 120, seed = 456)
  )

  obs_list <- list(
    study1 = generate_obs_data(n = 150, seed = 789),
    study2 = generate_obs_data(n = 180, seed = 111)
  )

  batch_results <- batch_analysis(
    rct_list,
    obs_list,
    Y ~ A + X1 + X2,
    B = 50,
    validate = FALSE
  )

  expect_s3_class(batch_results, "batch_unconf_test")
  expect_equal(length(batch_results), 2)
  expect_output(print(batch_results))
})

# Test export functions
test_that("export functions work", {
  set.seed(777)

  rct_data <- generate_rct_data(n = 100, seed = 123)
  obs_data <- generate_obs_data(n = 150, seed = 456)

  result <- unconfoundedness_test(
    rct_data,
    obs_data,
    Y ~ A + X1 + X2,
    B = 50,
    validate = FALSE,
    seed = 789
  )

  # Test CSV export
  temp_csv <- tempfile(fileext = ".csv")
  expect_output(export_to_csv(result, temp_csv), "exported")
  expect_true(file.exists(temp_csv))

  # Clean up
  unlink(temp_csv)

  # Test LaTeX table creation
  expect_output(
    latex_table <- create_latex_table(result),
    "\\\\begin\\{table\\}"
  )
  expect_true(is.character(latex_table))
})

# Test installation validation
test_that("installation validation works", {
  expect_output(validation_result <- validate_installation(), "Validating")
  expect_true(is.logical(validation_result))
})

# Test power analysis simulation (quick version)
test_that("power analysis simulation works", {
  # Very small simulation for testing
  power_result <- simulate_power_analysis(
    n_sim = 10, # Very small for testing
    n_rct = 50,
    n_obs = 75,
    treatment_effect_rct = 1.0,
    treatment_effect_obs = 1.2,
    seed = 888
  )

  expect_true(is.list(power_result))
  expect_true(all(
    c("settings", "results", "raw_results") %in% names(power_result)
  ))
  expect_true(is.numeric(power_result$results$power))
  expect_true(
    power_result$results$power >= 0 && power_result$results$power <= 1
  )
})

# Test sensitivity analysis components
test_that("sensitivity analysis components work", {
  set.seed(888)

  rct_data <- generate_rct_data(n = 200, seed = 123)
  obs_data <- generate_obs_data(n = 300, confounding_strength = 0.2, seed = 456)

  result <- unconfoundedness_test(
    rct_data,
    obs_data,
    Y ~ A + X1 + X2,
    B = 100,
    validate = FALSE,
    seed = 789
  )

  # Check that sensitivity analysis was computed
  expect_true("sensitivity" %in% names(result))
  expect_true(is.list(result$sensitivity))

  # Check E-value
  if (!is.null(result$sensitivity$e_value)) {
    expect_true(is.numeric(result$sensitivity$e_value$e_value))
    expect_true(result$sensitivity$e_value$e_value > 0)
  }
})

# Test diagnostic quality assessments
test_that("diagnostic assessments work", {
  set.seed(999)

  rct_data <- generate_rct_data(n = 200, seed = 123)
  obs_data <- generate_obs_data(n = 300, seed = 456)

  result <- unconfoundedness_test(
    rct_data,
    obs_data,
    Y ~ A + X1 + X2,
    B = 50,
    validate = FALSE,
    seed = 789
  )

  # Check diagnostic structure
  expect_true("diagnostics" %in% names(result))
  diag <- result$diagnostics

  # Overlap diagnostics
  expect_true("overlap" %in% names(diag))
  expect_true("overlap_quality" %in% names(diag$overlap))
  expect_true(
    diag$overlap$overlap_quality %in% c("Poor", "Fair", "Good", "Excellent")
  )

  # Balance diagnostics
  expect_true("balance" %in% names(diag))

  # Sample size diagnostics
  expect_true("sample_size" %in% names(diag))
  expect_equal(diag$sample_size$raw_n$rct, 200)
  expect_equal(diag$sample_size$raw_n$obs, 300)
})
