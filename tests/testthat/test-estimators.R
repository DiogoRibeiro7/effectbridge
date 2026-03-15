# Test suite for estimator implementations

test_that("G-computation estimator works", {
  set.seed(100)

  rct_data <- generate_rct_data(n = 200, seed = 100)
  obs_data <- generate_obs_data(n = 300, confounding_strength = 0.1, seed = 200)

  result <- unconfoundedness_test(
    rct_data,
    obs_data,
    Y ~ A + X1 + X2,
    estimator = "gcomp",
    inference_method = "analytical",
    validate = FALSE,
    seed = 300
  )

  expect_s3_class(result, "unconf_test")
  expect_equal(result$methods$estimator, "gcomp")
  expect_true(is.numeric(result$estimates$rct))
  expect_true(is.numeric(result$estimates$observational))
  expect_true(is.numeric(result$estimates$difference))
  expect_true(is.numeric(result$inference$p_value))
})

test_that("G-computation works with binary outcomes", {
  set.seed(101)

  rct_data <- generate_rct_data(n = 200, outcome_type = "binary", seed = 101)
  obs_data <- generate_obs_data(
    n = 300,
    outcome_type = "binary",
    confounding_strength = 0.1,
    seed = 201
  )

  result <- unconfoundedness_test(
    rct_data,
    obs_data,
    Y ~ A + X1 + X2,
    estimator = "gcomp",
    effect_measure = "rd",
    family_y = "binomial",
    inference_method = "analytical",
    validate = FALSE,
    seed = 301
  )

  expect_s3_class(result, "unconf_test")
  expect_equal(result$methods$estimator, "gcomp")
  expect_equal(result$methods$effect_measure, "rd")
})

test_that("IPW works with unstabilized weights", {
  set.seed(102)

  rct_data <- generate_rct_data(n = 200, seed = 102)
  obs_data <- generate_obs_data(n = 300, seed = 202)

  result <- unconfoundedness_test(
    rct_data,
    obs_data,
    Y ~ A + X1 + X2,
    estimator = "ipw",
    stabilize = FALSE,
    inference_method = "analytical",
    validate = FALSE,
    seed = 302
  )

  expect_s3_class(result, "unconf_test")
  expect_equal(result$methods$estimator, "ipw")
  expect_true(is.numeric(result$estimates$difference))
})

test_that("AIPW works with no trimming", {
  set.seed(103)

  rct_data <- generate_rct_data(n = 200, seed = 103)
  obs_data <- generate_obs_data(n = 300, seed = 203)

  result <- unconfoundedness_test(
    rct_data,
    obs_data,
    Y ~ A + X1 + X2,
    estimator = "aipw",
    trim = NULL,
    inference_method = "analytical",
    validate = FALSE,
    seed = 303
  )

  expect_s3_class(result, "unconf_test")
  expect_true(is.numeric(result$estimates$difference))
})

test_that("robust inference method works", {
  set.seed(104)

  rct_data <- generate_rct_data(n = 200, seed = 104)
  obs_data <- generate_obs_data(n = 300, seed = 204)

  result <- unconfoundedness_test(
    rct_data,
    obs_data,
    Y ~ A + X1 + X2,
    estimator = "aipw",
    inference_method = "robust",
    validate = FALSE,
    seed = 304
  )

  expect_s3_class(result, "unconf_test")
  expect_equal(result$methods$inference_method, "robust")
  expect_true(is.numeric(result$inference$p_value))
  expect_true(is.numeric(result$inference$standard_error))
})

test_that("effect difference computed correctly for risk difference", {
  set.seed(105)

  rct_data <- generate_rct_data(n = 200, treatment_effect = 1.0, seed = 105)
  obs_data <- generate_obs_data(
    n = 300,
    treatment_effect = 1.0,
    confounding_strength = 0,
    seed = 205
  )

  result <- unconfoundedness_test(
    rct_data,
    obs_data,
    Y ~ A + X1 + X2,
    estimator = "aipw",
    effect_measure = "rd",
    inference_method = "analytical",
    validate = FALSE,
    seed = 305
  )

  # Difference should be obs - rct
  expected_diff <- result$estimates$observational - result$estimates$rct
  expect_equal(result$estimates$difference, expected_diff)
})

test_that("positivity-edge weights remain finite and nonnegative", {
  set.seed(106)

  n <- 400
  A <- rbinom(n, 1, 0.5)
  e <- c(
    rep(1e-3, n / 4),
    rep(0.999, n / 4),
    rep(0.05, n / 4),
    rep(0.95, n / 4)
  )

  ipw_weights <- effectbridge:::compute_ipw_weights(
    A = A,
    e = e,
    weights = rep(1, n),
    stabilize = FALSE
  )

  expect_true(all(is.finite(ipw_weights)))
  expect_true(all(ipw_weights >= 0))
  expect_gt(max(ipw_weights), 100)
})

test_that("perfect-separation-like propensity model falls back safely", {
  set.seed(107)

  n <- 300
  X1 <- rnorm(n)
  X2 <- rbinom(n, 1, 0.5)
  A <- as.integer(X1 > 0) # deterministic treatment assignment (separation)
  X <- cbind("(Intercept)" = 1, X1 = X1, X2 = X2)

  expect_warning(
    ps_fit <- effectbridge:::fit_propensity_model(A, X, weights = rep(1, n)),
    "Propensity score model fitting failed; using marginal treated proportion as PS"
  )
  e <- ps_fit$fitted_values

  expect_true(all(e >= 1e-3))
  expect_true(all(e <= 1 - 1e-3))
  expect_equal(mean(e), mean(A), tolerance = 1e-12)
  expect_equal(stats::sd(e), 0, tolerance = 1e-12)

  ipw_weights <- effectbridge:::compute_ipw_weights(
    A = A,
    e = e,
    weights = rep(1, n),
    stabilize = FALSE
  )
  expect_true(all(is.finite(ipw_weights)))
  expect_true(all(ipw_weights >= 0))
})
