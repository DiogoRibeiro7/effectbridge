# Test suite for transport weighting and shift detection

test_that("transport weighting produces valid weights", {
  set.seed(200)

  rct_data <- generate_rct_data(n = 200, seed = 200)
  obs_data <- generate_obs_data(n = 300, covariate_shift = TRUE, seed = 300)

  result <- unconfoundedness_test(
    rct_data,
    obs_data,
    Y ~ A + X1 + X2,
    transport = "rct_to_obs",
    B = 50,
    validate = FALSE,
    seed = 400
  )

  expect_true(result$methods$transport_applied)
  expect_true(!is.null(result$diagnostics$transport))
  expect_true(result$diagnostics$transport$effective_sample_size > 0)
  expect_true(result$diagnostics$transport$coefficient_of_variation >= 0)
})

test_that("transport diagnostics have correct structure", {
  set.seed(201)

  rct_data <- generate_rct_data(n = 200, seed = 201)
  obs_data <- generate_obs_data(n = 300, covariate_shift = TRUE, seed = 301)

  result <- unconfoundedness_test(
    rct_data,
    obs_data,
    Y ~ A + X1 + X2,
    transport = "rct_to_obs",
    B = 50,
    validate = FALSE,
    seed = 401
  )

  transport_diag <- result$diagnostics$transport
  expect_true(is.list(transport_diag))
  expect_true("weight_range" %in% names(transport_diag))
  expect_true("weight_quantiles" %in% names(transport_diag))
  expect_true("effective_sample_size" %in% names(transport_diag))
  expect_true("coefficient_of_variation" %in% names(transport_diag))
  expect_equal(length(transport_diag$weight_range), 2)
  expect_true(transport_diag$weight_range[1] <= transport_diag$weight_range[2])
})

test_that("auto transport with energy method works", {
  set.seed(202)

  # Create data with clear covariate shift
  gen_rct <- function(n) {
    X1 <- rnorm(n)
    X2 <- rbinom(n, 1, 0.3)
    A <- rbinom(n, 1, 0.5)
    Y0 <- 0.5 + 0.5 * X1 + 0.3 * X2 + rnorm(n)
    Y1 <- Y0 + 1.0
    data.frame(Y = ifelse(A == 1, Y1, Y0), A, X1, X2)
  }

  gen_obs <- function(n) {
    X1 <- rnorm(n, 0.8)
    X2 <- rbinom(n, 1, 0.7)
    A <- rbinom(n, 1, plogis(-0.2 + 0.8 * X1 + 0.6 * X2))
    Y0 <- 0.5 + 0.5 * X1 + 0.3 * X2 + rnorm(n)
    Y1 <- Y0 + 1.0
    data.frame(Y = ifelse(A == 1, Y1, Y0), A, X1, X2)
  }

  d_rct <- gen_rct(300)
  d_obs <- gen_obs(400)

  # Only run energy test if package available
  skip_if_not_installed("energy")

  result <- unconfoundedness_test(
    d_rct,
    d_obs,
    Y ~ A + X1 + X2,
    estimator = "aipw",
    transport = "auto",
    auto_method = "energy",
    auto_alpha = 0.01,
    B = 50,
    validate = FALSE,
    seed = 402
  )

  # With this degree of shift, auto should detect it
  expect_true(result$diagnostics$auto_detection$shift_detected)
  expect_true(result$methods$transport_applied)
})

test_that("no transport leaves weights at 1", {
  set.seed(203)

  rct_data <- generate_rct_data(n = 200, seed = 203)
  obs_data <- generate_obs_data(n = 300, seed = 303)

  result <- unconfoundedness_test(
    rct_data,
    obs_data,
    Y ~ A + X1 + X2,
    transport = "none",
    B = 50,
    validate = FALSE,
    seed = 403
  )

  expect_false(result$methods$transport_applied)
  expect_null(result$diagnostics$transport)
})
