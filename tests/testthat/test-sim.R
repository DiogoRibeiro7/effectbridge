# tests/testthat/test-sim.R

test_that("data generators return expected columns", {
  r <- generate_rct_data(n = 50, seed = 1)
  o <- generate_obs_data(n = 50, seed = 1)
  expect_true(all(c("Y", "A") %in% names(r)))
  expect_true(all(c("Y", "A") %in% names(o)))
})

test_that("unconfoundedness_test runs with AIPW small n", {
  r <- generate_rct_data(n = 80, seed = 1)
  o <- generate_obs_data(n = 120, seed = 2)
  res <- unconfoundedness_test(
    data_rct = r, data_obs = o, formula = Y ~ A + X1 + X2,
    estimator = "aipw", B = 10, validate = FALSE
  )
  expect_s3_class(res, "unconf_test")
  expect_true(is.list(res$raw_effects))
})

test_that("generate_rct_data respects sample size and covariate count", {
  r <- generate_rct_data(n = 200, p_covariates = 5, seed = 10)
  expect_equal(nrow(r), 200)
  expect_true(all(paste0("X", 1:5) %in% names(r)))
  expect_equal(ncol(r), 7) # Y, A, X1..X5
})

test_that("generate_obs_data respects sample size and covariate count", {
  o <- generate_obs_data(n = 300, p_covariates = 4, seed = 20)
  expect_equal(nrow(o), 300)
  expect_true(all(paste0("X", 1:4) %in% names(o)))
  expect_equal(ncol(o), 6) # Y, A, X1..X4
})

test_that("treatment assignment is binary in generated data", {
  r <- generate_rct_data(n = 100, seed = 30)
  o <- generate_obs_data(n = 100, seed = 31)
  expect_true(all(r$A %in% c(0, 1)))
  expect_true(all(o$A %in% c(0, 1)))
})

test_that("RCT treatment is approximately balanced", {
  r <- generate_rct_data(n = 1000, seed = 40)
  # RCT assigns treatment with p = 0.5
  expect_true(abs(mean(r$A) - 0.5) < 0.1)
})

test_that("binary outcome type produces 0/1 outcomes", {
  r <- generate_rct_data(n = 200, outcome_type = "binary", seed = 50)
  o <- generate_obs_data(n = 200, outcome_type = "binary", seed = 51)
  expect_true(all(r$Y %in% c(0, 1)))
  expect_true(all(o$Y %in% c(0, 1)))
})

test_that("seed produces reproducible data", {
  r1 <- generate_rct_data(n = 50, seed = 99)
  r2 <- generate_rct_data(n = 50, seed = 99)
  expect_identical(r1, r2)

  o1 <- generate_obs_data(n = 50, seed = 99)
  o2 <- generate_obs_data(n = 50, seed = 99)
  expect_identical(o1, o2)
})

test_that("covariate_shift changes the covariate distribution", {
  o_no_shift <- generate_obs_data(n = 500, covariate_shift = FALSE, seed = 60)
  o_shift <- generate_obs_data(n = 500, covariate_shift = TRUE, seed = 60)
  # Shifted covariates should have higher means
  expect_true(mean(o_shift$X1) > mean(o_no_shift$X1))
})

test_that("confounding_strength = 0 produces near-random treatment", {

  o <- generate_obs_data(n = 1000, confounding_strength = 0, seed = 70)
  # With no confounding, treatment prevalence should be near baseline
  expect_true(abs(mean(o$A) - 0.38) < 0.1) # logistic(-0.5) ~ 0.38
})
