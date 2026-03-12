# tests/testthat/test-sim.R
test_that("data generators return expected columns", {
  r <- generate_rct_data(n = 50, seed = 1)
  o <- generate_obs_data(n = 50, seed = 1)
  expect_true(all(c("Y","A") %in% names(r)))
  expect_true(all(c("Y","A") %in% names(o)))
})

# tests/testthat/test-main.R
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
