# Test suite for input validation edge cases

test_that("validation catches non-data.frame inputs", {
  expect_error(
    unconfoundedness_test(
      matrix(1:10, ncol = 2),
      data.frame(Y = 1, A = 1, X1 = 1),
      Y ~ A + X1
    ),
    "data.frames"
  )
})

test_that("validation catches missing variables in obs data", {
  set.seed(500)
  rct_data <- generate_rct_data(n = 100, seed = 500)
  obs_data <- generate_obs_data(n = 150, seed = 600)

  expect_error(
    unconfoundedness_test(rct_data, obs_data[, -3], Y ~ A + X1 + X2),
    "Missing variables"
  )
})

test_that("validation catches formula with too few terms", {
  set.seed(501)
  rct_data <- generate_rct_data(n = 100, seed = 501)
  obs_data <- generate_obs_data(n = 150, seed = 601)

  expect_error(
    unconfoundedness_test(rct_data, obs_data, Y ~ 1),
    class = "error"
  )
})

test_that("validation catches RR/OR with gaussian family", {
  set.seed(502)
  rct_data <- generate_rct_data(n = 100, seed = 502)
  obs_data <- generate_obs_data(n = 150, seed = 602)

  expect_error(
    unconfoundedness_test(
      rct_data, obs_data, Y ~ A + X1 + X2,
      effect_measure = "rr",
      family_y = "gaussian"
    ),
    "binary outcomes"
  )
})

test_that("validation warns on small sample sizes", {
  set.seed(503)
  rct_data <- generate_rct_data(n = 30, seed = 503)
  obs_data <- generate_obs_data(n = 50, seed = 603)

  expect_warning(
    unconfoundedness_test(
      rct_data, obs_data, Y ~ A + X1 + X2,
      B = 50, seed = 703
    ),
    "Small"
  )
})

test_that("coerce_binary handles logical vectors", {
  result <- coerce_binary(c(TRUE, FALSE, TRUE, FALSE))
  expect_equal(result, c(1L, 0L, 1L, 0L))
})

test_that("coerce_binary handles factor variables", {
  f <- factor(c("control", "treated", "control", "treated"))
  result <- coerce_binary(f)
  expect_true(all(result %in% c(0, 1)))
})

test_that("coerce_binary handles character variables", {
  ch <- c("no", "yes", "no", "yes")
  result <- coerce_binary(ch)
  expect_true(all(result %in% c(0, 1)))
  # "yes" should map to 1 (alphabetically second)
  expect_equal(result, c(0L, 1L, 0L, 1L))
})

test_that("coerce_binary rejects non-binary numeric", {
  expect_error(coerce_binary(c(1, 2, 3)), "exactly 2")
})

test_that("coerce_binary rejects factor with 3 levels", {
  f <- factor(c("a", "b", "c"))
  expect_error(coerce_binary(f), "exactly 2")
})

test_that("validation alpha boundary cases", {
  set.seed(504)
  rct_data <- generate_rct_data(n = 100, seed = 504)
  obs_data <- generate_obs_data(n = 150, seed = 604)

  expect_error(
    unconfoundedness_test(rct_data, obs_data, Y ~ A + X1 + X2, alpha = 0),
    "alpha must be"
  )

  expect_error(
    unconfoundedness_test(rct_data, obs_data, Y ~ A + X1 + X2, alpha = 1),
    "alpha must be"
  )

  expect_error(
    unconfoundedness_test(rct_data, obs_data, Y ~ A + X1 + X2, alpha = -0.1),
    "alpha must be"
  )
})
