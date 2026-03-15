# Test suite for S3 methods and plotting

test_that("print method produces expected sections", {
  set.seed(300)

  rct_data <- generate_rct_data(n = 150, seed = 300)
  obs_data <- generate_obs_data(n = 200, seed = 400)

  result <- unconfoundedness_test(
    rct_data,
    obs_data,
    Y ~ A + X1 + X2,
    B = 50,
    validate = FALSE,
    seed = 500
  )

  output <- capture.output(print(result))
  output_text <- paste(output, collapse = "\n")

  expect_true(grepl("Unconfoundedness Test Results", output_text))
  expect_true(grepl("Methods", output_text))
  expect_true(grepl("Results", output_text))
  expect_true(grepl("Interpretation", output_text))
  expect_true(grepl("Key Diagnostics", output_text))
})

test_that("print method with custom digits works", {
  set.seed(301)

  rct_data <- generate_rct_data(n = 150, seed = 301)
  obs_data <- generate_obs_data(n = 200, seed = 401)

  result <- unconfoundedness_test(
    rct_data,
    obs_data,
    Y ~ A + X1 + X2,
    B = 50,
    validate = FALSE,
    seed = 501
  )

  # Should not error with different digit counts
  expect_output(print(result, digits = 2))
  expect_output(print(result, digits = 6))
})

test_that("summary method returns correct class and components", {
  set.seed(302)

  rct_data <- generate_rct_data(n = 150, seed = 302)
  obs_data <- generate_obs_data(n = 200, seed = 402)

  result <- unconfoundedness_test(
    rct_data,
    obs_data,
    Y ~ A + X1 + X2,
    B = 50,
    validate = FALSE,
    seed = 502
  )

  summ <- summary(result)
  expect_s3_class(summ, "summary.unconf_test")

  expect_true("call" %in% names(summ))
  expect_true("methods" %in% names(summ))
  expect_true("estimates" %in% names(summ))
  expect_true("inference" %in% names(summ))
  expect_true("diagnostics" %in% names(summ))
  expect_true("sensitivity" %in% names(summ))
  expect_true("recommendations" %in% names(summ))
})

test_that("summary print method produces output", {
  set.seed(303)

  rct_data <- generate_rct_data(n = 150, seed = 303)
  obs_data <- generate_obs_data(n = 200, seed = 403)

  result <- unconfoundedness_test(
    rct_data,
    obs_data,
    Y ~ A + X1 + X2,
    B = 50,
    validate = FALSE,
    seed = 503
  )

  summ <- summary(result)
  output <- capture.output(print(summ))
  output_text <- paste(output, collapse = "\n")

  expect_true(grepl("Comprehensive Unconfoundedness Test Summary", output_text))
  expect_true(grepl("Detailed Results", output_text))
  expect_true(grepl("Diagnostic Summary", output_text))
})

test_that("recommendations are generated based on diagnostics", {
  set.seed(304)

  rct_data <- generate_rct_data(n = 150, seed = 304)
  obs_data <- generate_obs_data(n = 200, seed = 404)

  result <- unconfoundedness_test(
    rct_data,
    obs_data,
    Y ~ A + X1 + X2,
    estimator = "ipw",
    B = 50,
    validate = FALSE,
    seed = 504
  )

  summ <- summary(result)

  # IPW should trigger recommendation to use AIPW
  expect_true(any(grepl("AIPW", summ$recommendations)))
  # All results should get either significance or non-significance recommendation
  expect_true(length(summ$recommendations) > 0)
})

test_that("plot effects method runs without error", {
  set.seed(305)

  rct_data <- generate_rct_data(n = 150, seed = 305)
  obs_data <- generate_obs_data(n = 200, seed = 405)

  result <- unconfoundedness_test(
    rct_data,
    obs_data,
    Y ~ A + X1 + X2,
    B = 50,
    validate = FALSE,
    seed = 505
  )

  expect_no_error(with_plot_sandbox(plot(result, type = "effects")))
})

test_that("plot overlap method runs without error", {
  set.seed(306)

  rct_data <- generate_rct_data(n = 150, seed = 306)
  obs_data <- generate_obs_data(n = 200, seed = 406)

  result <- unconfoundedness_test(
    rct_data,
    obs_data,
    Y ~ A + X1 + X2,
    B = 50,
    validate = FALSE,
    seed = 506
  )

  expect_no_error(with_plot_sandbox(plot(result, type = "overlap")))
})

test_that("plot balance method runs without error", {
  set.seed(307)

  rct_data <- generate_rct_data(n = 150, seed = 307)
  obs_data <- generate_obs_data(n = 200, seed = 407)

  result <- unconfoundedness_test(
    rct_data,
    obs_data,
    Y ~ A + X1 + X2,
    B = 50,
    validate = FALSE,
    seed = 507
  )

  expect_no_error(with_plot_sandbox(plot(result, type = "balance")))
})

test_that("plot bootstrap method runs without error", {
  set.seed(308)

  rct_data <- generate_rct_data(n = 150, seed = 308)
  obs_data <- generate_obs_data(n = 200, seed = 408)

  result <- unconfoundedness_test(
    rct_data,
    obs_data,
    Y ~ A + X1 + X2,
    inference_method = "bootstrap",
    B = 50,
    validate = FALSE,
    seed = 508
  )

  expect_no_error(with_plot_sandbox(plot(result, type = "bootstrap")))
})

test_that("plot bootstrap handles non-bootstrap inference gracefully", {
  set.seed(309)

  rct_data <- generate_rct_data(n = 150, seed = 309)
  obs_data <- generate_obs_data(n = 200, seed = 409)

  result <- unconfoundedness_test(
    rct_data,
    obs_data,
    Y ~ A + X1 + X2,
    inference_method = "analytical",
    validate = FALSE,
    seed = 509
  )

  expect_output(with_plot_sandbox(plot(result, type = "bootstrap")), "not available")
})

test_that("plot weights handles no-transport case gracefully", {
  set.seed(310)

  rct_data <- generate_rct_data(n = 150, seed = 310)
  obs_data <- generate_obs_data(n = 200, seed = 410)

  result <- unconfoundedness_test(
    rct_data,
    obs_data,
    Y ~ A + X1 + X2,
    transport = "none",
    B = 50,
    validate = FALSE,
    seed = 510
  )

  expect_output(with_plot_sandbox(plot(result, type = "weights")), "No transport")
})
