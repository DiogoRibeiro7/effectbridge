# Test suite for utility functions, data generation, and exports

test_that("generate_rct_data with custom covariates works", {
  data <- generate_rct_data(n = 100, p_covariates = 5, seed = 1)
  expect_equal(nrow(data), 100)
  expect_true(all(c("Y", "A", "X1", "X2", "X3", "X4", "X5") %in% names(data)))
  expect_true(all(data$A %in% c(0, 1)))
})

test_that("generate_obs_data with covariate shift works", {
  data_no_shift <- generate_obs_data(n = 500, covariate_shift = FALSE, seed = 1)
  data_shift <- generate_obs_data(n = 500, covariate_shift = TRUE, seed = 1)

  # Shifted data should have different X1 mean
  expect_true(abs(mean(data_shift$X1) - mean(data_no_shift$X1)) > 0.1)
})

test_that("generate_obs_data binary outcomes work", {
  data <- generate_obs_data(
    n = 200,
    outcome_type = "binary",
    confounding_strength = 0.3,
    seed = 2
  )
  expect_true(all(data$Y %in% c(0, 1)))
  expect_true(all(data$A %in% c(0, 1)))
})

test_that("adjust_p_values works with all methods", {
  p_vals <- c(0.01, 0.03, 0.05, 0.10, 0.50)

  for (method in c("bonferroni", "holm", "BH", "BY")) {
    adj_p <- adjust_p_values(p_vals, method = method)
    expect_equal(length(adj_p), length(p_vals))
    expect_true(all(adj_p >= p_vals))
    expect_true(all(adj_p <= 1))
  }
})

test_that("compute_required_sample_size returns valid results", {
  ss <- compute_required_sample_size(
    power = 0.8,
    alpha = 0.05,
    effect_size = 0.3,
    ratio_obs_to_rct = 3
  )

  expect_true(is.list(ss))
  expect_true(ss$n_rct > 0)
  expect_true(ss$n_obs > 0)
  expect_true(ss$total_n == ss$n_rct + ss$n_obs)
  # Ratio is approximate due to ceiling
  expect_true(abs(ss$n_obs / ss$n_rct - 3) < 0.1)
})

test_that("compute_required_sample_size scales with power", {
  ss_low <- compute_required_sample_size(power = 0.7, effect_size = 0.5)
  ss_high <- compute_required_sample_size(power = 0.9, effect_size = 0.5)

  expect_true(ss_high$n_rct > ss_low$n_rct)
})

test_that("compute_required_sample_size scales with effect size", {
  ss_small <- compute_required_sample_size(effect_size = 0.2)
  ss_large <- compute_required_sample_size(effect_size = 0.8)

  expect_true(ss_small$n_rct > ss_large$n_rct)
})

test_that("export_to_csv works with single result", {
  set.seed(400)

  rct_data <- generate_rct_data(n = 100, seed = 400)
  obs_data <- generate_obs_data(n = 150, seed = 500)

  result <- unconfoundedness_test(
    rct_data, obs_data, Y ~ A + X1 + X2,
    B = 50, validate = FALSE, seed = 600
  )

  temp_csv <- tempfile(fileext = ".csv")
  on.exit(unlink(temp_csv), add = TRUE)

  expect_output(export_to_csv(result, temp_csv), "exported")
  expect_true(file.exists(temp_csv))

  # Read back and validate
  exported <- read.csv(temp_csv)
  expect_equal(nrow(exported), 3)
  expect_true("Estimate" %in% names(exported))
})

test_that("export_to_csv works with list of results", {
  set.seed(401)

  rct1 <- generate_rct_data(n = 100, seed = 401)
  obs1 <- generate_obs_data(n = 150, seed = 501)
  rct2 <- generate_rct_data(n = 100, seed = 402)
  obs2 <- generate_obs_data(n = 150, seed = 502)

  r1 <- unconfoundedness_test(
    rct1, obs1, Y ~ A + X1 + X2, B = 50, validate = FALSE, seed = 601
  )
  r2 <- unconfoundedness_test(
    rct2, obs2, Y ~ A + X1 + X2, B = 50, validate = FALSE, seed = 602
  )

  results <- list(study1 = r1, study2 = r2)
  temp_csv <- tempfile(fileext = ".csv")
  on.exit(unlink(temp_csv), add = TRUE)

  expect_output(export_to_csv(results, temp_csv), "exported")
  exported <- read.csv(temp_csv)
  expect_equal(nrow(exported), 6) # 3 rows per result
})

test_that("export_to_csv rejects invalid input", {
  expect_error(export_to_csv("not a result", tempfile()), "must be")
})

test_that("create_latex_table works with single result", {
  set.seed(402)

  rct_data <- generate_rct_data(n = 100, seed = 402)
  obs_data <- generate_obs_data(n = 150, seed = 502)

  result <- unconfoundedness_test(
    rct_data, obs_data, Y ~ A + X1 + X2,
    B = 50, validate = FALSE, seed = 602
  )

  output <- capture.output(latex <- create_latex_table(result))
  expect_true(is.character(latex))
  expect_true(any(grepl("begin\\{table\\}", latex)))
  expect_true(any(grepl("end\\{table\\}", latex)))
})

test_that("create_latex_table works with multiple results", {
  set.seed(403)

  rct1 <- generate_rct_data(n = 100, seed = 403)
  obs1 <- generate_obs_data(n = 150, seed = 503)
  rct2 <- generate_rct_data(n = 100, seed = 404)
  obs2 <- generate_obs_data(n = 150, seed = 504)

  r1 <- unconfoundedness_test(
    rct1, obs1, Y ~ A + X1 + X2, B = 50, validate = FALSE, seed = 603
  )
  r2 <- unconfoundedness_test(
    rct2, obs2, Y ~ A + X1 + X2, B = 50, validate = FALSE, seed = 604
  )

  results <- list(study1 = r1, study2 = r2)
  output <- capture.output(latex <- create_latex_table(results))
  expect_true(any(grepl("study1", latex)))
  expect_true(any(grepl("study2", latex)))
})

test_that("create_latex_table rejects invalid input", {
  expect_error(create_latex_table(42), "must be")
})

test_that("validate_installation returns logical", {
  output <- capture.output(result <- validate_installation())
  expect_true(is.logical(result))
})

test_that("batch analysis handles errors in individual analyses gracefully", {
  set.seed(404)

  rct_list <- list(
    good = generate_rct_data(n = 100, seed = 404),
    bad = data.frame(Y = 1:5, A = c(1, 1, 1, 1, 1), X1 = 1:5, X2 = 1:5)
  )
  obs_list <- list(
    good = generate_obs_data(n = 150, seed = 504),
    bad = generate_obs_data(n = 150, seed = 505)
  )

  # Should not throw; the bad analysis should be captured internally
  expect_output(
    batch_results <- batch_analysis(
      rct_list, obs_list, Y ~ A + X1 + X2,
      B = 50, validate = FALSE
    ),
    "Running"
  )

  expect_s3_class(batch_results, "batch_unconf_test")
  expect_equal(length(batch_results), 2)
})

test_that("batch analysis rejects mismatched list lengths", {
  rct_list <- list(generate_rct_data(n = 100, seed = 1))
  obs_list <- list(
    generate_obs_data(n = 100, seed = 2),
    generate_obs_data(n = 100, seed = 3)
  )

  expect_error(
    batch_analysis(rct_list, obs_list, Y ~ A + X1 + X2),
    "same length"
  )
})
