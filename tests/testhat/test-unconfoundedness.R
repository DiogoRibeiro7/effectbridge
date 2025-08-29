test_that("auto transport stays off with no covariate shift", {
  set.seed(100)
  gen <- function(n){
    X1 <- rnorm(n); X2 <- rbinom(n,1,0.4)
    A  <- rbinom(n,1,plogis(-0.2 + 0.8*X1 + 0.6*X2))
    Y0 <- 0.5 + 0.5*X1 + 0.3*X2 + rnorm(n)
    Y1 <- Y0 + 1.0
    data.frame(Y=ifelse(A==1,Y1,Y0),A,X1,X2)
  }
  d_rct <- gen(600); d_obs <- gen(1500)
  out <- unconfoundedness_test(d_rct, d_obs, Y ~ A + X1 + X2,
                               estimator="aipw", family_y="gaussian",
                               transport="auto", auto_method="both",
                               auto_alpha=0.01, B=200, seed=100)
  expect_identical(out$transport_applied, "none")
  expect_false(out$diagnostics$auto$detected)
})

test_that("auto transport turns on with covariate shift", {
  set.seed(101)
  gen_rct <- function(n){
    X1 <- rnorm(n); X2 <- rbinom(n,1,0.3)
    A  <- rbinom(n,1,0.5)
    Y0 <- 0.5 + 0.5*X1 + 0.3*X2 + rnorm(n)
    Y1 <- Y0 + 1.0
    data.frame(Y=ifelse(A==1,Y1,Y0),A,X1,X2)
  }
  gen_obs <- function(n){
    X1 <- rnorm(n, 0.7); X2 <- rbinom(n,1,0.7)  # shift
    A  <- rbinom(n,1,plogis(-0.2 + 0.8*X1 + 0.6*X2))
    Y0 <- 0.5 + 0.5*X1 + 0.3*X2 + rnorm(n)
    Y1 <- Y0 + 1.0
    data.frame(Y=ifelse(A==1,Y1,Y0),A,X1,X2)
  }
  d_rct <- gen_rct(600); d_obs <- gen_obs(1500)
  out <- unconfoundedness_test(d_rct, d_obs, Y ~ A + X1 + X2,
                               estimator="aipw", family_y="gaussian",
                               transport="auto", auto_method="ks",
                               auto_alpha=0.01, B=200, seed=101)
  expect_identical(out$transport_applied, "rct_to_obs")
  expect_true(out$diagnostics$auto$detected)
})
