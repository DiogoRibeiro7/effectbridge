# unconfoundedr

> Test (un)confoundedness by comparing an effect from an RCT-like dataset to the same estimand from an observational dataset. Supports IPW/AIPW, bootstrap CIs, a Wald test, and optional transportability weighting (manual or auto-detected via KS/energy tests).

[![R-CMD-check](https://img.shields.io/github/actions/workflow/status/diogoribeiro7/unconfoundedr/R-CMD-check.yaml?label=R-CMD-check)](https://github.com/diogoribeiro7/unconfoundedr/actions/workflows/R-CMD-check.yaml) [![pkgdown](https://img.shields.io/github/actions/workflow/status/diogoribeiro7/unconfoundedr/pkgdown.yaml?label=pkgdown)](https://github.com/diogoribeiro7/unconfoundedr/actions/workflows/pkgdown.yaml) [![Codecov](https://codecov.io/gh/diogoribeiro7/unconfoundedr/branch/main/graph/badge.svg)](https://codecov.io/gh/diogoribeiro7/unconfoundedr)

## Why this package

Before running a full causal workflow on observational data, it's useful to check whether the **same estimand** agrees between a randomized trial (or trial-like slice) and your observational cohort. If the effects differ beyond sampling error, you have evidence against ignorability for the observational analysis (or a transportability failure). This package implements that comparison with practical estimators and robust inference.

## Features

- **Estimators**: `IPW` and **`AIPW` (doubly robust)** for the marginal ATE (risk difference for binary outcomes).
- **Inference**: percentile **bootstrap** CI for the difference Δ = ω_obs − ω_rct; **Wald test** (AIPW uses IF variance when no transport weighting is applied).
- **Transportability weighting**:

  - `transport = "none"`: no reweighting (default).
  - `transport = "rct_to_obs"`: reweight RCT to the OBS covariate mix via a density-ratio logistic model.
  - `transport = "auto"`: run shift tests (univariate **KS** with BH correction and/or multivariate **energy** test) and reweight only if shift is detected.

- **Stability tools**: stabilized weights, quantile trimming, positivity diagnostics, and transport weights' effective sample size (ESS).

## Installation

```r
# install.packages("devtools")
devtools::install_github("yourname/unconfoundedr")

# Optional for multivariate shift testing in transport = "auto"
install.packages("energy")
```

## Quick start

```r
library(unconfoundedr)
set.seed(42)

# Simulate an RCT-like dataset
n_r <- 800
X1 <- rnorm(n_r); X2 <- rbinom(n_r, 1, 0.4)
A_r <- rbinom(n_r, 1, 0.5)    # randomized
Y0_r <- 0.5 + 0.5*X1 + 0.3*X2 + rnorm(n_r)
Y1_r <- Y0_r + 1.0
rct_df <- data.frame(Y = ifelse(A_r==1, Y1_r, Y0_r), A = A_r, X1, X2)

# Simulate an observational dataset (with covariate shift)
n_o <- 2000
X1 <- rnorm(n_o, 0.4); X2 <- rbinom(n_o, 1, 0.7)
A_o <- rbinom(n_o, 1, plogis(-0.2 + 0.8*X1 + 0.6*X2))
U   <- rnorm(n_o)
Y0_o <- 0.5 + 0.5*X1 + 0.3*X2 + 0.3*U + rnorm(n_o)
Y1_o <- Y0_o + 1.0
obs_df <- data.frame(Y = ifelse(A_o==1, Y1_o, Y0_o), A = A_o, X1, X2)

# Compare effects with auto transport detection
out <- unconfoundedness_test(
  data_rct  = rct_df,
  data_obs  = obs_df,
  formula   = Y ~ A + X1 + X2,  # A must be the first RHS term
  estimator = "aipw",
  family_y  = "gaussian",
  transport = "auto",           # "none" | "rct_to_obs" | "auto"
  auto_method = "both",         # "ks" | "energy" | "both"
  auto_alpha  = 0.01,
  B = 500,
  seed = 42
)

print(out)
out$summary
out$diagnostics$auto       # KS/energy decision details
out$diagnostics$transport  # weight ranges and ESS when applied
```

## API(high-level)

```r
unconfoundedness_test(
  data_rct, data_obs, formula,
  estimator = c("aipw","ipw"),
  stabilize = TRUE,
  trim = c(0.01, 0.99),
  family_y = c("gaussian","binomial"),
  transport = c("none","rct_to_obs","auto"),
  auto_method = c("both","ks","energy"),
  auto_alpha = 0.01,
  auto_energy_R = 199L,
  B = 1000L,
  alpha = 0.05,
  seed = NULL
)
```

## Interpretation & tips

- **Null.** H0: ω_obs = ω_rct (risk-difference scale). Failure to reject supports ignorability under comparability assumptions.
- **Rejection.** Could be unmeasured confounding, bad controls, or transportability failure. Inspect diagnostics.
- **Positivity.** Large extreme-PS shares or tiny transport ESS → increase trimming or adjust covariates.
- **Speed.** Lower `B` to prototype; use `B >= 500` for final.
- **Reproducibility.** Set `seed`.

## License

MIT. See `LICENSE`.

## Citation

Please cite **unconfoundedr** and the underlying paper you benchmark against. See `CITATION`.
