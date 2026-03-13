# estimators.R - Estimator implementations and related helpers

#' Compute effects using specified estimator
#' @param parsed_data Parsed data object
#' @param estimator Estimation method
#' @param effect_measure Effect measure
#' @param family_y Outcome family
#' @param w_rct RCT weights
#' @param w_obs Observational weights
#' @param stabilize Stabilize weights
#' @param trim Trimming bounds
#' @return Effects object
compute_effects <- function(
  parsed_data,
  estimator,
  effect_measure,
  family_y,
  w_rct,
  w_obs,
  stabilize,
  trim
) {
  # Compute RCT effect
  rct_effect <- estimate_single_effect(
    Y = parsed_data$Y_rct,
    A = parsed_data$A_rct,
    X = parsed_data$X_rct,
    weights = w_rct,
    estimator = estimator,
    effect_measure = effect_measure,
    family_y = family_y,
    stabilize = stabilize,
    trim = trim
  )

  # Compute observational effect
  obs_effect <- estimate_single_effect(
    Y = parsed_data$Y_obs,
    A = parsed_data$A_obs,
    X = parsed_data$X_obs,
    weights = w_obs,
    estimator = estimator,
    effect_measure = effect_measure,
    family_y = family_y,
    stabilize = stabilize,
    trim = trim
  )

  list(
    rct_effect = rct_effect,
    obs_effect = obs_effect
  )
}

#' Estimate effect for single dataset
#' @param Y Outcome vector
#' @param A Treatment vector
#' @param X Covariate matrix
#' @param weights Sample weights
#' @param estimator Estimation method
#' @param effect_measure Effect measure
#' @param family_y Outcome family
#' @param stabilize Stabilize weights
#' @param trim Trimming bounds
#' @return Single effect estimate
estimate_single_effect <- function(
  Y,
  A,
  X,
  weights,
  estimator,
  effect_measure,
  family_y,
  stabilize,
  trim
) {
  switch(
    estimator,
    "ipw" = estimate_ipw(
      Y,
      A,
      X,
      weights,
      effect_measure,
      family_y,
      stabilize,
      trim
    ),
    "aipw" = estimate_aipw(
      Y,
      A,
      X,
      weights,
      effect_measure,
      family_y,
      stabilize,
      trim
    ),
    "tmle" = estimate_tmle(Y, A, X, weights, effect_measure, family_y),
    "gcomp" = estimate_gcomp(Y, A, X, weights, effect_measure, family_y),
    "matching" = estimate_matching(Y, A, X, weights, effect_measure, family_y),
    stop("Unknown estimator: ", estimator)
  )
}

#' IPW Estimator
#' @param Y Outcome vector
#' @param A Treatment vector
#' @param X Covariate matrix
#' @param weights Sample weights
#' @param effect_measure Effect measure
#' @param family_y Outcome family
#' @param stabilize Stabilize weights
#' @param trim Trimming bounds
#' @return IPW estimate
estimate_ipw <- function(
  Y,
  A,
  X,
  weights,
  effect_measure,
  family_y,
  stabilize,
  trim
) {
  # Fit propensity score model
  ps_fit <- fit_propensity_model(A, X, weights)
  e <- ps_fit$fitted_values

  # Compute IPW weights
  ipw_weights <- compute_ipw_weights(A, e, weights, stabilize)

  # Apply trimming if specified
  if (!is.null(trim)) {
    keep_idx <- apply_weight_trimming(ipw_weights, trim)
    Y <- Y[keep_idx]
    A <- A[keep_idx]
    ipw_weights <- ipw_weights[keep_idx]
    e <- e[keep_idx]
  }

  # Estimate potential outcomes
  mu1 <- weighted.mean(Y[A == 1], ipw_weights[A == 1])
  mu0 <- weighted.mean(Y[A == 0], ipw_weights[A == 0])

  # Compute effect measure
  effect <- compute_effect_measure(mu1, mu0, effect_measure)

  # Compute influence function for variance
  IF <- compute_ipw_influence_function(
    Y,
    A,
    e,
    ipw_weights,
    mu1,
    mu0,
    effect_measure
  )

  list(
    estimate = effect,
    mu1 = mu1,
    mu0 = mu0,
    variance = var(IF, na.rm = TRUE),
    influence_function = IF,
    propensity_scores = e,
    weights = ipw_weights,
    details = list(
      ps_model = ps_fit,
      method = "ipw"
    )
  )
}

#' AIPW (Doubly Robust) Estimator
#' @param Y Outcome vector
#' @param A Treatment vector
#' @param X Covariate matrix
#' @param weights Sample weights
#' @param effect_measure Effect measure
#' @param family_y Outcome family
#' @param stabilize Stabilize weights
#' @param trim Trimming bounds
#' @return AIPW estimate
estimate_aipw <- function(
  Y,
  A,
  X,
  weights,
  effect_measure,
  family_y,
  stabilize,
  trim
) {
  # Fit propensity score model
  ps_fit <- fit_propensity_model(A, X, weights)
  e <- ps_fit$fitted_values

  # Fit outcome regression models
  outcome_fits <- fit_outcome_models(Y, A, X, weights, family_y)
  mu1_hat <- outcome_fits$mu1_pred
  mu0_hat <- outcome_fits$mu0_pred

  # Compute AIPW pseudo-outcomes
  pseudo_outcomes <- compute_aipw_pseudo_outcomes(
    Y,
    A,
    e,
    mu1_hat,
    mu0_hat,
    weights,
    stabilize
  )

  # Apply trimming if specified
  if (!is.null(trim)) {
    base_weights <- weights * (A / e + (1 - A) / (1 - e))
    keep_idx <- apply_weight_trimming(base_weights, trim)
    Y <- Y[keep_idx]
    A <- A[keep_idx]
    e <- e[keep_idx]
    mu1_hat <- mu1_hat[keep_idx]
    mu0_hat <- mu0_hat[keep_idx]
    weights <- weights[keep_idx]
    pseudo_outcomes <- compute_aipw_pseudo_outcomes(
      Y,
      A,
      e,
      mu1_hat,
      mu0_hat,
      weights,
      stabilize
    )
  }

  # Estimate potential outcomes
  mu1 <- weighted.mean(pseudo_outcomes$mu1, weights)
  mu0 <- weighted.mean(pseudo_outcomes$mu0, weights)

  # Compute effect measure
  effect <- compute_effect_measure(mu1, mu0, effect_measure)

  # Compute influence function
  IF <- compute_aipw_influence_function(
    pseudo_outcomes,
    weights,
    mu1,
    mu0,
    effect_measure
  )

  list(
    estimate = effect,
    mu1 = mu1,
    mu0 = mu0,
    variance = var(IF, na.rm = TRUE),
    influence_function = IF,
    propensity_scores = e,
    outcome_predictions = list(mu1 = mu1_hat, mu0 = mu0_hat),
    details = list(
      ps_model = ps_fit,
      outcome_models = outcome_fits,
      method = "aipw"
    )
  )
}

#' TMLE Estimator
#' @param Y Outcome vector
#' @param A Treatment vector
#' @param X Covariate matrix
#' @param weights Sample weights
#' @param effect_measure Effect measure
#' @param family_y Outcome family
#' @return TMLE estimate
estimate_tmle <- function(Y, A, X, weights, effect_measure, family_y) {
  if (!requireNamespace("tmle", quietly = TRUE)) {
    stop("Package 'tmle' required for TMLE estimator")
  }

  # Convert matrix to data frame for tmle
  if (ncol(X) > 1 || !all(X[, 1] == 1)) {
    W <- X[, !colnames(X) %in% "(Intercept)", drop = FALSE]
    if (ncol(W) == 0) W <- NULL
  } else {
    W <- NULL
  }

  # Fit TMLE
  tmle_fit <- tryCatch(
    {
      tmle::tmle(
        Y = Y,
        A = A,
        W = W,
        family = if (family_y == "binomial") "binomial" else "gaussian",
        V = 5 # 5-fold cross-validation
      )
    },
    error = function(e) {
      # Fallback to simpler TMLE if cross-validation fails
      tmle::tmle(
        Y = Y,
        A = A,
        W = W,
        family = if (family_y == "binomial") "binomial" else "gaussian"
      )
    }
  )

  # Extract estimates
  mu1 <- tmle_fit$estimates$EY1$psi
  mu0 <- tmle_fit$estimates$EY0$psi
  effect_rd <- tmle_fit$estimates$ATE$psi

  # Convert to requested effect measure
  if (effect_measure == "rd") {
    effect <- effect_rd
  } else if (effect_measure == "rr") {
    effect <- if (mu0 != 0) mu1 / mu0 else NA_real_
  } else if (effect_measure == "or") {
    if (mu0 != 0 && mu0 != 1 && mu1 != 0 && mu1 != 1) {
      effect <- (mu1 / (1 - mu1)) / (mu0 / (1 - mu0))
    } else {
      effect <- NA_real_
    }
  }

  # Variance (for risk difference)
  variance <- tmle_fit$estimates$ATE$var.psi

  # Adjust variance for other measures using delta method
  if (effect_measure == "rr" && mu0 != 0 && !is.na(effect)) {
    # Delta method for log risk ratio
    grad <- c(1 / mu0, -mu1 / (mu0^2))
    var_mu <- matrix(
      c(tmle_fit$estimates$EY1$var.psi, 0, 0, tmle_fit$estimates$EY0$var.psi),
      2,
      2
    )
    variance <- as.numeric(t(grad) %*% var_mu %*% grad * (effect^2))
  } else if (effect_measure == "or" && !is.na(effect)) {
    # Delta method for log odds ratio (approximate)
    if (mu1 != 0 && mu1 != 1 && mu0 != 0 && mu0 != 1) {
      variance <- tmle_fit$estimates$ATE$var.psi *
        (effect^2) /
        ((mu1 * (1 - mu1) + mu0 * (1 - mu0))^2)
    }
  }

  list(
    estimate = effect,
    mu1 = mu1,
    mu0 = mu0,
    variance = variance,
    influence_function = tmle_fit$estimates$ATE$IC,
    details = list(
      tmle_fit = tmle_fit,
      method = "tmle"
    )
  )
}

#' G-computation Estimator
#' @param Y Outcome vector
#' @param A Treatment vector
#' @param X Covariate matrix
#' @param weights Sample weights
#' @param effect_measure Effect measure
#' @param family_y Outcome family
#' @return G-computation estimate
estimate_gcomp <- function(Y, A, X, weights, effect_measure, family_y) {
  # Fit outcome models
  outcome_fits <- fit_outcome_models(Y, A, X, weights, family_y)

  # Predict potential outcomes for all subjects
  mu1_hat <- outcome_fits$mu1_pred
  mu0_hat <- outcome_fits$mu0_pred

  # G-computation estimates
  mu1 <- weighted.mean(mu1_hat, weights)
  mu0 <- weighted.mean(mu0_hat, weights)

  # Compute effect measure
  effect <- compute_effect_measure(mu1, mu0, effect_measure)

  # Bootstrap-based variance estimation
  n <- length(Y)
  W_scaled <- weights / sum(weights)

  # Influence function approximation for G-computation
  IF <- W_scaled * (mu1_hat - mu0_hat - effect)
  variance <- var(IF, na.rm = TRUE)

  list(
    estimate = effect,
    mu1 = mu1,
    mu0 = mu0,
    variance = variance,
    influence_function = IF,
    outcome_predictions = list(mu1 = mu1_hat, mu0 = mu0_hat),
    details = list(
      outcome_models = outcome_fits,
      method = "gcomp"
    )
  )
}

#' Matching Estimator
#' @param Y Outcome vector
#' @param A Treatment vector
#' @param X Covariate matrix
#' @param weights Sample weights
#' @param effect_measure Effect measure
#' @param family_y Outcome family
#' @return Matching estimate
estimate_matching <- function(Y, A, X, weights, effect_measure, family_y) {
  if (!requireNamespace("MatchIt", quietly = TRUE)) {
    stop("Package 'MatchIt' required for matching estimator")
  }

  # Prepare data
  W <- X[, !colnames(X) %in% "(Intercept)", drop = FALSE]
  if (ncol(W) == 0) {
    stop("Matching requires covariates")
  }

  data_match <- data.frame(Y = Y, A = A, W, weights = weights)

  # Perform matching
  match_fit <- tryCatch(
    {
      MatchIt::matchit(
        A ~ .,
        data = data_match[, !names(data_match) %in% c("Y", "weights")],
        method = "nearest",
        distance = "glm",
        replace = FALSE,
        ratio = 1
      )
    },
    error = function(e) {
      # Fallback to exact matching if available
      if (ncol(W) <= 3) {
        MatchIt::matchit(
          A ~ .,
          data = data_match[, !names(data_match) %in% c("Y", "weights")],
          method = "exact"
        )
      } else {
        stop("Matching failed: ", e$message)
      }
    }
  )

  # Extract matched data
  matched_data <- MatchIt::match.data(match_fit)

  # Compute effect on matched sample
  Y_matched <- matched_data$Y
  A_matched <- matched_data$A
  match_weights <- matched_data$weights

  # Simple difference in means on matched sample
  mu1 <- weighted.mean(Y_matched[A_matched == 1], match_weights[A_matched == 1])
  mu0 <- weighted.mean(Y_matched[A_matched == 0], match_weights[A_matched == 0])

  # Compute effect measure
  effect <- compute_effect_measure(mu1, mu0, effect_measure)

  # Approximate variance using matched sample
  n1 <- sum(A_matched == 1)
  n0 <- sum(A_matched == 0)

  if (n1 > 1 && n0 > 1) {
    s1_sq <- wtd.var(Y_matched[A_matched == 1], match_weights[A_matched == 1])
    s0_sq <- wtd.var(Y_matched[A_matched == 0], match_weights[A_matched == 0])

    if (effect_measure == "rd") {
      variance <- s1_sq / n1 + s0_sq / n0
    } else {
      # Delta method approximation for other measures
      if (effect_measure == "rr" && mu0 != 0) {
        variance <- (s1_sq / n1) / (mu0^2) + (s0_sq / n0) * (mu1^2) / (mu0^4)
      } else if (
        effect_measure == "or" && mu0 != 0 && mu0 != 1 && mu1 != 0 && mu1 != 1
      ) {
        variance <- (s1_sq / n1) /
          (mu1^2 * (1 - mu1)^2) +
          (s0_sq / n0) / (mu0^2 * (1 - mu0)^2)
      } else {
        variance <- NA_real_
      }
    }
  } else {
    variance <- NA_real_
  }

  list(
    estimate = effect,
    mu1 = mu1,
    mu0 = mu0,
    variance = variance,
    influence_function = NULL, # Not easily available for matching
    details = list(
      match_fit = match_fit,
      matched_data = matched_data,
      n_matched = nrow(matched_data),
      method = "matching"
    )
  )
}

# Helper functions for estimators

#' Extract covariate names from a model formula
#'
#' Given a formula like `Y ~ A + X1 + X2`, return the covariate names
#' to be used in the propensity model: all symbols on the RHS except `treat`.
#' Works with interactions (X1:X2), functions (log(X1)), and backticked names.
#' If a variable isn't a column in `data`, it's dropped.
#'
#' @param formula stats::formula like `Y ~ A + X1 + X2`
#' @param data data.frame used to validate which names actually exist
#' @param treat character(1) name of the treatment variable (default "A")
#' @return character() vector of covariate names present in `data`
#' @examples
#' df <- data.frame(Y=0, A=0, X1=1, X2=2, X3=3)
#' extract_covariates(Y ~ A + X1 + X2, df)          # c("X1","X2")
#' extract_covariates(Y ~ A + log(X1) + X2:X3, df)  # c("X1","X2","X3")
#'
#' @keywords internal
extract_covariates <- function(formula, data, treat = "A") {
  # --- type checks ------------------------------------------------------------
  if (!inherits(formula, "formula")) {
    stop("`formula` must be a stats::formula.", call. = FALSE)
  }
  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame.", call. = FALSE)
  }
  if (!is.character(treat) || length(treat) != 1L || is.na(treat) || treat == "") {
    stop("`treat` must be a non-empty character(1).", call. = FALSE)
  }

  # --- pull symbols from RHS only --------------------------------------------
  rhs <- formula[[3L]]
  rhs_vars <- all.vars(rhs)

  # --- drop treatment and keep only columns that exist ------------------------
  covars <- setdiff(rhs_vars, treat)
  covars <- intersect(covars, colnames(data))

  covars
}

#' Build a propensity-score formula safely
#'
#' Constructs a formula of the form \code{A ~ x1 + x2 + ...} from a set of
#' predictor names. If \code{predictors} is empty, it falls back to the
#' intercept-only model \code{A ~ 1} and emits a warning.
#'
#' @param predictors character()
#'   Vector of predictor (covariate) column names. May be length 0.
#'
#' @return formula
#'   A valid formula object for use in \code{glm()}.
#'
#' @examples
#' build_ps_formula(c("X1","X2"))  # A ~ X1 + X2
#' build_ps_formula(character(0))  # A ~ 1 (with warning)
#'
#' @keywords internal
build_ps_formula <- function(predictors) {
  # ---- Type checks -----------------------------------------------------------
  if (!is.character(predictors)) {
    stop("`predictors` must be a character vector of column names.", call. = FALSE)
  }
  if (anyNA(predictors)) {
    stop("`predictors` contains NA values.", call. = FALSE)
  }

  # ---- Empty set => intercept-only ------------------------------------------
  if (length(predictors) == 0L) {
    warning(
      "No covariates available for the propensity model; using intercept-only model (A ~ 1).",
      call. = FALSE
    )
    return(A ~ 1)
  }

  # ---- Construct A ~ x1 + x2 + ... ------------------------------------------
  stats::as.formula(paste("A ~", paste(predictors, collapse = " + ")))
}

#' Fit propensity score model
#'
#' Builds a logistic regression for \code{A} given covariates \code{X}.
#' Handles the cases where \code{X} has only an intercept or no columns
#' by falling back to \code{A ~ 1}. Also stabilizes fitted probabilities to
#' avoid exact 0/1 weights downstream.
#'
#' @param A numeric|integer|logical
#'   Treatment indicator (0/1). Coerced to integer 0/1 if logical.
#' @param X matrix|data.frame
#'   Covariate design (may include an \code{"(Intercept)"} column). Row count
#'   must match \code{length(A)}. Column names are required (added if missing).
#' @param weights numeric|NULL
#'   Optional sampling weights of length \code{length(A)}. Defaults to 1.
#'
#' @return list
#'   \itemize{
#'     \item \code{model}: the \code{glm} fit (or a placeholder list on fallback)
#'     \item \code{fitted_values}: numeric vector of stabilized propensities
#'     \item \code{formula}: the formula used
#'   }
#'
#' @examples
#' set.seed(1)
#' A <- rbinom(100, 1, 0.4)
#' X <- cbind(`(Intercept)` = 1, X1 = rnorm(100), X2 = rnorm(100))
#' fit <- fit_propensity_model(A, X, weights = NULL)
#' head(fit$fitted_values)
#'
#' @keywords internal
fit_propensity_model <- function(A, X, weights = NULL) {
  # ---- Type checks -----------------------------------------------------------
  n <- length(A)
  if (is.logical(A)) A <- as.integer(A)
  if (!is.numeric(A) || any(!A %in% c(0, 1))) {
    stop("`A` must be a binary vector (0/1 or logical).", call. = FALSE)
  }

  if (!is.matrix(X) && !is.data.frame(X)) {
    stop("`X` must be a matrix or data.frame.", call. = FALSE)
  }
  if (nrow(X) != n) {
    stop("`nrow(X)` must equal length(A).", call. = FALSE)
  }

  # Normalize colnames (required for formula terms)
  if (is.null(colnames(X))) {
    colnames(X) <- paste0("X", seq_len(ncol(X)))
  }

  # Weights default to 1
  if (is.null(weights)) {
    weights <- rep(1, n)
  } else {
    if (!is.numeric(weights) || length(weights) != n) {
      stop("`weights` must be numeric and the same length as `A`.", call. = FALSE)
    }
  }

  # ---- Build predictors data.frame (drop intercept-like columns) ------------
  X_df <- as.data.frame(X, check.names = TRUE)
  is_intercept_col <- function(col, nm) {
    (nm == "(Intercept)") || (is.numeric(col) && all(col == 1))
  }
  keep <- !mapply(is_intercept_col, X_df, names(X_df))
  predictors_df <- X_df[, keep, drop = FALSE]

  # ---- Build safe formula ----------------------------------------------------
  predictor_names <- colnames(predictors_df)
  f_ps <- build_ps_formula(predictor_names)

  # ---- Assemble modeling data.frame -----------------------------------------
  df <- if (length(predictor_names)) {
    data.frame(A = A, predictors_df, check.names = TRUE)
  } else {
    data.frame(A = A, check.names = TRUE)
  }

  # ---- Fit model (robust fallback) ------------------------------------------
  fit <- tryCatch(
    stats::glm(f_ps, data = df, family = stats::binomial(), weights = weights),
    error = function(e) {
      # Separation / singularities / other failures: fallback to marginal p
      warning(
        "Propensity score model fitting failed; using marginal treated proportion as PS.",
        call. = FALSE
      )
      list(fitted.values = rep(mean(A), n))
    }
  )

  # ---- Stabilize fitted probabilities ---------------------------------------
  if (inherits(fit, "glm")) {
    p <- stats::fitted(fit)
  } else {
    p <- fit$fitted.values
  }
  # clamp to (1e-3, 1-1e-3) to avoid extreme weights
  p <- pmin(pmax(p, 1e-3), 1 - 1e-3)

  list(
    model = fit,
    fitted_values = p,
    formula = f_ps
  )
}

#' Fit outcome regression models
#' @param Y Outcome vector
#' @param A Treatment vector
#' @param X Covariate matrix
#' @param weights Sample weights
#' @param family_y Outcome family
#' @return Outcome model fits
fit_outcome_models <- function(Y, A, X, weights, family_y) {
  # Prepare predictors
  if (ncol(X) > 1 || !all(X[, 1] == 1)) {
    predictors <- X[, !colnames(X) %in% "(Intercept)", drop = FALSE]
    if (ncol(predictors) == 0) {
      use_predictors <- FALSE
    } else {
      use_predictors <- TRUE
    }
  } else {
    use_predictors <- FALSE
  }

  # Fit model for treated (A=1)
  idx_1 <- A == 1
  if (sum(idx_1) < 2) {
    stop("Insufficient treated observations for outcome modeling")
  }

  if (use_predictors) {
    df_1 <- data.frame(Y = Y[idx_1], predictors[idx_1, , drop = FALSE])
    formula_1 <- as.formula(paste(
      "Y ~",
      paste(colnames(predictors), collapse = " + ")
    ))
  } else {
    df_1 <- data.frame(Y = Y[idx_1])
    formula_1 <- Y ~ 1
  }

  family_obj <- if (family_y == "binomial") binomial() else gaussian()

  fit_1 <- tryCatch(
    {
      glm(formula_1, data = df_1, family = family_obj, weights = weights[idx_1])
    },
    error = function(e) {
      # Fallback to intercept-only model
      glm(
        Y ~ 1,
        data = data.frame(Y = Y[idx_1]),
        family = family_obj,
        weights = weights[idx_1]
      )
    }
  )

  # Fit model for control (A=0)
  idx_0 <- A == 0
  if (sum(idx_0) < 2) {
    stop("Insufficient control observations for outcome modeling")
  }

  if (use_predictors) {
    df_0 <- data.frame(Y = Y[idx_0], predictors[idx_0, , drop = FALSE])
    formula_0 <- as.formula(paste(
      "Y ~",
      paste(colnames(predictors), collapse = " + ")
    ))
  } else {
    df_0 <- data.frame(Y = Y[idx_0])
    formula_0 <- Y ~ 1
  }

  fit_0 <- tryCatch(
    {
      glm(formula_0, data = df_0, family = family_obj, weights = weights[idx_0])
    },
    error = function(e) {
      # Fallback to intercept-only model
      glm(
        Y ~ 1,
        data = data.frame(Y = Y[idx_0]),
        family = family_obj,
        weights = weights[idx_0]
      )
    }
  )

  # Predict for all observations
  if (use_predictors && inherits(fit_1, "glm") && inherits(fit_0, "glm")) {
    newdata <- data.frame(predictors)
    names(newdata) <- colnames(predictors)
  } else {
    newdata <- data.frame(row.names = seq_along(Y))
  }

  mu1_pred <- tryCatch(
    {
      predict(fit_1, newdata = newdata, type = "response")
    },
    error = function(e) {
      rep(mean(Y[idx_1]), length(Y))
    }
  )

  mu0_pred <- tryCatch(
    {
      predict(fit_0, newdata = newdata, type = "response")
    },
    error = function(e) {
      rep(mean(Y[idx_0]), length(Y))
    }
  )

  list(
    fit_1 = fit_1,
    fit_0 = fit_0,
    mu1_pred = as.numeric(mu1_pred),
    mu0_pred = as.numeric(mu0_pred),
    formulas = list(treated = formula_1, control = formula_0)
  )
}

#' Compute IPW weights
#' @param A Treatment vector
#' @param e Propensity scores
#' @param weights Sample weights
#' @param stabilize Whether to stabilize
#' @return IPW weights
compute_ipw_weights <- function(A, e, weights, stabilize) {
  if (stabilize) {
    # Stabilized weights
    p_a <- mean(A) # Marginal treatment probability
    numerator <- A * p_a + (1 - A) * (1 - p_a)
  } else {
    numerator <- 1
  }

  ipw_weights <- weights * numerator / (A * e + (1 - A) * (1 - e))

  # Check for extreme weights
  if (any(is.infinite(ipw_weights)) || any(ipw_weights < 0)) {
    warning("Extreme IPW weights detected, check propensity score model")
    ipw_weights[is.infinite(ipw_weights)] <- max(
      ipw_weights[is.finite(ipw_weights)],
      na.rm = TRUE
    )
    ipw_weights[ipw_weights < 0] <- 0
  }

  ipw_weights
}

#' Apply weight trimming
#' @param weights Weight vector
#' @param trim Trimming bounds
#' @return Indices to keep
apply_weight_trimming <- function(weights, trim) {
  if (is.null(trim)) {
    return(seq_along(weights))
  }

  lower_bound <- quantile(weights, trim[1], na.rm = TRUE)
  upper_bound <- quantile(weights, trim[2], na.rm = TRUE)

  which(weights >= lower_bound & weights <= upper_bound & !is.na(weights))
}

#' Compute effect measure from potential outcome means
#' @param mu1 Mean under treatment
#' @param mu0 Mean under control
#' @param effect_measure Effect measure
#' @return Effect estimate
compute_effect_measure <- function(mu1, mu0, effect_measure) {
  switch(
    effect_measure,
    "rd" = mu1 - mu0, # Risk difference
    "rr" = {
      # Risk ratio
      if (mu0 == 0) {
        warning("Control mean is 0, risk ratio undefined")
        NA_real_
      } else {
        mu1 / mu0
      }
    },
    "or" = {
      # Odds ratio
      if (mu0 == 0 || mu0 == 1 || mu1 == 0 || mu1 == 1) {
        warning("Probabilities at boundary, odds ratio undefined")
        NA_real_
      } else {
        (mu1 / (1 - mu1)) / (mu0 / (1 - mu0))
      }
    },
    stop("Unknown effect measure: ", effect_measure)
  )
}

#' Compute AIPW pseudo-outcomes
#' @param Y Outcome vector
#' @param A Treatment vector
#' @param e Propensity scores
#' @param mu1_hat Outcome predictions under treatment
#' @param mu0_hat Outcome predictions under control
#' @param weights Sample weights
#' @param stabilize Whether to stabilize
#' @return Pseudo-outcomes
compute_aipw_pseudo_outcomes <- function(
  Y,
  A,
  e,
  mu1_hat,
  mu0_hat,
  weights,
  stabilize
) {
  if (stabilize) {
    p_a <- mean(A)
    numerator_1 <- p_a
    numerator_0 <- 1 - p_a
  } else {
    numerator_1 <- numerator_0 <- 1
  }

  # AIPW pseudo-outcomes
  mu1_pseudo <- mu1_hat + numerator_1 * A * (Y - mu1_hat) / e
  mu0_pseudo <- mu0_hat + numerator_0 * (1 - A) * (Y - mu0_hat) / (1 - e)

  list(mu1 = mu1_pseudo, mu0 = mu0_pseudo)
}

#' Compute IPW influence function
#' @param Y Outcome vector
#' @param A Treatment vector
#' @param e Propensity scores
#' @param weights IPW weights
#' @param mu1 Mean under treatment
#' @param mu0 Mean under control
#' @param effect_measure Effect measure
#' @return Influence function
compute_ipw_influence_function <- function(
  Y,
  A,
  e,
  weights,
  mu1,
  mu0,
  effect_measure
) {
  n <- length(Y)

  # IPW components
  w1 <- A / e
  w0 <- (1 - A) / (1 - e)

  # Normalize by sum of weights
  sum_w1 <- sum(w1)
  sum_w0 <- sum(w0)

  if1 <- w1 * (Y - mu1) / sum_w1
  if0 <- w0 * (Y - mu0) / sum_w0

  # Influence function for effect measure
  switch(
    effect_measure,
    "rd" = if1 - if0,
    "rr" = {
      if (mu0 != 0) {
        if1 / mu0 - mu1 * if0 / (mu0^2)
      } else {
        rep(NA_real_, n)
      }
    },
    "or" = {
      # Delta method for log odds ratio
      if (mu1 != 0 && mu1 != 1 && mu0 != 0 && mu0 != 1) {
        d_mu1 <- 1 / (mu1 * (1 - mu1))
        d_mu0 <- -1 / (mu0 * (1 - mu0))
        or_val <- (mu1 / (1 - mu1)) / (mu0 / (1 - mu0))
        or_val * (d_mu1 * if1 + d_mu0 * if0)
      } else {
        rep(NA_real_, n)
      }
    },
    stop("Unknown effect measure: ", effect_measure)
  )
}

#' Compute AIPW influence function
#' @param pseudo_outcomes Pseudo-outcomes
#' @param weights Sample weights
#' @param mu1 Mean under treatment
#' @param mu0 Mean under control
#' @param effect_measure Effect measure
#' @return Influence function
compute_aipw_influence_function <- function(
  pseudo_outcomes,
  weights,
  mu1,
  mu0,
  effect_measure
) {
  # Weighted mean influence functions
  W <- weights / sum(weights)

  if1 <- W * (pseudo_outcomes$mu1 - mu1)
  if0 <- W * (pseudo_outcomes$mu0 - mu0)

  # Influence function for effect measure
  switch(
    effect_measure,
    "rd" = if1 - if0,
    "rr" = {
      if (mu0 != 0) {
        if1 / mu0 - mu1 * if0 / (mu0^2)
      } else {
        rep(NA_real_, length(if1))
      }
    },
    "or" = {
      # Delta method for log odds ratio
      if (mu1 != 0 && mu1 != 1 && mu0 != 0 && mu0 != 1) {
        d_mu1 <- 1 / (mu1 * (1 - mu1))
        d_mu0 <- -1 / (mu0 * (1 - mu0))
        or_val <- (mu1 / (1 - mu1)) / (mu0 / (1 - mu0))
        or_val * (d_mu1 * if1 + d_mu0 * if0)
      } else {
        rep(NA_real_, length(if1))
      }
    },
    stop("Unknown effect measure: ", effect_measure)
  )
}

#' Weighted variance
#' @param x Vector
#' @param weights Weights
#' @return Weighted variance
wtd.var <- function(x, weights = NULL) {
  if (is.null(weights)) {
    return(var(x, na.rm = TRUE))
  }

  # Remove missing values
  valid_idx <- !is.na(x) & !is.na(weights)
  x <- x[valid_idx]
  weights <- weights[valid_idx]

  if (length(x) < 2) {
    return(NA_real_)
  }

  if (all(weights == weights[1])) {
    return(var(x))
  }

  w_mean <- weighted.mean(x, weights)
  sum(weights * (x - w_mean)^2) / (sum(weights) - 1)
}
