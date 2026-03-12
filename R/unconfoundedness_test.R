# unconfoundedness_test.R - Main function, input validation, and data parsing
# R >= 4.1

#' @title Test unconfoundedness with comprehensive diagnostics and multiple estimators
#' @description
#' Compare a marginal treatment effect in an RCT-like dataset with the same
#' estimand in an observational dataset. Enhanced version supports:
#' - Estimators: IPW, AIPW, TMLE, G-computation, matching
#' - Effect measures: risk difference, risk ratio, odds ratio
#' - Inference: bootstrap CI, analytical SE, robust SE
#' - Transport weighting with auto-detection
#' - Comprehensive diagnostics and visualization
#' - Sensitivity analysis tools
#'
#' @param data_rct,data_obs data.frame with same Y/A definitions and covariates
#' @param formula model formula `Y ~ A + X1 + X2`. A must be binary and the first RHS term
#' @param estimator character: "aipw" (default), "ipw", "tmle", "gcomp", "matching"
#' @param effect_measure character: "rd" (risk difference), "rr" (risk ratio), "or" (odds ratio)
#' @param stabilize logical; stabilized IPW (default TRUE)
#' @param trim length-2 numeric in \[0,1\] for weight trimming; NULL disables
#' @param family_y "gaussian" or "binomial"
#' @param transport "none", "rct_to_obs", or "auto" (default "none")
#' @param auto_method "both", "ks", "energy" for shift detection
#' @param auto_alpha significance level for shift tests (default 0.01)
#' @param auto_energy_R permutations for energy test (default 199)
#' @param inference_method "bootstrap", "analytical", "robust"
#' @param B bootstrap replicates (default 1000)
#' @param alpha CI level (default 0.05 for 95% CI)
#' @param parallel logical; use parallel processing for bootstrap
#' @param n_cores integer; cores for parallel processing (NULL = auto-detect)
#' @param seed RNG seed or NULL
#' @param validate logical; perform extensive input validation
#' @return enhanced unconf_test object with comprehensive diagnostics
#' @examples
#' \dontrun{
#' # Enhanced example with multiple estimators
#' set.seed(42)
#' data_rct <- generate_rct_data(n = 800)
#' data_obs <- generate_obs_data(n = 2000, confounding = TRUE)
#'
#' # Comprehensive test
#' result <- unconfoundedness_test(
#'   data_rct, data_obs,
#'   Y ~ A + X1 + X2,
#'   estimator = "aipw",
#'   effect_measure = "rd",
#'   transport = "auto",
#'   inference_method = "bootstrap",
#'   B = 1000,
#'   parallel = TRUE,
#'   validate = TRUE
#' )
#'
#' print(result)
#' plot(result)
#' summary(result)
#' }
#' @export
unconfoundedness_test <- function(
  data_rct,
  data_obs,
  formula,
  estimator = c("aipw", "ipw", "tmle", "gcomp", "matching"),
  effect_measure = c("rd", "rr", "or"),
  stabilize = TRUE,
  trim = c(0.01, 0.99),
  family_y = c("gaussian", "binomial"),
  transport = c("none", "rct_to_obs", "auto"),
  auto_method = c("both", "ks", "energy"),
  auto_alpha = 0.01,
  auto_energy_R = 199L,
  inference_method = c("bootstrap", "analytical", "robust"),
  B = 1000L,
  alpha = 0.05,
  parallel = FALSE,
  n_cores = NULL,
  seed = NULL,
  validate = TRUE
) {
  # Match arguments
  estimator <- match.arg(estimator)
  effect_measure <- match.arg(effect_measure)
  family_y <- match.arg(family_y)
  transport <- match.arg(transport)
  auto_method <- match.arg(auto_method)
  inference_method <- match.arg(inference_method)

  # Set seed early if provided
  if (!is.null(seed)) {
    set.seed(as.integer(seed))
  }

  # Enhanced input validation
  if (validate) {
    validation_result <- validate_inputs(
      data_rct,
      data_obs,
      formula,
      estimator,
      effect_measure,
      family_y,
      alpha,
      trim
    )
    if (!validation_result$valid) {
      stop("Input validation failed: ", validation_result$message)
    }
    # Display warnings
    if (length(validation_result$warnings) > 0) {
      for (w in validation_result$warnings) {
        warning(w, call. = FALSE)
      }
    }
  }

  # Parse formula and prepare data
  parsed_data <- parse_formula_and_data(formula, data_rct, data_obs)

  # Auto shift detection and transport decision
  transport_result <- decide_transport(
    parsed_data$X_rct,
    parsed_data$X_obs,
    transport,
    auto_method,
    auto_alpha,
    auto_energy_R
  )

  # Apply transport weights
  weights <- compute_transport_weights(
    parsed_data$X_rct,
    parsed_data$X_obs,
    transport_result$apply_transport
  )

  # Compute effects with chosen estimator
  effects <- compute_effects(
    parsed_data,
    estimator,
    effect_measure,
    family_y,
    weights$w_rct,
    weights$w_obs,
    stabilize,
    trim
  )

  # Compute difference
  diff_estimate <- compute_effect_difference(
    effects$rct_effect,
    effects$obs_effect,
    effect_measure
  )

  # Inference
  inference_result <- compute_inference(
    parsed_data,
    estimator,
    effect_measure,
    family_y,
    weights,
    effects,
    diff_estimate,
    inference_method,
    B,
    alpha,
    parallel,
    n_cores
  )

  # Comprehensive diagnostics
  diagnostics <- compute_diagnostics(
    parsed_data,
    effects,
    weights,
    transport_result,
    estimator,
    family_y
  )

  # Sensitivity analysis
  sensitivity <- compute_sensitivity_analysis(
    effects,
    diagnostics,
    effect_measure
  )

  # Create result object
  result <- create_result_object(
    effects,
    diff_estimate,
    inference_result,
    diagnostics,
    sensitivity,
    transport_result,
    estimator,
    effect_measure,
    family_y,
    inference_method,
    match.call()
  )

  class(result) <- c("unconf_test", "list")
  result
}

#' Input validation function
#' @param data_rct RCT data
#' @param data_obs Observational data
#' @param formula Model formula
#' @param estimator Chosen estimator
#' @param effect_measure Effect measure
#' @param family_y Outcome family
#' @param alpha Alpha level
#' @param trim Trimming bounds
#' @return List with validation results
validate_inputs <- function(
  data_rct,
  data_obs,
  formula,
  estimator,
  effect_measure,
  family_y,
  alpha,
  trim
) {
  warnings <- character(0)

  # Basic data checks
  if (!is.data.frame(data_rct) || !is.data.frame(data_obs)) {
    return(list(
      valid = FALSE,
      message = "data_rct and data_obs must be data.frames"
    ))
  }

  if (nrow(data_rct) < 50) {
    warnings <- c(
      warnings,
      "Small RCT sample size (n < 50) may lead to unstable estimates"
    )
  }

  if (nrow(data_obs) < 100) {
    warnings <- c(
      warnings,
      "Small observational sample size (n < 100) may lead to unstable estimates"
    )
  }

  # Formula validation
  tryCatch(
    {
      terms(formula)
    },
    error = function(e) {
      return(list(valid = FALSE, message = "Invalid formula specification"))
    }
  )

  # Check required variables
  vars <- all.vars(formula)
  if (length(vars) < 2) {
    return(list(
      valid = FALSE,
      message = "Formula must include outcome and treatment"
    ))
  }

  missing_rct <- setdiff(vars, names(data_rct))
  missing_obs <- setdiff(vars, names(data_obs))

  if (length(missing_rct) > 0) {
    return(list(
      valid = FALSE,
      message = paste(
        "Missing variables in RCT data:",
        paste(missing_rct, collapse = ", ")
      )
    ))
  }

  if (length(missing_obs) > 0) {
    return(list(
      valid = FALSE,
      message = paste(
        "Missing variables in observational data:",
        paste(missing_obs, collapse = ", ")
      )
    ))
  }

  # Check for missing values
  rct_missing <- sum(is.na(data_rct[vars]))
  obs_missing <- sum(is.na(data_obs[vars]))

  if (rct_missing > 0) {
    warnings <- c(
      warnings,
      paste("RCT data has", rct_missing, "missing values")
    )
  }

  if (obs_missing > 0) {
    warnings <- c(
      warnings,
      paste("Observational data has", obs_missing, "missing values")
    )
  }

  # Parameter validation
  if (alpha <= 0 || alpha >= 1) {
    return(list(valid = FALSE, message = "alpha must be between 0 and 1"))
  }

  if (!is.null(trim)) {
    if (
      length(trim) != 2 || trim[1] >= trim[2] || any(trim < 0) || any(trim > 1)
    ) {
      return(list(
        valid = FALSE,
        message = "trim must be c(lower, upper) with 0 <= lower < upper <= 1"
      ))
    }
  }

  # Estimator-specific warnings
  if (estimator == "tmle") {
    if (!requireNamespace("tmle", quietly = TRUE)) {
      return(list(
        valid = FALSE,
        message = "Package 'tmle' required for TMLE estimator"
      ))
    }
  }

  if (estimator == "matching") {
    if (!requireNamespace("MatchIt", quietly = TRUE)) {
      return(list(
        valid = FALSE,
        message = "Package 'MatchIt' required for matching estimator"
      ))
    }
  }

  # Effect measure compatibility
  if (effect_measure %in% c("rr", "or") && family_y != "binomial") {
    return(list(
      valid = FALSE,
      message = "Risk ratio and odds ratio only available for binary outcomes"
    ))
  }

  list(valid = TRUE, warnings = warnings)
}

#' Parse formula and prepare data matrices
#' @param formula Model formula
#' @param data_rct RCT data
#' @param data_obs Observational data
#' @return List with parsed components
parse_formula_and_data <- function(formula, data_rct, data_obs) {
  # Parse formula components
  tt <- terms(formula)
  vars <- all.vars(formula)
  y_name <- vars[1]
  rhs_terms <- attr(tt, "term.labels")
  a_name <- strsplit(rhs_terms[1], ":|\\*")[[1]][1]

  # Create model frames
  mf_rct <- model.frame(formula, data_rct)
  mf_obs <- model.frame(formula, data_obs)

  # Extract variables
  Y_rct <- mf_rct[[y_name]]
  A_rct <- mf_rct[[a_name]]
  Y_obs <- mf_obs[[y_name]]
  A_obs <- mf_obs[[a_name]]

  # Coerce treatment to 0/1
  A_rct <- coerce_binary(A_rct, "Treatment in RCT data")
  A_obs <- coerce_binary(A_obs, "Treatment in observational data")

  # Create covariate matrices
  X_terms <- setdiff(rhs_terms, a_name)
  if (length(X_terms) > 0) {
    form_X <- reformulate(X_terms)
    X_rct <- model.matrix(form_X, data = data_rct)
    X_obs <- model.matrix(form_X, data = data_obs)
  } else {
    X_rct <- matrix(
      1,
      nrow = nrow(data_rct),
      dimnames = list(NULL, "(Intercept)")
    )
    X_obs <- matrix(
      1,
      nrow = nrow(data_obs),
      dimnames = list(NULL, "(Intercept)")
    )
  }

  list(
    Y_rct = Y_rct,
    A_rct = A_rct,
    X_rct = X_rct,
    Y_obs = Y_obs,
    A_obs = A_obs,
    X_obs = X_obs,
    y_name = y_name,
    a_name = a_name,
    data_rct = data_rct,
    data_obs = data_obs,
    formula = formula
  )
}

#' Coerce variable to binary 0/1
#' @param x Variable to coerce
#' @param var_name Variable name for error messages
#' @return Binary 0/1 vector
coerce_binary <- function(x, var_name = "Variable") {
  if (is.logical(x)) {
    return(as.integer(x))
  }

  if (is.factor(x)) {
    x <- droplevels(x)
    if (nlevels(x) != 2) {
      stop(paste(var_name, "must have exactly 2 levels"))
    }
    return(as.integer(x) - 1L)
  }

  if (is.character(x)) {
    unique_vals <- unique(x)
    if (length(unique_vals) != 2) {
      stop(paste(var_name, "must have exactly 2 unique values"))
    }
    return(as.integer(x == sort(unique_vals)[2]))
  }

  if (is.numeric(x)) {
    unique_vals <- sort(unique(x))
    if (length(unique_vals) != 2) {
      stop(paste(var_name, "must have exactly 2 unique values"))
    }
    if (all(unique_vals == c(0, 1))) {
      return(as.integer(x))
    } else {
      return(as.integer(x == unique_vals[2]))
    }
  }

  stop(paste(
    var_name,
    "must be logical, factor, character, or numeric with 2 levels"
  ))
}
