# unconfoundedness_test.R
# R >= 4.1

#' @title Test unconfoundedness with optional transport weighting (auto/KS/energy)
#' @description
#' Compare a marginal treatment effect ω in an RCT-like dataset with the same
#' estimand in an observational dataset. Supports:
#' - Estimators: IPW or AIPW (doubly robust).
#' - Inference: bootstrap CI for Δ = ω_obs − ω_rct; Wald-style z test.
#' - Transport weighting:
#'   * "none": no reweighting;
#'   * "rct_to_obs": always reweight RCT to OBS covariates via density ratio logistic;
#'   * "auto": run shift tests (KS/energy) and reweight only if shift is detected.
#'
#' @param data_rct,data_obs data.frame with same Y/A definitions and covariates.
#' @param formula model formula `Y ~ A + X1 + X2`. A must be binary and the first RHS term.
#' @param estimator "aipw" (default) or "ipw".
#' @param stabilize logical; stabilized IPW (default TRUE).
#' @param trim length-2 numeric in [0,1] for weight trimming (e.g., c(0.01,0.99)); NULL disables.
#' @param family_y "gaussian" or "binomial" (for AIPW outcome regression and checks).
#' @param transport "none", "rct_to_obs", or "auto" (default "none").
#' @param auto_method which tests to use under `transport="auto"`: "both" (default), "ks", or "energy".
#' @param auto_alpha significance level for shift tests (default 0.01).
#' @param auto_energy_R number of permutations for energy test (default 199). Ignored if `energy` not installed.
#' @param B bootstrap replicates for Δ CI (default 1000).
#' @param alpha CI tail (default 0.05 for 95% CI).
#' @param seed RNG seed or NULL.
#' @return list with estimates, CI, p-values, settings, and diagnostics (incl. auto decision).
#' @examples
#' \dontrun{
#' set.seed(1)
#' gen_rct <- function(n){
#'   X1 <- rnorm(n); X2 <- rbinom(n,1,0.4)
#'   A  <- rbinom(n,1,0.5)
#'   Y0 <- 0.5 + 0.5*X1 + 0.3*X2 + rnorm(n)
#'   Y1 <- Y0 + 1.0
#'   data.frame(Y=ifelse(A==1,Y1,Y0),A,X1,X2)
#' }
#' gen_obs <- function(n){
#'   X1 <- rnorm(n, 0.4); X2 <- rbinom(n,1,0.7)  # shift
#'   A  <- rbinom(n,1,plogis(-0.2 + 0.8*X1 + 0.6*X2))
#'   U  <- rnorm(n)
#'   Y0 <- 0.5 + 0.5*X1 + 0.3*X2 + 0.3*U + rnorm(n)
#'   Y1 <- Y0 + 1.0
#'   data.frame(Y=ifelse(A==1,Y1,Y0),A,X1,X2)
#' }
#' d_rct <- gen_rct(800); d_obs <- gen_obs(2000)
#' out <- unconfoundedness_test(d_rct, d_obs, Y ~ A + X1 + X2,
#'   estimator="aipw", family_y="gaussian",
#'   transport="auto", auto_method="both", auto_alpha=0.01,
#'   B=500, seed=42)
#' print(out); out$diagnostics$auto
#' }
unconfoundedness_test <- function(
  data_rct,
  data_obs,
  formula,
  estimator = c("aipw", "ipw"),
  stabilize = TRUE,
  trim = c(0.01, 0.99),
  family_y = c("gaussian", "binomial"),
  transport = c("none", "rct_to_obs", "auto"),
  auto_method = c("both", "ks", "energy"),
  auto_alpha = 0.01,
  auto_energy_R = 199L,
  B = 1000L,
  alpha = 0.05,
  seed = NULL
) {
  # ---------------- Args & basic checks ----------------
  stopifnot(is.data.frame(data_rct), is.data.frame(data_obs))
  estimator <- match.arg(estimator)
  family_y <- match.arg(family_y)
  transport <- match.arg(transport)
  auto_method <- match.arg(auto_method)
  stopifnot(is.numeric(alpha) && alpha > 0 && alpha < 1)
  stopifnot(is.numeric(auto_alpha) && auto_alpha > 0 && auto_alpha < 1)
  stopifnot(is.null(seed) || (is.numeric(seed) && length(seed) == 1))
  if (!is.null(trim)) {
    stopifnot(
      is.numeric(trim),
      length(trim) == 2,
      trim[1] < trim[2],
      all(trim >= 0),
      all(trim <= 1)
    )
  }

  # ---------------- Parse formula ----------------
  tt <- terms(formula)
  vars <- all.vars(formula)
  if (length(vars) < 2L) {
    stop("Formula must be Y ~ A + Xs.")
  }
  y_name <- vars[1L]
  rhs_terms <- attr(tt, "term.labels")
  if (length(rhs_terms) < 1L) {
    stop("Formula RHS must include treatment A (binary).")
  }
  a_name <- strsplit(rhs_terms[1L], ":|\\*")[[1L]][1L] # first symbol treated as A

  req <- unique(c(vars, a_name, y_name))
  if (!all(req %in% names(data_rct))) {
    stop("Missing formula columns in data_rct.")
  }
  if (!all(req %in% names(data_obs))) {
    stop("Missing formula columns in data_obs.")
  }

  # ---------------- Coerce A to {0,1} ----------------
  coerce_A <- function(df, a) {
    v <- df[[a]]
    if (is.factor(v)) {
      v <- droplevels(v)
    }
    if (is.logical(v)) {
      v <- as.integer(v)
    }
    if (is.character(v)) {
      uv <- unique(v)
      if (length(uv) != 2) {
        stop("Character A must have 2 levels.")
      }
      v <- as.integer(v == uv[2L])
    }
    if (is.factor(v)) {
      v <- as.integer(v == levels(v)[2L])
    }
    if (!is.numeric(v)) {
      stop("A must be numeric/logic/factor/character 2-level.")
    }
    if (!all(v %in% c(0, 1))) {
      u <- sort(unique(v))
      stop(sprintf("A must be 0/1. Found: %s", paste(u, collapse = ", ")))
    }
    df[[a]] <- v
    df
  }
  data_rct <- coerce_A(data_rct, a_name)
  data_obs <- coerce_A(data_obs, a_name)

  # ---------------- Frames & checks ----------------
  mf_rct <- model.frame(formula, data_rct)
  mf_obs <- model.frame(formula, data_obs)
  Y_r <- mf_rct[[y_name]]
  A_r <- mf_rct[[a_name]]
  Y_o <- mf_obs[[y_name]]
  A_o <- mf_obs[[a_name]]
  if (family_y == "binomial") {
    if (!all(Y_r %in% c(0, 1)) || !all(Y_o %in% c(0, 1))) {
      stop("With family_y='binomial', Y must be 0/1 in both datasets.")
    }
  } else {
    if (!is.numeric(Y_r) || !is.numeric(Y_o)) {
      stop("With family_y='gaussian', Y must be numeric.")
    }
  }

  # ---------------- Covariates (exclude A) ----------------
  rhs_form <- reformulate(rhs_terms, response = y_name)
  X_terms <- setdiff(attr(terms(rhs_form), "term.labels"), a_name)
  form_X <- if (length(X_terms)) reformulate(X_terms) else ~1
  X_r <- model.matrix(form_X, data = data_rct)
  X_o <- model.matrix(form_X, data = data_obs)

  # ---------------- Auto shift detection (if requested) ----------------
  auto_diag <- NULL
  transport_applied <- transport

  # KS: univariate two-sample across each column (except intercept)
  auto_detect_shift_ks <- function(Xr, Xo, alpha) {
    cols <- setdiff(colnames(Xr), "(Intercept)")
    if (length(cols) == 0L) {
      return(list(
        available = TRUE,
        pvals = NULL,
        padj = NULL,
        reject_any = FALSE,
        rejected_cols = character(0)
      ))
    }
    pvals <- vapply(
      cols,
      function(cn) {
        suppressWarnings(stats::ks.test(Xr[, cn], Xo[, cn])$p.value)
      },
      numeric(1)
    )
    padj <- stats::p.adjust(pvals, method = "BH")
    reject <- padj < alpha
    list(
      available = TRUE,
      pvals = pvals,
      padj = padj,
      reject_any = any(reject),
      rejected_cols = names(padj)[reject]
    )
  }

  # Energy: multivariate 2-sample E-test with permutations
  auto_detect_shift_energy <- function(Xr, Xo, alpha, R) {
    have_energy <- requireNamespace("energy", quietly = TRUE)
    if (!have_energy) {
      return(list(available = FALSE, p_value = NA_real_, reject = FALSE, R = R))
    }
    Z <- rbind(Xr, Xo)
    sizes <- c(nrow(Xr), nrow(Xo))
    et <- energy::eqdist.etest(Z, sizes = sizes, R = as.integer(R))
    list(
      available = TRUE,
      p_value = et$p.value,
      reject = (et$p.value < alpha),
      R = R
    )
  }

  if (transport == "auto") {
    ks_res <- if (auto_method %in% c("both", "ks")) {
      auto_detect_shift_ks(X_r, X_o, auto_alpha)
    } else {
      NULL
    }
    en_res <- if (auto_method %in% c("both", "energy")) {
      auto_detect_shift_energy(X_r, X_o, auto_alpha, auto_energy_R)
    } else {
      NULL
    }

    detected <- FALSE
    used <- character(0)
    if (!is.null(ks_res) && ks_res$available) {
      if (ks_res$reject_any) {
        detected <- TRUE
        used <- c(used, "ks")
      }
    }
    if (!is.null(en_res) && en_res$available) {
      if (en_res$reject) {
        detected <- TRUE
        used <- c(used, "energy")
      }
    }
    # If method="energy" but not available, fallback to "ks" only (document in diagnostics)
    if (auto_method == "energy" && (is.null(en_res) || !en_res$available)) {
      # Nothing else to run; detected stays FALSE unless ks was also requested.
      used <- c(used, "energy_unavailable")
    }

    transport_applied <- if (detected) "rct_to_obs" else "none"
    auto_diag <- list(
      requested_method = auto_method,
      alpha = auto_alpha,
      detected = detected,
      tests_used = used,
      ks = ks_res,
      energy = en_res
    )
  }

  # ---------------- Transport weights (applied or not) ----------------
  density_ratio_weights <- function(Xr, Xo) {
    stopifnot(is.matrix(Xr), is.matrix(Xo))
    Z <- rbind(Xr, Xo)
    lab <- c(rep(0L, nrow(Xr)), rep(1L, nrow(Xo))) # 0=RCT, 1=OBS
    df <- data.frame(S = lab, Z)
    if ("(Intercept)" %in% names(df)) {
      df[["(Intercept)"]] <- NULL
    }
    fit <- stats::glm(S ~ ., data = df, family = stats::binomial())
    pr <- stats::predict(fit, newdata = data.frame(Z = Xr), type = "response")
    pr <- pmin(pmax(as.numeric(pr), 1e-6), 1 - 1e-6)
    w <- pr / (1 - pr)
    w / mean(w)
  }
  w_rct <- if (transport_applied == "rct_to_obs") {
    density_ratio_weights(X_r, X_o)
  } else {
    rep(1, nrow(X_r))
  }
  w_obs <- rep(1, nrow(X_o))

  # ---------------- Estimators (weighted) ----------------
  fit_ps <- function(A, X, samp_w) {
    df <- data.frame(A = A, X)
    if ("(Intercept)" %in% names(df)) {
      df[["(Intercept)"]] <- NULL
    }
    suppressWarnings({
      fit <- stats::glm(
        A ~ .,
        data = df,
        family = stats::binomial(),
        weights = samp_w
      )
    })
    p <- pmin(pmax(as.numeric(stats::fitted(fit)), 1e-6), 1 - 1e-6)
    list(p = p, fit = fit)
  }
  fit_outcome_models <- function(Y, A, X, family, samp_w) {
    fam <- if (family == "binomial") stats::binomial() else stats::gaussian()
    df1 <- data.frame(Y = Y[A == 1], X[A == 1, , drop = FALSE])
    w1 <- samp_w[A == 1]
    df0 <- data.frame(Y = Y[A == 0], X[A == 0, , drop = FALSE])
    w0 <- samp_w[A == 0]
    if ("(Intercept)" %in% names(df1)) {
      df1[["(Intercept)"]] <- NULL
      df0[["(Intercept)"]] <- NULL
    }
    f1 <- suppressWarnings(stats::glm(
      Y ~ .,
      data = df1,
      family = fam,
      weights = w1
    ))
    f0 <- suppressWarnings(stats::glm(
      Y ~ .,
      data = df0,
      family = fam,
      weights = w0
    ))
    type <- "response"
    m1 <- as.numeric(stats::predict(f1, newdata = data.frame(X), type = type))
    m0 <- as.numeric(stats::predict(f0, newdata = data.frame(X), type = type))
    list(m1 = m1, m0 = m0, fit1 = f1, fit0 = f0)
  }
  stab_num <- function(A) {
    pA <- mean(A)
    A * pA + (1 - A) * (1 - pA)
  }
  ate_ipw <- function(Y, A, e, samp_w, stabilize = TRUE, trim = NULL) {
    numer <- if (stabilize) stab_num(A) else 1
    w_t <- samp_w * numer * A / e
    w_c <- samp_w * numer * (1 - A) / (1 - e)
    w_all <- w_t + w_c
    if (!is.null(trim)) {
      ql <- stats::quantile(w_all, trim[1])
      qh <- stats::quantile(w_all, trim[2])
      keep <- (w_all >= ql) & (w_all <= qh)
      Y <- Y[keep]
      w_t <- w_t[keep]
      w_c <- w_c[keep]
    }
    mu1 <- sum(w_t * Y) / sum(w_t)
    mu0 <- sum(w_c * Y) / sum(w_c)
    list(tau = mu1 - mu0)
  }
  ate_aipw <- function(Y, A, e, m1, m0, samp_w, stabilize = TRUE, trim = NULL) {
    base_w <- samp_w * (A / e + (1 - A) / (1 - e))
    if (!is.null(trim)) {
      ql <- stats::quantile(base_w, trim[1])
      qh <- stats::quantile(base_w, trim[2])
      keep <- (base_w >= ql) & (base_w <= qh)
      Y <- Y[keep]
      A <- A[keep]
      e <- e[keep]
      m1 <- m1[keep]
      m0 <- m0[keep]
      samp_w <- samp_w[keep]
    }
    numer <- if (stabilize) stab_num(A) else 1
    pseudo <- m1 -
      m0 +
      numer * (A * (Y - m1) / e - (1 - A) * (Y - m0) / (1 - e))
    W <- samp_w
    Sw <- sum(W)
    tau <- sum(W * pseudo) / Sw
    IF <- (W / Sw) * (pseudo - tau) # weighted influence for mean
    list(tau = tau, IF = IF)
  }
  one_effect <- function(df, family, est, samp_w) {
    mf <- model.frame(formula, df)
    Y <- mf[[y_name]]
    A <- mf[[a_name]]
    X <- model.matrix(form_X, df)
    ps <- fit_ps(A, X, samp_w)
    e <- ps$p
    if (est == "ipw") {
      ipw <- ate_ipw(
        Y,
        A,
        e,
        samp_w = samp_w,
        stabilize = stabilize,
        trim = trim
      )
      list(
        ate = ipw$tau,
        var = NA_real_,
        IF = NULL,
        details = list(ps_model = ps$fit, e = e, method = "ipw")
      )
    } else {
      outs <- fit_outcome_models(Y, A, X, family, samp_w)
      aipw <- ate_aipw(
        Y,
        A,
        e,
        outs$m1,
        outs$m0,
        samp_w = samp_w,
        stabilize = stabilize,
        trim = trim
      )
      v <- sum(aipw$IF^2) # var(Σ IF_i) since IF already scaled by W/Sw
      list(
        ate = aipw$tau,
        var = v,
        IF = aipw$IF,
        details = list(
          ps_model = ps$fit,
          out1 = outs$fit1,
          out0 = outs$fit0,
          e = e,
          m1 = outs$m1,
          m0 = outs$m0,
          method = "aipw"
        )
      )
    }
  }

  if (!is.null(seed)) {
    set.seed(as.integer(seed))
  }

  # ---------------- Point estimates ----------------
  est_r <- one_effect(data_rct, family_y, estimator, samp_w = w_rct)
  est_o <- one_effect(data_obs, family_y, estimator, samp_w = w_obs)
  diff_hat <- est_o$ate - est_r$ate

  # ---------------- Bootstrap Δ (decision is fixed across reps) ----------------
  if (B < 200L) {
    warning("Low bootstrap reps; consider B >= 500 for stable CIs.")
  }
  boot_diffs <- numeric(B)
  n_r <- nrow(data_rct)
  n_o <- nrow(data_obs)
  for (b in seq_len(B)) {
    idx_r <- sample.int(n_r, n_r, replace = TRUE)
    idx_o <- sample.int(n_o, n_o, replace = TRUE)
    Xr_b <- X_r[idx_r, , drop = FALSE]
    Xo_b <- X_o[idx_o, , drop = FALSE]
    wr_b <- if (transport_applied == "rct_to_obs") {
      density_ratio_weights(Xr_b, Xo_b)
    } else {
      rep(1, length(idx_r))
    }
    eb_r <- one_effect(
      data_rct[idx_r, , drop = FALSE],
      family_y,
      estimator,
      samp_w = wr_b
    )$ate
    eb_o <- one_effect(
      data_obs[idx_o, , drop = FALSE],
      family_y,
      estimator,
      samp_w = rep(1, length(idx_o))
    )$ate
    boot_diffs[b] <- eb_o - eb_r
  }
  ci <- stats::quantile(
    boot_diffs,
    probs = c(alpha / 2, 1 - alpha / 2),
    names = FALSE
  )

  # ---------------- Wald z (SE choice) ----------------
  # If we applied transport weighting OR estimator==ipw -> bootstrap SE.
  # Else (no transport, AIPW) -> IF-based SE.
  if (transport_applied != "none" || estimator == "ipw") {
    se_diff <- stats::sd(boot_diffs)
  } else {
    se_diff <- sqrt(est_r$var + est_o$var)
  }
  z <- as.numeric(diff_hat / se_diff)
  p_wald <- 2 * (1 - stats::pnorm(abs(z)))

  # ---------------- Diagnostics ----------------
  diag_overlap <- function(A, e) {
    rng <- range(e)
    extreme <- mean(e < 0.01 | e > 0.99)
    list(ps_min = rng[1], ps_max = rng[2], prop_extreme = extreme)
  }
  diag_r <- diag_overlap(
    model.frame(formula, data_rct)[[a_name]],
    est_r$details$e
  )
  diag_o <- diag_overlap(
    model.frame(formula, data_obs)[[a_name]],
    est_o$details$e
  )

  diag_transport <- NULL
  if (transport_applied == "rct_to_obs") {
    ESS <- (sum(w_rct)^2) / sum(w_rct^2)
    diag_transport <- list(
      type = "rct_to_obs",
      w_mean = mean(w_rct),
      w_min = min(w_rct),
      w_max = max(w_rct),
      w_q = as.numeric(stats::quantile(w_rct, c(0.01, 0.5, 0.99))),
      ESS = ESS
    )
  }

  out <- list(
    estimator = estimator,
    family_y = family_y,
    transport = transport,
    transport_applied = transport_applied, # "none" or "rct_to_obs"
    rct_effect = est_r$ate,
    obs_effect = est_o$ate,
    diff_obs_minus_rct = diff_hat,
    bootstrap_CI = c(lower = ci[1], upper = ci[2]),
    wald = list(z = z, p_value = p_wald, se_diff = se_diff),
    diagnostics = list(
      rct_ps = diag_r,
      obs_ps = diag_o,
      transport = diag_transport,
      auto = auto_diag,
      n_rct = nrow(data_rct),
      n_obs = nrow(data_obs),
      trim = trim,
      stabilize = stabilize
    ),
    call = match.call()
  )

  out$summary <- data.frame(
    Estimator = estimator,
    Transport_Requested = transport,
    Transport_Applied = transport_applied,
    Effect_RCT = est_r$ate,
    Effect_OBS = est_o$ate,
    Diff_OBS_minus_RCT = diff_hat,
    CI_lower = ci[1],
    CI_upper = ci[2],
    Wald_z = z,
    Wald_p = p_wald,
    n_RCT = nrow(data_rct),
    n_OBS = nrow(data_obs),
    row.names = NULL
  )

  class(out) <- c("unconf_test", "list")
  out
}

#' Print method for unconf_test objects
#'
#' @param x An object of class `unconf_test`.
#' @param ... Passed to methods (unused).
#' @return Invisibly returns `x`.
#' @export
print.unconf_test <- function(x, ...) {
  cat("# Unconfoundedness test (RCT vs OBS)\n")
  cat(sprintf(
    "- Estimator: %s; Outcome family: %s; Transport requested: %s; Applied: %s\n",
    x$estimator,
    x$family_y,
    x$transport,
    x$transport_applied
  ))
  cat(sprintf(
    "- ω̂_RCT = %.6f; ω̂_OBS = %.6f; Δ = %.6f\n",
    x$rct_effect,
    x$obs_effect,
    x$diff_obs_minus_rct
  ))
  cat(sprintf(
    "- Bootstrap CI for Δ: [%.6f, %.6f]\n",
    x$bootstrap_CI[1],
    x$bootstrap_CI[2]
  ))
  cat(sprintf(
    "- Wald test: z = %.3f, p = %.4g (H0: Δ = 0)\n",
    x$wald$z,
    x$wald$p_value
  ))
  cat("- Positivity diagnostics (PS in [min,max], extreme PS <1% or >99%):\n")
  cat(sprintf(
    "  RCT: [%.3f, %.3f], extreme=%.2f%%;  OBS: [%.3f, %.3f], extreme=%.2f%%\n",
    x$diagnostics$rct_ps$ps_min,
    x$diagnostics$rct_ps$ps_max,
    100 * x$diagnostics$rct_ps$prop_extreme,
    x$diagnostics$obs_ps$ps_min,
    x$diagnostics$obs_ps$ps_max,
    100 * x$diagnostics$obs_ps$prop_extreme
  ))
  if (!is.null(x$diagnostics$transport)) {
    tr <- x$diagnostics$transport
    cat("- Transport (RCT->OBS) weights summary:\n")
    cat(sprintf(
      "  mean=%.3f, min=%.3f, q01=%.3f, median=%.3f, q99=%.3f, max=%.3f, ESS=%.1f\n",
      tr$w_mean,
      tr$w_min,
      tr$w_q[1],
      tr$w_q[2],
      tr$w_q[3],
      tr$w_max,
      tr$ESS
    ))
  }
  if (!is.null(x$diagnostics$auto)) {
    ad <- x$diagnostics$auto
    cat(sprintf(
      "- Auto shift detection: detected=%s, method=%s, alpha=%.3f\n",
      ad$detected,
      ad$requested_method,
      ad$alpha
    ))
    if (!is.null(ad$ks)) {
      ks <- ad$ks
      if (ks$available && length(ks$padj)) {
        bad <- if (length(ks$rejected_cols)) {
          paste(ks$rejected_cols, collapse = ", ")
        } else {
          "none"
        }
        cat(sprintf(
          "  KS: any_reject=%s; rejected_cols=%s\n",
          ks$reject_any,
          bad
        ))
      } else {
        cat("  KS: not run or no covariates beyond intercept.\n")
      }
    }
    if (!is.null(ad$energy)) {
      en <- ad$energy
      if (en$available) {
        cat(sprintf("  Energy: p=%.4g (R=%d)\n", en$p_value, en$R))
      } else {
        cat("  Energy: package 'energy' not available; skipped.\n")
      }
    }
  }
  invisible(x)
}

# Not exported; small helper for print safety
`%||%` <- function(a, b) {
  if (is.null(a)) b else a
}
