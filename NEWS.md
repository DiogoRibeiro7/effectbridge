# effectbridge 0.2.1 (development)

## Improvements

- Refactored monolithic source into focused modules for maintainability:

  - `unconfoundedness_test.R`: main function, validation, formula parsing.
  - `estimators.R`: IPW, AIPW, TMLE, G-computation, matching implementations.
  - `inference_and_diagnostics.R`: inference methods and diagnostic suite.
  - `methods.R`: S3 print/summary methods.
  - `plotting.R`: all visualization functions.
  - `transport.R`: transport weighting and covariate shift detection.
  - `utils_and_helpers.R`: data generators, batch analysis, reporting utilities.

- Organized pkgdown reference into categorized sections (Core Testing,
  Data Generation, Power & Sample Size, S3 Methods, Reporting, Utilities).

- Standardized code style to 2-space indentation across all source files.

- Removed unused internal helpers (`%nin%`, `safe_extract`).

- Expanded test suite to 7 files covering estimators, validation, transport,
  methods/plotting, and utilities.

- Added `.lintr` configuration for consistent linting across contributors.

- Rewrote `getting-started` vignette: fixed malformed code chunks (broken
  line breaks and extra backticks) and properly formatted all R examples.

- Expanded `test-sim.R` with thorough data generator tests (sample sizes,
  covariate counts, binary outcomes, reproducibility, covariate shift, and
  confounding strength).

- Updated `CITATION.cff`: synced email with DESCRIPTION, added ORCID,
  updated abstract to reflect all five estimators, fixed package URL.

- Added `Language: en` field to DESCRIPTION.

## Maintenance

- Improved repository hygiene by cleaning and expanding `.gitignore`:

  - removed duplicate patterns,
  - added common R build/check artifacts (`*.Rcheck/`, `*.tar.gz`, `check/`),
  - added additional local/session ignores (`.Ruserdata`, `renv/library/`, `renv/python/`),
  - added local tooling cache ignores (`.pytest_cache/`, `.claude/`),
  - kept documentation/build output ignores (`docs/`, `pkgdown/`, `.quarto/`).

- Fixed incorrect variable reference in README example code.

- Fixed Codecov badge to reference `develop` branch instead of `main`.

- Removed broken `DEVELOPMENT.md` link from README (content already in
  CONTRIBUTING.md).

# effectbridge 0.2.0

## Fixes

- Fixed critical bug in propensity score fitting:

  - Added `build_ps_formula()` to safely construct formulas.
  - Added `extract_covariates()` to pull covariates from a user formula.
  - `fit_propensity_model()` now:

    - drops intercept-only columns,
    - falls back to `A ~ 1` when no predictors are available,
    - clamps fitted propensities to avoid 0/1 weights,
    - handles perfect separation by reverting to marginal treatment probability.

- Fixed print method for batch results:

  - Replaced invalid `cat("=" * 50)` with `cat(strrep("=", 50))`.

- Removed stale reference to non-existent `inference_and_diagnostics.R`.

## Improvements

- Stronger input validation with clearer error messages for formulas, weights, and covariates.
- More robust handling of empty or misspecified formulas in test cases.
- Internal helpers are now keyworded as `internal` to clarify user-facing API.

# effectbridge 0.1.0

- Initial CRAN-style release.
- Core function `unconfoundedness_test()` implemented with support for:

  - IPW, AIPW, TMLE, G-computation, and matching estimators.
  - Effect measures: risk difference, risk ratio, odds ratio.
  - Inference methods: bootstrap, analytical, and robust SE.
  - Optional transport weighting with auto-detection (KS, energy tests).
  - Basic diagnostics and sensitivity analysis.

