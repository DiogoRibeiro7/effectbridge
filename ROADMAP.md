# effectbridge Roadmap

This roadmap focuses on stabilizing the package for release readiness while
improving reliability, usability, and documentation quality.

## Now (next 2-4 weeks)

1. Release consistency and package identity
- Ensure package rename to `effectbridge` is complete across docs, tests, CI, and metadata.
- Regenerate and validate `NAMESPACE` and `man/*.Rd` from roxygen comments.
- Verify `DESCRIPTION`, `CITATION.cff`, and `inst/CITATION` version alignment.

2. CI reliability and quality gates
- Enforce `R CMD check --as-cran` passing in GitHub Actions.
- Add lintr/style checks in CI with clear failure messages.
- Improve test coverage reporting and set a minimum threshold.

3. Test artifact hygiene
- Stop committing generated files from test runs (for example `tests/testthat/Rplots.pdf`).
- Ensure test helpers create artifacts in temp dirs and clean up after execution.

## Next (1-2 months)

1. Estimator correctness and diagnostics hardening
- Add focused tests for edge cases: positivity violations, perfect separation, and tiny samples.
- Validate robustness of inference methods (bootstrap vs analytical) under known simulations.
- Improve user-facing diagnostics for overlap and model failure modes.

2. API ergonomics
- Clarify argument contracts and defaults for `unconfoundedness_test()`.
- Standardize error/warning messages across estimators and transport modules.
- Add high-signal examples for common workflows (AIPW baseline, transport auto mode, sensitivity analysis).

3. Performance and scalability
- Profile hotspots in bootstrap and transport routines.
- Reduce repeated model fits where results can be reused.
- Document recommended settings for fast exploratory runs vs full inference runs.

## Later (2-4 months)

1. Release and distribution
- Prepare CRAN submission checklist and revdep-style compatibility checks.
- Tag stable release, publish release notes, and verify pkgdown deployment.

2. Extensions
- Add support for additional effect types and subgroup analyses.
- Add optional integration adapters for common causal workflows (for example `targets` pipelines).

3. Governance
- Define semantic versioning/release cadence.
- Add lightweight maintainer checklist for release, triage, and docs updates.
