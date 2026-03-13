# effectbridge Issue Backlog

This backlog is prioritized and designed to be converted into GitHub issues.
Use the issue templates in `.github/ISSUE_TEMPLATE/` when creating each one.

## P0 (do first)

### 1) `chore: complete repository rename validation to effectbridge`
- Type: `chore`
- Labels: `maintenance`, `documentation`
- Why:
- Ensure no stale references to the old repository/package identity remain.
- Acceptance criteria:
- `rg -n "unconfoundedr|unconfoudedr"` returns no active references except historical NEWS text where intended.
- Links in README, pkgdown config, and citation files resolve to `effectbridge`.
- CI badges point to the new repository.

### 2) `fix: prevent generated test artifacts from being committed`
- Type: `fix`
- Labels: `bug`, `testing`
- Why:
- Generated artifacts (for example `Rplots.pdf`) create noisy diffs and accidental commits.
- Acceptance criteria:
- Tests do not leave files under `tests/testthat/` after `devtools::test()`.
- `.gitignore` covers known generated test artifacts.
- A regression test or helper ensures temporary graphics devices are cleaned up.

### 3) `ci: enforce R CMD check --as-cran and fail on warnings`
- Type: `ci`
- Labels: `ci`, `quality`
- Why:
- Tight quality gates reduce release risk.
- Acceptance criteria:
- CI workflow runs `R CMD check --as-cran`.
- PRs fail when check warnings/errors occur.
- CI summary shows clear failing step and remediation hint.

## P1 (next wave)

### 4) `test: add edge-case tests for positivity and perfect separation`
- Type: `test`
- Labels: `testing`, `estimators`
- Why:
- These are key failure modes for causal estimators.
- Acceptance criteria:
- New tests cover near-zero/near-one propensity scenarios.
- New tests cover logistic separation fallback behavior.
- Expected warnings/errors are asserted explicitly.

### 5) `test: validate bootstrap and analytical inference against simulation baselines`
- Type: `test`
- Labels: `testing`, `inference`
- Why:
- Inference reliability must be validated empirically.
- Acceptance criteria:
- Simulation tests compare bias/coverage for bootstrap and analytical SE paths.
- Tolerances are documented and justified.
- Tests are deterministic with fixed seeds.

### 6) `refactor: standardize error and warning messages across modules`
- Type: `refactor`
- Labels: `maintenance`, `api`
- Why:
- Consistent diagnostics improve user experience and debugging.
- Acceptance criteria:
- Common validation and messaging helpers are centralized.
- Core workflows use consistent wording and structure.
- Snapshot tests cover critical warning/error messages.

### 7) `docs: expand quick-start examples for common workflows`
- Type: `docs`
- Labels: `documentation`, `examples`
- Why:
- Users need clear paths for typical tasks.
- Acceptance criteria:
- README/vignette include concise examples for AIPW baseline, transport auto mode, and sensitivity interpretation.
- Examples run successfully under CI documentation checks.

## P2 (later)

### 8) `perf: profile bootstrap pipeline and remove avoidable repeated work`
- Type: `perf`
- Labels: `performance`, `inference`
- Why:
- Bootstrap can dominate runtime in real analyses.
- Acceptance criteria:
- Profiling identifies top runtime hotspots.
- At least one measurable optimization is implemented and benchmarked.
- Benchmark notes are documented in PR.

### 9) `feat: add subgroup effect comparison workflow`
- Type: `feat`
- Labels: `enhancement`, `estimators`
- Why:
- Users often need heterogeneity checks by strata.
- Acceptance criteria:
- API for subgroup analysis is documented.
- Output includes subgroup-specific estimates and inference.
- Tests cover happy path and invalid subgroup definitions.

### 10) `chore: create CRAN readiness checklist and release runbook`
- Type: `chore`
- Labels: `release`, `documentation`
- Why:
- Repeatable releases reduce errors and onboarding cost.
- Acceptance criteria:
- Checklist includes preflight checks, version bump steps, NEWS, and tag workflow.
- Release runbook is added to repo and linked from CONTRIBUTING.

## Optional tracking issues

### 11) `meta: 0.3.0 release tracker`
- Type: `meta`
- Labels: `release`, `tracking`
- Why:
- Central place to coordinate milestone scope.
- Acceptance criteria:
- Links to all must-ship issues.
- Explicit in/out scope and target release date.

### 12) `meta: testing quality dashboard`
- Type: `meta`
- Labels: `testing`, `tracking`
- Why:
- Keep visibility on stability over time.
- Acceptance criteria:
- Tracks coverage trend, flaky tests, and open test debt.
- Updated at least once per sprint/release cycle.
