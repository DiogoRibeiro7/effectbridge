# effectbridge

> **Test unconfoundedness by comparing treatment effects from RCT-like and observational datasets**

[![R-CMD-check](https://img.shields.io/github/actions/workflow/status/diogoribeiro7/effectbridge/R-CMD-check.yaml?label=R-CMD-check)](https://github.com/diogoribeiro7/effectbridge/actions/workflows/R-CMD-check.yaml) 
[![pkgdown](https://img.shields.io/github/actions/workflow/status/diogoribeiro7/effectbridge/pkgdown.yaml?label=pkgdown)](https://github.com/diogoribeiro7/effectbridge/actions/workflows/pkgdown.yaml) 
[![Codecov](https://codecov.io/gh/diogoribeiro7/effectbridge/branch/develop/graph/badge.svg)](https://codecov.io/gh/diogoribeiro7/effectbridge)
[![CRAN status](https://www.r-pkg.org/badges/version/effectbridge)](https://CRAN.R-project.org/package=effectbridge)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)

## Overview

The `effectbridge` package provides a comprehensive toolkit for testing the **unconfoundedness assumption** (ignorability) by comparing marginal treatment effects estimated from an RCT-like dataset with those from an observational dataset. When these effects differ significantly, it suggests potential unmeasured confounding in the observational analysis.

### Why This Matters

Before trusting causal inferences from observational data, it's crucial to assess whether the **"no unmeasured confounders"** assumption holds. This package implements a practical approach:

1. **Estimate the same causal estimand** in both RCT and observational data
2. **Compare the estimates** - large differences suggest unmeasured confounding
3. **Account for transportability** - populations may differ between studies
4. **Provide comprehensive diagnostics** - assess overlap, balance, and robustness

## Key Features

### 🔧 **Multiple Estimators**
- **AIPW (Augmented IPW)**: Doubly robust estimation
- **IPW**: Inverse probability weighting
- **TMLE**: Targeted maximum likelihood estimation  
- **G-computation**: Outcome regression approach
- **Matching**: Propensity score matching

### 📊 **Effect Measures**
- **Risk Difference (RD)**: μ₁ - μ₀
- **Risk Ratio (RR)**: μ₁ / μ₀  
- **Odds Ratio (OR)**: [μ₁/(1-μ₁)] / [μ₀/(1-μ₀)]

### 🌉 **Transportability Weighting**
- **Manual**: Always reweight RCT to observational population
- **Auto-detection**: Use KS tests and/or energy statistics to detect covariate shift
- **Comprehensive diagnostics**: Effective sample sizes, weight distributions

### 📈 **Robust Inference**
- **Bootstrap confidence intervals** with parallel processing
- **Analytical standard errors** using influence functions
- **Robust sandwich estimators** for model misspecification

### 🔍 **Comprehensive Diagnostics**
- **Overlap assessment**: Propensity score distributions and positivity
- **Balance evaluation**: Standardized mean differences
- **Model diagnostics**: Convergence, fit statistics
- **Sample size analysis**: Effective sample sizes, minimum group sizes

### 🛡️ **Sensitivity Analysis**
- **E-values**: Quantify robustness to unmeasured confounding  
- **Fragility assessment**: Identify vulnerable aspects of analysis
- **Robustness evaluation**: Model specification, estimation method, data quality

### 📊 **Visualization & Reporting**
- **Diagnostic plots**: Effects, overlap, balance, weights, bootstrap distributions
- **Automated reports**: HTML/PDF reports with comprehensive diagnostics
- **Export capabilities**: LaTeX tables, CSV summaries
- **Batch analysis**: Process multiple dataset comparisons

## Installation

```r
# From GitHub (development version)
devtools::install_github("DiogoRibeiro7/effectbridge")

# Optional dependencies for full functionality
install.packages(c("tmle", "MatchIt", "energy", "parallel", "rmarkdown"))

# Validate installation
library(effectbridge)
validate_installation()
```

## Quick Start

```r
library(effectbridge)

# Generate example data
set.seed(42)
rct_data <- generate_rct_data(n = 800, treatment_effect = 1.0)
obs_data <- generate_obs_data(n = 2000, treatment_effect = 1.2, 
                              confounding_strength = 0.3, covariate_shift = TRUE)

# Run comprehensive unconfoundedness test  
result <- unconfoundedness_test(
  data_rct = rct_data,
  data_obs = obs_data, 
  formula = Y ~ A + X1 + X2,           # A must be first RHS term
  estimator = "aipw",                  # Doubly robust
  effect_measure = "rd",               # Risk difference
  transport = "auto",                  # Auto-detect covariate shift
  auto_method = "both",                # KS + energy tests
  inference_method = "bootstrap",      # Bootstrap inference
  B = 1000,                           # Bootstrap replicates
  parallel = TRUE,                     # Parallel processing
  seed = 42
)

# View results
print(result)
summary(result)

# Diagnostic plots
plot(result, type = "effects")       # Effect comparison
plot(result, type = "overlap")       # Propensity score overlap  
plot(result, type = "balance")       # Covariate balance
plot(result, type = "bootstrap")     # Bootstrap distribution
```

## Advanced Usage

### Multiple Estimators Comparison

```r
estimators <- c("aipw", "ipw", "tmle", "gcomp")
results <- list()

for (est in estimators) {
  results[[est]] <- unconfoundedness_test(
    rct_data, obs_data, Y ~ A + X1 + X2,
    estimator = est, B = 500, seed = 42
  )
}

# Compare results
sapply(results, function(r) r$estimates$difference)
export_to_csv(results, "comparison_results.csv")
```

### Power Analysis

```r
# Simulate power across different scenarios
power_analysis <- simulate_power_analysis(
  n_sim = 1000,
  n_rct = 500, n_obs = 2000,
  treatment_effect_rct = 1.0,
  treatment_effect_obs = 1.3,  # Confounding bias
  confounding_strength = 0.4,
  seed = 123
)

print(power_analysis)
# Power: 0.856, Coverage: 0.943, Bias: 0.297
```

### Batch Analysis

```r
# Multiple studies
rct_studies <- list(
  study_a = rct_data_a,
  study_b = rct_data_b,  
  study_c = rct_data_c
)

obs_studies <- list(
  study_a = obs_data_a,
  study_b = obs_data_b,
  study_c = obs_data_c  
)

batch_results <- batch_analysis(
  rct_studies, obs_studies, 
  formula = Y ~ A + X1 + X2 + X3,
  estimator = "aipw",
  transport = "auto"
)

print(batch_results)
create_latex_table(batch_results, file = "batch_results.tex")
```

### Generate Diagnostic Report

```r
# Comprehensive HTML report
create_diagnostic_report(result, file = "diagnostics.html", format = "html")

# PDF report  
create_diagnostic_report(result, file = "diagnostics.pdf", format = "pdf")
```

## Interpreting Results

### Key Output Components

```r
result$estimates
# $rct: 1.023
# $observational: 1.347  
# $difference: 0.324

result$inference
# $p_value: 0.003
# $confidence_interval: [0.112, 0.536]
# $standard_error: 0.108

result$sensitivity$e_value
# $e_value: 2.31
# $interpretation: "Moderately robust to unmeasured confounding"
```

### Decision Framework

| P-value | Difference | Interpretation | Recommendation |
|---------|------------|----------------|----------------|
| > 0.05  | Small      | ✅ Supports unconfoundedness | Proceed with observational analysis |
| < 0.05  | Large      | ⚠️ Suggests confounding | Investigate sources, sensitivity analysis |
| < 0.01  | Very large | 🚨 Strong evidence of bias | Consider alternative approaches |

### Diagnostic Quality Indicators

- **Overlap Quality**: "Good" or "Excellent" preferred
- **Balance Quality**: SMD < 0.1 preferred  
- **E-value**: > 2.0 indicates reasonable robustness
- **Effective Sample Size**: > 50% of original preferred

## Methodological Details

### Estimand

The package targets the **Average Treatment Effect (ATE)**:
```
τ = E[Y(1) - Y(0)]
```

For binary outcomes with different effect measures:
- **RD**: P(Y=1|A=1) - P(Y=1|A=0)  
- **RR**: P(Y=1|A=1) / P(Y=1|A=0)
- **OR**: [P(Y=1|A=1)/(1-P(Y=1|A=1))] / [P(Y=1|A=0)/(1-P(Y=1|A=0))]

### Test Statistic

**Null Hypothesis**: H₀: τ_obs = τ_rct (unconfoundedness holds)

**Test Statistic**: 
```
Z = (τ̂_obs - τ̂_rct) / SE(τ̂_obs - τ̂_rct)
```

### Transportability

When populations differ, the package can reweight the RCT to match the observational population using density ratio estimation:

```
w(x) = P(S=obs|X=x) / P(S=rct|X=x)
```

## Dependencies

### Required
- R (≥ 4.1.0)
- stats, utils, graphics

### Optional (for full functionality)
- **tmle**: TMLE estimator
- **MatchIt**: Matching methods  
- **energy**: Multivariate shift detection
- **parallel**: Parallel bootstrap
- **rmarkdown**: Diagnostic reports

## Contributing

We welcome contributions! Please see:
- [Contributing Guidelines](.github/CONTRIBUTING.md)
- [Code of Conduct](.github/CODE_OF_CONDUCT.md)

### Development Workflow

```r
# Setup development environment
devtools::load_all()
devtools::document()
devtools::test()
devtools::check()

# Add new tests
usethis::use_test("new_feature")

# Update documentation
devtools::document()
pkgdown::build_site()
```

## Citation

If you use this package, please cite:

```bibtex
@software{effectbridge2025,
  title = {effectbridge: Test Unconfoundedness by Comparing RCT and Observational Effects},
  author = {Diogo Ribeiro},
  year = {2025},
  version = {0.2.0},
  url = {https://github.com/DiogoRibeiro7/effectbridge},
}
```

## Related Work

- **Hartman & Hidalgo (2018)**: "An Assessment of the Augmented Inverse Propensity Weighted Estimator"
- **D'Amour et al. (2017)**: "Overlap in observational studies with high-dimensional covariates"  
- **VanderWeele & Ding (2017)**: "Sensitivity Analysis in Observational Research"

## License

MIT License. See [LICENSE](LICENSE) for details.

## Support

- 📖 **Documentation**: [https://diogoribeiro7.github.io/packages/effectbridge/](https://diogoribeiro7.github.io/packages/effectbridge/)
- 🐛 **Bug Reports**: [GitHub Issues](https://github.com/DiogoRibeiro7/effectbridge/issues)  
- 💬 **Questions**: [GitHub Discussions](https://github.com/DiogoRibeiro7/effectbridge/discussions)
- 📧 **Email**: diogo.debastos.ribeiro@gmail.com

---

*Built with ❤️ for robust causal inference*

