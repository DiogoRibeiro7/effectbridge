# Contributing to unconfoundedr

Thank you for your interest in contributing! 🎉

## Quick Start

1. **Fork** the repository
2. **Clone** your fork: `git clone https://github.com/yourusername/unconfoundedr.git`
3. **Create** a branch: `git checkout -b feature-name`
4. **Make** your changes
5. **Test**: `devtools::test()`
6. **Check**: `devtools::check()`
7. **Commit** and **push**
8. **Submit** a pull request

## Development Setup

```r
# Install development dependencies
devtools::install_dev_deps()

# Load package for development
devtools::load_all()

# Run tests
devtools::test()

# Check package
devtools::check()

# Build documentation
devtools::document()
```

## Types of Contributions

### 🐛 Bug Reports

Use the bug report template and include:

- Reproducible example using `reprex`
- Session info from `sessionInfo()`
- Expected vs actual behavior

### ✨ Feature Requests

Use the feature request template and describe:

- Problem statement
- Proposed solution
- Use cases and examples

### 🔧 Code Contributions

- Follow the existing code style
- Add tests for new functionality
- Update documentation
- Run `devtools::check()` before submitting

### 📝 Documentation

- Fix typos or improve clarity
- Add examples to function documentation
- Write or improve vignettes

## Code Style

We follow the [tidyverse style guide](https://style.tidyverse.org/):

```r
# Good
result <- unconfoundedness_test(
  data_rct = rct_data,
  data_obs = obs_data,
  formula = Y ~ A + X1 + X2
)

# Bad
result<-unconfoundedness_test(data_rct=rct_data,data_obs=obs_data,formula=Y~A+X1+X2)
```

## Testing

- Write tests for all new functions
- Place tests in `tests/testthat/`
- Name test files `test-functionality.R`
- Use descriptive test names

```r
test_that("unconfoundedness_test works with basic input", {
  # Test code here
  expect_equal(result$estimate, expected_value, tolerance = 1e-3)
})
```

## Documentation

- Document all exported functions with roxygen2
- Include `@param`, `@return`, and `@examples`
- Keep examples simple and fast-running

```r
#' Test unconfoundedness assumption
#'
#' @param data_rct RCT-like dataset
#' @param data_obs Observational dataset  
#' @return Object of class "unconf_test"
#' @export
#' @examples
#' # Simple example
#' result <- unconfoundedness_test(rct_data, obs_data, Y ~ A + X)
```

## Commit Messages

Use conventional commits:

- `feat:` new features
- `fix:` bug fixes
- `docs:` documentation changes
- `test:` adding tests
- `refactor:` code refactoring

Example: `feat: add TMLE estimator support`

## Questions?

- 💬 [GitHub Discussions](https://github.com/DiogoRibeiro7/unconfoundedr/discussions)
- 📧 Email: <diogo.debastos.ribeiro@gmail.com>
- 📖 Check the [documentation](https://diogoribeiro7.github.io/packages/unconfoundedr/)

We appreciate all contributions, large and small! 🙏
