# utils_and_helpers.R - Utility functions and data generators

#' Generate simulated RCT data
#' @param n Sample size
#' @param p_covariates Number of covariates
#' @param treatment_effect True treatment effect
#' @param outcome_type "continuous" or "binary"
#' @param seed Random seed
#' @return Simulated RCT data frame
#' @export
#' @examples
#' \dontrun{
#' rct_data <- generate_rct_data(n = 500, treatment_effect = 1.0)
#' head(rct_data)
#' }
generate_rct_data <- function(
    n = 500,
    p_covariates = 2,
    treatment_effect = 1.0,
    outcome_type = c("continuous", "binary"),
    seed = NULL
) {
    if (!is.null(seed)) {
        set.seed(seed)
    }

    outcome_type <- match.arg(outcome_type)

    # Generate covariates
    X <- matrix(rnorm(n * p_covariates), n, p_covariates)
    colnames(X) <- paste0("X", 1:p_covariates)

    # Randomized treatment assignment
    A <- rbinom(n, 1, 0.5)

    # Generate potential outcomes
    if (outcome_type == "continuous") {
        # Linear outcome model
        Y0 <- 0.5 + X %*% rep(0.5, p_covariates) + rnorm(n)
        Y1 <- Y0 + treatment_effect
        Y <- ifelse(A == 1, Y1, Y0)
    } else {
        # Binary outcome (logistic model)
        logit_Y0 <- -1 + X %*% rep(0.5, p_covariates)
        logit_Y1 <- logit_Y0 + treatment_effect
        Y0 <- rbinom(n, 1, plogis(logit_Y0))
        Y1 <- rbinom(n, 1, plogis(logit_Y1))
        Y <- ifelse(A == 1, Y1, Y0)
    }

    data.frame(Y = Y, A = A, X)
}

#' Generate simulated observational data
#' @param n Sample size
#' @param p_covariates Number of covariates
#' @param treatment_effect True treatment effect
#' @param confounding_strength Strength of confounding (0 = none, 1 = strong)
#' @param covariate_shift Whether to introduce covariate shift
#' @param outcome_type "continuous" or "binary"
#' @param seed Random seed
#' @return Simulated observational data frame
#' @export
#' @examples
#' \dontrun{
#' obs_data <- generate_obs_data(n = 1000, confounding_strength = 0.5,
#'                               covariate_shift = TRUE)
#' head(obs_data)
#' }
generate_obs_data <- function(
    n = 1000,
    p_covariates = 2,
    treatment_effect = 1.0,
    confounding_strength = 0.3,
    covariate_shift = FALSE,
    outcome_type = c("continuous", "binary"),
    seed = NULL
) {
    if (!is.null(seed)) {
        set.seed(seed)
    }

    outcome_type <- match.arg(outcome_type)

    # Generate covariates (potentially with shift)
    if (covariate_shift) {
        # Introduce shift in covariate distribution
        X <- matrix(rnorm(n * p_covariates, mean = 0.3), n, p_covariates)
    } else {
        X <- matrix(rnorm(n * p_covariates), n, p_covariates)
    }
    colnames(X) <- paste0("X", 1:p_covariates)

    # Generate unmeasured confounder
    U <- rnorm(n)

    # Non-randomized treatment assignment (depends on X and U)
    treatment_model <- -0.5 +
        X %*% rep(0.3, p_covariates) +
        confounding_strength * U
    A <- rbinom(n, 1, plogis(treatment_model))

    # Generate potential outcomes
    if (outcome_type == "continuous") {
        # Outcome depends on X, U, and treatment
        Y0 <- 0.5 +
            X %*% rep(0.5, p_covariates) +
            confounding_strength * U +
            rnorm(n)
        Y1 <- Y0 + treatment_effect
        Y <- ifelse(A == 1, Y1, Y0)
    } else {
        # Binary outcome (logistic model)
        logit_Y0 <- -1 + X %*% rep(0.5, p_covariates) + confounding_strength * U
        logit_Y1 <- logit_Y0 + treatment_effect
        Y0 <- rbinom(n, 1, plogis(logit_Y0))
        Y1 <- rbinom(n, 1, plogis(logit_Y1))
        Y <- ifelse(A == 1, Y1, Y0)
    }

    data.frame(Y = Y, A = A, X)
}

#' Simulate multiple datasets for power analysis
#' @param n_sim Number of simulations
#' @param n_rct RCT sample size
#' @param n_obs Observational sample size
#' @param treatment_effect_rct True treatment effect in RCT
#' @param treatment_effect_obs True treatment effect in observational study
#' @param confounding_strength Strength of confounding in observational study
#' @param alpha Significance level
#' @param seed Random seed
#' @return Power analysis results
#' @export
#' @examples
#' \dontrun{
#' power_results <- simulate_power_analysis(
#'   n_sim = 100, n_rct = 500, n_obs = 1000,
#'   treatment_effect_rct = 1.0, treatment_effect_obs = 1.2
#' )
#' print(power_results)
#' }
simulate_power_analysis <- function(
    n_sim = 500,
    n_rct = 500,
    n_obs = 1000,
    treatment_effect_rct = 1.0,
    treatment_effect_obs = 1.0,
    confounding_strength = 0.3,
    alpha = 0.05,
    seed = NULL
) {
    if (!is.null(seed)) {
        set.seed(seed)
    }

    # Storage for results
    p_values <- numeric(n_sim)
    differences <- numeric(n_sim)
    coverage <- logical(n_sim)

    true_difference <- treatment_effect_obs - treatment_effect_rct

    cat(sprintf("Running %d simulations...\n", n_sim))

    for (i in 1:n_sim) {
        if (i %% 100 == 0) {
            cat(sprintf("Simulation %d/%d\n", i, n_sim))
        }

        # Generate data
        rct_data <- generate_rct_data(
            n = n_rct,
            treatment_effect = treatment_effect_rct,
            seed = i
        )

        obs_data <- generate_obs_data(
            n = n_obs,
            treatment_effect = treatment_effect_obs,
            confounding_strength = confounding_strength,
            seed = i + n_sim
        )

        # Run test (with reduced bootstrap for speed)
        tryCatch(
            {
                result <- unconfoundedness_test(
                    rct_data,
                    obs_data,
                    formula = Y ~ A + X1 + X2,
                    estimator = "aipw",
                    transport = "none", # Skip transport for power analysis
                    B = 250, # Reduced bootstrap
                    validate = FALSE, # Skip validation for speed
                    seed = i + 2 * n_sim
                )

                p_values[i] <- result$inference$p_value
                differences[i] <- result$estimates$difference

                # Check coverage
                ci <- result$inference$confidence_interval
                coverage[i] <- (true_difference >= ci[1] &&
                    true_difference <= ci[2])
            },
            error = function(e) {
                p_values[i] <- NA
                differences[i] <- NA
                coverage[i] <- FALSE
            }
        )
    }

    # Compute power and other statistics
    valid_sims <- !is.na(p_values)
    n_valid <- sum(valid_sims)

    power <- mean(p_values[valid_sims] < alpha)
    type_1_error <- if (true_difference == 0) power else NA
    coverage_rate <- mean(coverage[valid_sims])

    # Bias and MSE
    bias <- mean(differences[valid_sims]) - true_difference
    mse <- mean((differences[valid_sims] - true_difference)^2)

    list(
        settings = list(
            n_sim = n_sim,
            n_valid = n_valid,
            n_rct = n_rct,
            n_obs = n_obs,
            treatment_effect_rct = treatment_effect_rct,
            treatment_effect_obs = treatment_effect_obs,
            true_difference = true_difference,
            confounding_strength = confounding_strength,
            alpha = alpha
        ),
        results = list(
            power = power,
            type_1_error = type_1_error,
            coverage_rate = coverage_rate,
            bias = bias,
            mse = mse,
            mean_difference = mean(differences[valid_sims]),
            sd_difference = sd(differences[valid_sims])
        ),
        raw_results = list(
            p_values = p_values[valid_sims],
            differences = differences[valid_sims],
            coverage = coverage[valid_sims]
        )
    )
}

#' Multiple comparison correction for subgroup analyses
#' @param p_values Vector of p-values
#' @param method Correction method
#' @return Adjusted p-values
#' @export
adjust_p_values <- function(
    p_values,
    method = c("bonferroni", "holm", "BH", "BY")
) {
    method <- match.arg(method)
    p.adjust(p_values, method = method)
}

#' Compute sample size for desired power
#' @param power Desired power (default 0.8)
#' @param alpha Significance level (default 0.05)
#' @param effect_size Expected effect size difference
#' @param ratio_obs_to_rct Ratio of observational to RCT sample size
#' @return Required sample sizes
#' @export
#' @examples
#' \dontrun{
#' sample_sizes <- compute_required_sample_size(
#'   power = 0.8, effect_size = 0.3, ratio_obs_to_rct = 2
#' )
#' print(sample_sizes)
#' }
compute_required_sample_size <- function(
    power = 0.8,
    alpha = 0.05,
    effect_size = 0.5,
    ratio_obs_to_rct = 2
) {
    # Approximate sample size calculation based on two-sample t-test
    # This is approximate - exact calculation would require simulation

    z_alpha <- qnorm(1 - alpha / 2)
    z_beta <- qnorm(power)

    # Assuming equal variances (approximation)
    # For two independent samples with different sizes
    pooled_factor <- 1 + 1 / ratio_obs_to_rct

    n_rct <- 2 * ((z_alpha + z_beta)^2) * pooled_factor / (effect_size^2)
    n_obs <- n_rct * ratio_obs_to_rct

    list(
        n_rct = ceiling(n_rct),
        n_obs = ceiling(n_obs),
        total_n = ceiling(n_rct + n_obs),
        assumptions = list(
            equal_variances = TRUE,
            two_sided_test = TRUE,
            normal_approximation = TRUE
        )
    )
}

#' Create latex table from results
#' @param x unconf_test object or list of results
#' @param file Output file path (NULL for console output)
#' @param caption Table caption
#' @return LaTeX table code
#' @export
create_latex_table <- function(
    x,
    file = NULL,
    caption = "Unconfoundedness Test Results"
) {
    if (inherits(x, "unconf_test")) {
        # Single result
        latex_code <- format_single_result_latex(x, caption)
    } else if (
        is.list(x) && all(sapply(x, function(obj) inherits(obj, "unconf_test")))
    ) {
        # Multiple results
        latex_code <- format_multiple_results_latex(x, caption)
    } else {
        stop("Input must be unconf_test object or list of unconf_test objects")
    }

    if (!is.null(file)) {
        writeLines(latex_code, file)
        cat(sprintf("LaTeX table written to %s\n", file))
    }

    cat(latex_code, sep = "\n")
    invisible(latex_code)
}

#' Format single result for LaTeX
#' @param x unconf_test object
#' @param caption Table caption
#' @return LaTeX code
format_single_result_latex <- function(x, caption) {
    c(
        "\\begin{table}[htbp]",
        "\\centering",
        paste0("\\caption{", caption, "}"),
        "\\begin{tabular}{lcccc}",
        "\\toprule",
        "Dataset & Estimate & SE & 95\\% CI & P-value \\\\",
        "\\midrule",
        sprintf("RCT & %.3f & --- & --- & --- \\\\", x$estimates$rct),
        sprintf(
            "Observational & %.3f & --- & --- & --- \\\\",
            x$estimates$observational
        ),
        sprintf(
            "Difference & %.3f & %.3f & [%.3f, %.3f] & %s \\\\",
            x$estimates$difference,
            x$inference$standard_error,
            x$inference$confidence_interval[1],
            x$inference$confidence_interval[2],
            format_p_value_latex(x$inference$p_value)
        ),
        "\\bottomrule",
        "\\end{tabular}",
        sprintf(
            "\\label{tab:unconfoundedness_%s}",
            format(Sys.time(), "%Y%m%d_%H%M")
        ),
        "\\end{table}"
    )
}

#' Format multiple results for LaTeX
#' @param results List of unconf_test objects
#' @param caption Table caption
#' @return LaTeX code
format_multiple_results_latex <- function(results, caption) {
    # Extract key information from each result
    n_results <- length(results)
    result_names <- names(results) %||% paste("Analysis", 1:n_results)

    # Create table rows
    table_rows <- character(n_results)
    for (i in 1:n_results) {
        x <- results[[i]]
        table_rows[i] <- sprintf(
            "%s & %.3f & %.3f & %.3f & [%.3f, %.3f] & %s \\\\",
            result_names[i],
            x$estimates$rct,
            x$estimates$observational,
            x$estimates$difference,
            x$inference$confidence_interval[1],
            x$inference$confidence_interval[2],
            format_p_value_latex(x$inference$p_value)
        )
    }

    c(
        "\\begin{table}[htbp]",
        "\\centering",
        paste0("\\caption{", caption, "}"),
        "\\begin{tabular}{lcccccc}",
        "\\toprule",
        "Analysis & RCT Effect & Obs Effect & Difference & 95\\% CI & P-value \\\\",
        "\\midrule",
        table_rows,
        "\\bottomrule",
        "\\end{tabular}",
        sprintf(
            "\\label{tab:unconfoundedness_multiple_%s}",
            format(Sys.time(), "%Y%m%d_%H%M")
        ),
        "\\end{table}"
    )
}

#' Format p-value for LaTeX
#' @param p P-value
#' @return Formatted p-value
format_p_value_latex <- function(p) {
    if (is.na(p)) {
        "---"
    } else if (p < 0.001) {
        "$<$ 0.001"
    } else {
        sprintf("%.3f", p)
    }
}

#' Export results to CSV
#' @param x unconf_test object or list of results
#' @param file Output CSV file path
#' @return Invisibly returns the data frame
#' @export
export_to_csv <- function(x, file) {
    if (inherits(x, "unconf_test")) {
        df <- x$summary
    } else if (
        is.list(x) && all(sapply(x, function(obj) inherits(obj, "unconf_test")))
    ) {
        # Combine multiple results
        df_list <- lapply(names(x) %||% seq_along(x), function(i) {
            result <- if (is.numeric(i)) x[[i]] else x[[i]]
            result_df <- result$summary
            result_df$Analysis <- if (is.character(i)) {
                i
            } else {
                paste("Analysis", i)
            }
            result_df
        })
        df <- do.call(rbind, df_list)
    } else {
        stop("Input must be unconf_test object or list of unconf_test objects")
    }

    write.csv(df, file, row.names = FALSE)
    cat(sprintf("Results exported to %s\n", file))
    invisible(df)
}

#' Create diagnostic report
#' @param x unconf_test object
#' @param file Output file path (HTML or PDF)
#' @param format Output format ("html" or "pdf")
#' @return File path of created report
#' @export
create_diagnostic_report <- function(
    x,
    file = "unconfoundedness_report.html",
    format = c("html", "pdf")
) {
    format <- match.arg(format)

    if (!requireNamespace("rmarkdown", quietly = TRUE)) {
        stop("Package 'rmarkdown' required for creating reports")
    }

    # Create temporary R Markdown file
    rmd_content <- create_rmd_content(x)
    temp_rmd <- tempfile(fileext = ".Rmd")
    writeLines(rmd_content, temp_rmd)

    # Render report
    output_format <- if (format == "html") {
        rmarkdown::html_document(theme = "default", toc = TRUE)
    } else {
        rmarkdown::pdf_document(toc = TRUE)
    }

    rmarkdown::render(
        temp_rmd,
        output_format = output_format,
        output_file = file,
        quiet = TRUE,
        envir = new.env()
    )

    cat(sprintf("Diagnostic report created: %s\n", file))
    invisible(file)
}

#' Create R Markdown content for diagnostic report
#' @param x unconf_test object
#' @return R Markdown content as character vector
create_rmd_content <- function(x) {
    c(
        "---",
        "title: \"Unconfoundedness Test Diagnostic Report\"",
        paste0("date: \"", Sys.Date(), "\""),
        "output:",
        "  html_document:",
        "    theme: default",
        "    toc: true",
        "    toc_float: true",
        "---",
        "",
        "```{r setup, include=FALSE}",
        "knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)",
        "```",
        "",
        "## Summary",
        "",
        "This report presents the results of an unconfoundedness test comparing treatment effects between an RCT-like dataset and an observational dataset.",
        "",
        "### Key Results",
        "",
        sprintf("- **RCT Effect**: %.4f", x$estimates$rct),
        sprintf("- **Observational Effect**: %.4f", x$estimates$observational),
        sprintf("- **Difference**: %.4f", x$estimates$difference),
        sprintf(
            "- **95%% Confidence Interval**: [%.4f, %.4f]",
            x$inference$confidence_interval[1],
            x$inference$confidence_interval[2]
        ),
        sprintf("- **P-value**: %s", format_p_value(x$inference$p_value)),
        "",
        "### Interpretation",
        "",
        interpret_test_result(x$inference$p_value, x$estimates$difference),
        "",
        "## Methods",
        "",
        sprintf("- **Estimator**: %s", x$methods$estimator),
        sprintf("- **Effect Measure**: %s", x$methods$effect_measure),
        sprintf("- **Outcome Family**: %s", x$methods$outcome_family),
        sprintf("- **Inference Method**: %s", x$methods$inference_method),
        sprintf("- **Transport Method**: %s", x$methods$transport_method),
        "",
        "## Diagnostics",
        "",
        "### Overlap Assessment",
        "",
        sprintf(
            "- **Overlap Quality**: %s",
            x$diagnostics$overlap$overlap_quality
        ),
        sprintf(
            "- **Treatment Prevalence**: RCT=%.1f%%, Obs=%.1f%%",
            100 * x$diagnostics$overlap$treatment_prevalence$rct,
            100 * x$diagnostics$overlap$treatment_prevalence$obs
        ),
        "",
        "### Sample Size Information",
        "",
        sprintf(
            "- **RCT Sample Size**: %d",
            x$diagnostics$sample_size$raw_n$rct
        ),
        sprintf(
            "- **Observational Sample Size**: %d",
            x$diagnostics$sample_size$raw_n$obs
        ),
        sprintf(
            "- **Minimum Group Size**: %d",
            x$diagnostics$sample_size$minimum_group_size
        ),
        "",
        "## Sensitivity Analysis",
        "",
        if (!is.null(x$sensitivity$e_value)) {
            c(
                sprintf("- **E-value**: %.2f", x$sensitivity$e_value$e_value),
                sprintf(
                    "- **Interpretation**: %s",
                    x$sensitivity$e_value$interpretation
                )
            )
        } else {
            "E-value not computed."
        },
        "",
        "## Recommendations",
        "",
        paste("-", generate_recommendations(x)),
        "",
        "---",
        "",
        "*Report generated by the unconfoundedr package.*"
    )
}

#' Batch analysis for multiple datasets
#' @param rct_data_list List of RCT datasets
#' @param obs_data_list List of observational datasets
#' @param formula Model formula
#' @param estimator Estimator to use
#' @param ... Additional arguments passed to unconfoundedness_test
#' @return List of results
#' @export
batch_analysis <- function(
    rct_data_list,
    obs_data_list,
    formula,
    estimator = "aipw",
    ...
) {
    if (length(rct_data_list) != length(obs_data_list)) {
        stop("RCT and observational data lists must have same length")
    }

    n_analyses <- length(rct_data_list)
    results <- vector("list", n_analyses)
    names(results) <- names(rct_data_list) %||% paste("Analysis", 1:n_analyses)

    cat(sprintf("Running %d analyses...\n", n_analyses))

    for (i in 1:n_analyses) {
        cat(sprintf("Analysis %d/%d: %s\n", i, n_analyses, names(results)[i]))

        tryCatch(
            {
                results[[i]] <- unconfoundedness_test(
                    data_rct = rct_data_list[[i]],
                    data_obs = obs_data_list[[i]],
                    formula = formula,
                    estimator = estimator,
                    ...
                )
            },
            error = function(e) {
                cat(sprintf("Error in analysis %d: %s\n", i, e$message))
                results[[i]] <- list(error = e$message)
            }
        )
    }

    class(results) <- c("batch_unconf_test", "list")
    results
}

#' Print method for batch results
#' @param x batch_unconf_test object
#' @param ... Additional arguments
#' @export
print.batch_unconf_test <- function(x, ...) {
    cat("# Batch Unconfoundedness Test Results\n")
    cat("=" * 50, "\n\n")

    n_analyses <- length(x)
    n_successful <- sum(sapply(x, function(obj) !is.null(obj$estimates)))
    n_failed <- n_analyses - n_successful

    cat(sprintf("Total analyses: %d\n", n_analyses))
    cat(sprintf("Successful: %d\n", n_successful))
    cat(sprintf("Failed: %d\n", n_failed))
    cat("\n")

    if (n_successful > 0) {
        # Create summary table
        successful_results <- x[sapply(x, function(obj) {
            !is.null(obj$estimates)
        })]

        summary_table <- do.call(
            rbind,
            lapply(names(successful_results), function(name) {
                result <- successful_results[[name]]
                data.frame(
                    Analysis = name,
                    RCT_Effect = result$estimates$rct,
                    Obs_Effect = result$estimates$observational,
                    Difference = result$estimates$difference,
                    P_Value = result$inference$p_value,
                    Significant = result$inference$p_value < 0.05,
                    stringsAsFactors = FALSE
                )
            })
        )

        cat("## Summary of Results\n")
        print(summary_table, row.names = FALSE, digits = 4)

        # Overall statistics
        cat(sprintf(
            "\nSignificant results: %d/%d (%.1f%%)\n",
            sum(summary_table$Significant),
            n_successful,
            100 * mean(summary_table$Significant)
        ))
    }

    invisible(x)
}

#' Validate package installation and dependencies
#' @return Logical indicating if all dependencies are available
#' @export
validate_installation <- function() {
    cat("Validating unconfoundedr installation...\n\n")

    # Core dependencies
    core_deps <- c("stats", "utils")
    core_available <- sapply(core_deps, requireNamespace, quietly = TRUE)

    # Optional dependencies
    optional_deps <- c("tmle", "MatchIt", "energy", "parallel", "rmarkdown")
    optional_available <- sapply(
        optional_deps,
        requireNamespace,
        quietly = TRUE
    )

    # Report results
    cat("Core dependencies:\n")
    for (i in seq_along(core_deps)) {
        status <- if (core_available[i]) "✓" else "✗"
        cat(sprintf("  %s %s\n", status, core_deps[i]))
    }

    cat("\nOptional dependencies:\n")
    for (i in seq_along(optional_deps)) {
        status <- if (optional_available[i]) "✓" else "✗"
        functionality <- switch(
            optional_deps[i],
            "tmle" = "(for TMLE estimator)",
            "MatchIt" = "(for matching estimator)",
            "energy" = "(for multivariate shift detection)",
            "parallel" = "(for parallel bootstrap)",
            "rmarkdown" = "(for diagnostic reports)",
            ""
        )
        cat(sprintf("  %s %s %s\n", status, optional_deps[i], functionality))
    }

    all_core_available <- all(core_available)

    if (all_core_available) {
        cat("\n✓ Installation validated successfully!\n")
    } else {
        cat(
            "\n✗ Installation validation failed. Please install missing core dependencies.\n"
        )
    }

    # Missing optional suggestions
    missing_optional <- optional_deps[!optional_available]
    if (length(missing_optional) > 0) {
        cat(sprintf(
            "\nTo enable all features, install: %s\n",
            paste(missing_optional, collapse = ", ")
        ))
    }

    invisible(all_core_available)
}

# Additional utility operators and functions
`%||%` <- function(a, b) if (is.null(a)) b else a
`%nin%` <- function(x, table) !x %in% table

#' Safely extract element from list
#' @param x List or vector
#' @param element Element name or index
#' @param default Default value if element doesn't exist
#' @return Element value or default
safe_extract <- function(x, element, default = NA) {
    if (is.list(x) && element %in% names(x)) {
        x[[element]]
    } else if (!is.list(x) && is.numeric(element) && element <= length(x)) {
        x[element]
    } else {
        default
    }
}
