# =============================================================================
# Script:       models.R
# Author:       Logan A. Sauers
# Repository:   github.com/SauersALogan/<repo-name>
# Description:  Shared modelling functions. run_lmer_parallel() fits one lmer
#               per item (taxon/gene) in parallel, extracts pairwise emmeans
#               contrasts and fit diagnostics, and applies BH correction.
#               Definitions only; the future plan is set by the calling script.
# Output:       run_lmer_parallel()
# =============================================================================
library(lme4)
library(emmeans)
library(MuMIn)
library(dplyr)
library(future.apply)
library(progressr)

## ---- Differential expression function (can use for abundance also) ----
run_lmer_parallel <- function(data_list, formula, emmeans_term,
    fdr_group, singular_action = c("null", "flag_keep", "flag_exclude"),
    output_cols = list(), inner_split = NULL) {

    singular_action <- match.arg(singular_action)

    cat("\n=== Running parallel lmer models ===\n")
    cat("Items to model:", length(data_list), "\n")
    cat("Formula:", deparse(formula), "\n")
    cat("Singular action:", singular_action, "\n")

    start_time <- Sys.time()

    with_progress({
        p <- progressor(steps = length(data_list))

        results <- future_lapply(names(data_list), function(item_name) {
            d <- data_list[[item_name]]

            sub_list <- if (!is.null(inner_split)) {
                split(d, d[[inner_split]])
            } else {
                list(d)
            }

            sub_results <- lapply(sub_list, function(d) {

                tryCatch({
                    m <- lmer(formula, data = d,
                        control = lmerControl(optimizer = "bobyqa"))

                    singular_flag <- isSingular(m)

                    if (singular_flag && singular_action == "null")
                        return(NULL)
                    if (singular_flag && singular_action == "flag_exclude")
                        return(NULL)

                    r2_vals <- r.squaredGLMM(m)
                    residuals_m <- residuals(m)
                    shapiro_p <- if (length(residuals_m) >= 3 &&
                        length(residuals_m) <= 5000) {
                        shapiro.test(residuals_m)$p.value
                    } else NA

                    em   <- emmeans(m, emmeans_term)
                    cf   <- contrast(em, method = "pairwise", adjust = "none")
                    summ <- as.data.frame(summary(cf))

                    summ <- summ %>%
                        dplyr::mutate(r_squared = r2_vals[1, "R2m"],
                            adj_r_squared = r2_vals[1, "R2c"],
                            shapiro_p = shapiro_p,
                            singular_fit = singular_flag,
                            n_obs = nrow(d))

                    for (col_name in names(output_cols)) {
                        summ[[col_name]] <- eval(output_cols[[col_name]],
                            envir = d)
                    }

                    summ

                }, error = function(e) {
                    message("Error in ", item_name, ": ", e$message)
                    return(NULL)
                })
            })

            p()
            return(dplyr::bind_rows(sub_results))

        }, future.seed = TRUE,
           future.packages = c("lme4", "emmeans", "MuMIn", "dplyr"))
    })

    results_df <- dplyr::bind_rows(results)

    cat("Time:", round(difftime(Sys.time(),
        start_time, units = "mins"), 2), "minutes\n")
    if (length(output_cols) > 0) {
        cat("Items with results:",
            n_distinct(results_df[[names(output_cols)[1]]]), "\n")
    }

    results_df <- results_df %>%
        dplyr::group_by(across(all_of(fdr_group))) %>%
        dplyr::mutate(p_adjust = p.adjust(p.value, method = "BH"),
            significant = p_adjust < 0.05) %>%
        dplyr::ungroup()

    cat("Significant results (FDR < 0.05):",
        sum(results_df$significant, na.rm = TRUE), "\n")

    return(results_df)
}