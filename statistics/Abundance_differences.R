# =============================================================================
# Script:       Abundance_differences.R
# Author:       Logan A. Sauers
# Repository:   github.com/SauersALogan/<repo-name>
# Description:  Per-family mixed models testing pup Genotype and mom_genotype
#               (each crossed with sampling_age, random intercept for Cage),
#               fitted in parallel. BH correction across families per term.
# Dependencies: data_structuring.R (clr route), abundances.R (abs route)
# Input:        Data (from data_structuring.R) or Data_absolute_Both.rds
# Output:       family_tests_<data_type>.rds / .csv in OUT_DIR
# =============================================================================
library(phyloseq)
library(microbiome)
library(lme4)
library(car)
library(dplyr)
library(future)
library(future.apply)
library(progressr)

# ---- Paths (OUT_DIR set in .Renviron) ----
out_dir <- Sys.getenv("out_dir")
if (out_dir == "") stop("out_dir not set in .Renviron.")

# ---- Settings ----
data_type <- "abs_log10"   # "clr" = relative (CLR), "abs_log10" = absolute (log10)
min_prev  <- 0.1       # prevalence filter, applied before transformation
n_workers <- 8
pup_group_col <- "DOB"     # "Cage" or "DOB": level mice are nested within (pup model)
mouse_col     <- "mouse_index"   # [TO CONFIRM] mouse ID column


source("statistics/functions/data_structuring.R")
 
# =====================================================================
# Data preparation (family level, filtered, transformed)
# =====================================================================
remove_spikeins <- function(ps) {
    tx <- as.data.frame(tax_table(ps))
    spk <- rownames(tx)[grepl("imtechella|allobacillus", tx$species,
        ignore.case = TRUE)]
    prune_taxa(!taxa_names(ps) %in% spk, ps)
}
 
if (data_type == "clr") {
    ps <- subset_samples(Data, Purpose == "assemblage")
    ps <- remove_spikeins(ps)
    ps <- tax_glom(ps, taxrank = "family")
    ps <- filter_taxa(ps, function(x) mean(x > 0) >= min_prev, prune = TRUE)
    ps <- microbiome::transform(ps, "clr")
} else if (data_type == "abs_log10") {
    ps <- readRDS(file.path(out_dir, "Data_absolute_Both.rds"))
    ps <- subset_samples(ps, Purpose == "assemblage")
    ps <- tax_glom(ps, taxrank = "family")
    ps <- filter_taxa(ps, function(x) mean(x > 0) >= min_prev, prune = TRUE)
    # [TO DECIDE] zeros remain after filtering; log10 needs a pseudocount.
    # Here: half the smallest non-zero value in the table.
    abs_m <- as(otu_table(ps), "matrix")
    pc <- min(abs_m[abs_m > 0]) / 2
    ps <- transform_sample_counts(ps, function(x) log10(x + pc))
} else {
    stop("Set data_type to 'clr' or 'abs_log10'.")
}
 
long_df <- psmelt(ps) %>%
    dplyr::rename(response = Abundance)
 
family_list <- split(long_df, long_df$family)
cat("Families to test:", length(family_list), "\n")
 
# =====================================================================
# Pup Genotype: per-family lmer (parallel)
# =====================================================================
test_family_pup <- function(d, family_name) {
    # Workers are fresh R sessions: set sum contrasts here so Type III is valid
    options(contrasts = c("contr.sum", "contr.poly"))
 
    tryCatch({
        m <- lmer(response ~ Genotype * sampling_age + (1 | Cage), data = d,
            control = lmerControl(optimizer = "bobyqa"))
        a <- as.data.frame(Anova(m, type = "III"))
        a$term <- rownames(a)
        a <- a[a$term != "(Intercept)", ]
        data.frame(family   = family_name,
                   term     = a$term,
                   chisq    = a$Chisq,
                   df       = a$Df,
                   p        = a$`Pr(>Chisq)`,
                   singular = isSingular(m),
                   conv_msg = paste(m@optinfo$conv$lme4$messages,
                       collapse = "; "),
                   n_obs    = nrow(d))
    }, error = function(e) {
        message("Error in ", family_name, ": ", e$message)
        NULL
    })
}
 
plan(multisession, workers = n_workers)
if (interactive()) handlers(global = TRUE)
handlers("progress")
 
start_time <- Sys.time()
 
with_progress({
    p <- progressor(steps = length(family_list))
    pup_results <- future_lapply(names(family_list), function(fam) {
        r <- test_family_pup(family_list[[fam]], fam)
        p()
        r
    }, future.seed = TRUE,
       future.packages = c("lme4", "car", "dplyr"))
})
 
plan(sequential)
 
cat("Time:", round(difftime(Sys.time(), start_time, units = "mins"), 2),
    "minutes\n")
 
pup_df <- dplyr::bind_rows(pup_results) %>%
    dplyr::group_by(term) %>%
    dplyr::mutate(p_adjust = p.adjust(p, method = "BH"),
        significant = p_adjust < 0.05) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(term, p_adjust)
 
cat("Pup Genotype - significant (FDR < 0.05) by term:\n")
print(pup_df %>% dplyr::group_by(term) %>%
    dplyr::summarise(n_sig = sum(significant, na.rm = TRUE),
        n_singular = sum(singular),
        n_conv_warning = sum(conv_msg != ""), .groups = "drop"))
 
# =====================================================================
# mom_genotype: cage-level test on cage means
# =====================================================================
mom_df <- dplyr::bind_rows(lapply(names(family_list), function(fam) {
    cm <- family_list[[fam]] %>%
        dplyr::group_by(Cage, mom_genotype) %>%
        dplyr::summarise(response = mean(response), .groups = "drop")
    if (var(cm$response) < 1e-10) return(NULL)   # no variation, untestable
    fit <- summary(lm(response ~ mom_genotype, data = cm))
    data.frame(family  = fam,
               n_cages = nrow(cm),
               t       = fit$coefficients[2, 3],
               p       = fit$coefficients[2, 4])
})) %>%
    dplyr::mutate(p_adjust = p.adjust(p, method = "BH"),
        significant = p_adjust < 0.05) %>%
    dplyr::arrange(p_adjust)
 
cat("mom_genotype - significant (FDR < 0.05):",
    sum(mom_df$significant), "of", nrow(mom_df), "families\n")
 
# =====================================================================
# Save
# =====================================================================
saveRDS(pup_df, file.path(out_dir,
    paste0("family_tests_pup_", data_type, ".rds")))
write.csv(pup_df, file.path(out_dir,
    paste0("family_tests_pup_", data_type, ".csv")), row.names = FALSE)
 
saveRDS(mom_df, file.path(out_dir,
    paste0("family_tests_mom_", data_type, ".rds")))
write.csv(mom_df, file.path(out_dir,
    paste0("family_tests_mom_", data_type, ".csv")), row.names = FALSE)
 