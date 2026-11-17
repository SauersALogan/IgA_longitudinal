# =============================================================================
# Script:       XX_Abundances.R
# Author:       Logan A. Sauers
# Repository:   github.com/SauersALogan/<repo-name>
# Description:  [TO FILL] Spike-in abundance and total 16S abundance analysis
# Dependencies: data_structuring.R, plotting.R
# Input:        Data (from data_structuring.R)
#               [TO FILL] total 16S abundance source
# Output:       [TO FILL]
# =============================================================================
library(future)
library(future.apply)
library(progressr)
library(svglite)
library(glmmTMB)
library(DHARMa)
library(car)
library(tidyr)
library(tibble)

# ---- Paths (FIG_DIR / OUT_DIR set in .Renviron) ----
fig_dir <- Sys.getenv("fig_dir")
out_dir <- Sys.getenv("out_dir")
if (fig_dir == "" || out_dir == "") stop("fig_dir or out_dir not set in .Renviron.")

# ---- Shared objects and plotting setup ----
source("statistics/functions/data_structuring.R")
source("statistics/functions/plotting.R")
source("statistics/functions/models.R")

# =====================================================================
# Parallel and progress setup
# =====================================================================
n_cores <- parallelly::availableCores(omit = 1)
options(future.globals.maxSize = 2000 * 1024^2)
plan(multisession, workers = n_cores)

if (interactive()) handlers(global = TRUE)
handlers("progress")

# =====================================================================
# Spike-in relative abundance comparison
# =====================================================================
Data_assemblage <- subset_samples(Data, Purpose == "assemblage")

## ---- Extract spike-in ----
tax <- as.data.frame(tax_table(Data_assemblage))
spikein_taxa <- rownames(tax)[
    grepl("imtechella|allobacillus", tax$species, ignore.case = TRUE)
]

spikein_df <- psmelt(prune_taxa(spikein_taxa, Data_assemblage))
spikein_df$Spike <- ifelse(grepl("imtechella", spikein_df$species,
    ignore.case = TRUE), "Imtechella", "Allobacillus")
spikein_df

## ---- Model genotype - Imtechella ----
genotype_Imt_model <- glmmTMB(log10(Abundance) ~ Genotype * sampling_age + 
    (1 | Cage), family = gaussian(),
    data = subset(spikein_df, Spike == "Imtechella"))

Anova(genotype_Imt_model, type = "III", component = "cond")

## ---- Genotype diagnostics - Imtechella ----
genotype_Imt_sim <- simulateResiduals(genotype_Imt_model, n = 1000)
plot(genotype_Imt_sim)
testDispersion(genotype_Imt_sim)

## ---- Mom genotype - Imtechella ----
mom_Imt_model <- glmmTMB(log10(Abundance) ~ mom_genotype * sampling_age + 
    (1 | Cage), family = gaussian(),
    data = subset(spikein_df, Spike == "Imtechella"))

Anova(mom_Imt_model, type = "III", component = "cond")

## ---- Mom diagnostics - Imtechella ----
mom_Imt_sim <- simulateResiduals(mom_Imt_model, n = 1000)
plot(mom_Imt_sim)

## ---- Model genotype - Allobacillus ----
genotype_Allo_model <- glmmTMB(log10(Abundance) ~ Genotype * sampling_age + 
    (1 | Cage), family = gaussian(),
    data = subset(spikein_df, Spike == "Allobacillus"))

Anova(genotype_Allo_model, type = "III", component = "cond")

## ---- Genotype diagnostics - Allobacillus ----
genotype_Allo_sim <- simulateResiduals(genotype_Allo_model, n = 1000)
plot(genotype_Allo_sim)
testDispersion(genotype_Allo_sim)

## ---- Mom genotype - Allobacillus ----
mom_Allo_model <- glmmTMB(log10(Abundance) ~ mom_genotype * sampling_age + 
    (1 | Cage), family = gaussian(),
    data = subset(spikein_df, Spike == "Allobacillus"))

Anova(mom_Allo_model, type = "III", component = "cond")

## ---- Mom diagnostics - Allobacillus ----
mom_Allo_sim <- simulateResiduals(mom_Allo_model, n = 1000)
plot(mom_Allo_sim)

## ---- Emmeans ----
genotype_Imt_emm  <- as.data.frame(emmeans(genotype_Imt_model,
    ~ Genotype | sampling_age, type = "response"))
mom_Imt_emm       <- as.data.frame(emmeans(mom_Imt_model,
    ~ mom_genotype | sampling_age, type = "response"))
genotype_Allo_emm <- as.data.frame(emmeans(genotype_Allo_model,
    ~ Genotype | sampling_age, type = "response"))
mom_Allo_emm      <- as.data.frame(emmeans(mom_Allo_model,
    ~ mom_genotype | sampling_age, type = "response"))

## ---- Plotting ----
plot_spikein_emm <- function(emm_df, raw_df, group_var, colours, y_lab) {
    ggplot(emm_df, aes(x = sampling_age, y = response,
            colour = .data[[group_var]])) +
        geom_point(data = raw_df, aes(y = Abundance),
            alpha = 0.3, size = pt_small,
            position = position_jitterdodge(jitter.width = 0.13,
                dodge.width = 0.5, seed = 1)) +
        geom_pointrange(aes(ymin = lower.CL, ymax = upper.CL),
            linewidth = lw_main, size = pt_mean / 6, fatten = 2,
            position = position_dodge(width = 0.5)) +
        scale_colour_manual(values = colours) +
        scale_y_log10() +
        theme_for_figures(show_labels = TRUE) +
        labs(x = "Sampling age", y = y_lab)
}

Imt_df <- subset(spikein_df, Spike == "Imtechella")
Allo_df <- subset(spikein_df, Spike == "Allobacillus")

genotype_Imt_plot <- plot_spikein_emm(genotype_Imt_emm, Imt_df, "Genotype",
    genotype_colours, "Imtechella (%)")
mom_Imt_plot <- plot_spikein_emm(mom_Imt_emm, Imt_df, "mom_genotype",
    mom_genotype_colours, "Imtechella (%)")
genotype_Allo_plot <- plot_spikein_emm(genotype_Allo_emm, Allo_df, "Genotype",
    genotype_colours, "Allobacillus (%)")
mom_Allo_plot <- plot_spikein_emm(mom_Allo_emm, Allo_df, "mom_genotype",
    mom_genotype_colours, "Allobacillus (%)")

# =====================================================================
# Absolute abundance comparison
# =====================================================================
## ---- Spike-in scalers ----
# Option A: Imtechella only
Imt_scaling <- spikein_taxa[grepl("imtechella", tax[spikein_taxa, "species"], 
    ignore.case = TRUE)]
Imt_cells  <- 2.5e7
# Option B: Allobacillus only
Allo_scaling <- spikein_taxa[grepl("allobacillus", tax[spikein_taxa, "species"], 
    ignore.case = TRUE)]
Allo_cells  <- 2.5e7
# Option C: Both combined
Both_scaling <- spikein_taxa
Both_cells <- 5e7

## ---- Shared scalers ----
rel_mat <- as(otu_table(Data), "matrix")
feces_g <- setNames(sample_data(Data)$feces_weight,
    sample_names(Data))[colnames(rel_mat)]

## ---- Absolute abundance function ----
calc_absolute <- function(rel_mat, scaling_taxa, spike_cells, feces_g,
    spikein_taxa) {
    spike_rel <- colSums(rel_mat[scaling_taxa, , drop = FALSE])
    if (any(spike_rel == 0)) warning("Zero spike-in in: ",
        paste(names(spike_rel)[spike_rel == 0], collapse = ", "))
    scale_factor <- spike_cells / (spike_rel * feces_g)
    sweep(rel_mat[!rownames(rel_mat) %in% spikein_taxa, , drop = FALSE],
        2, scale_factor, "*")
}

## ---- Calculation ----
abs_Imt  <- calc_absolute(rel_mat, Imt_scaling,  Imt_cells,  feces_g, spikein_taxa)
abs_Allo <- calc_absolute(rel_mat, Allo_scaling, Allo_cells, feces_g, spikein_taxa)
abs_Both <- calc_absolute(rel_mat, Both_scaling, Both_cells, feces_g, spikein_taxa)

## ---- Merge into one dataframe ----
load_comparison <- data.frame(Sample = colnames(rel_mat),
    load_Imt  = colSums(abs_Imt),
    load_Allo = colSums(abs_Allo),
    load_Both = colSums(abs_Both))

## ---- Reformat data and add metadata ----
load_long <- load_comparison %>%
    pivot_longer(-Sample, names_to = "Scaler", names_prefix = "load_",
        values_to = "cells_per_g") %>%
    left_join(as(sample_data(Data), "data.frame") %>%
        rownames_to_column("Sample"), by = "Sample") %>%
    filter(Purpose == "assemblage")

## ---- Fitmodels per scaler ----
# GLM control is adjusted to ensure convergence, default optimizer fails
fit_load <- function(df, term) {
    f <- as.formula(paste0("log10(cells_per_g) ~ ", term,
        " * sampling_age + (1 | Cage)"))
    glmmTMB(f, family = gaussian(), data = df,
        control = glmmTMBControl(optCtrl = list(iter.max = 1e4, eval.max = 1e4)))
}

load_models <- list()
for (s in c("Imt", "Allo", "Both")) {
    d <- subset(load_long, Scaler == s)
    load_models[[paste0("genotype_", s)]] <- fit_load(d, "Genotype")
    load_models[[paste0("mom_", s)]] <- fit_load(d, "mom_genotype")
}

## ---- Diagnostics ----
load_sims <- lapply(load_models, simulateResiduals, n = 1000)
for (n in names(load_sims)) plot(load_sims[[n]], title = n)

## ---- Anova comparison across scalers ----
anova_comparison <- bind_rows(lapply(names(load_models), function(n) {
    a <- as.data.frame(Anova(load_models[[n]], type = "III",
        component = "cond"))
    a$Term  <- rownames(a)
    a$Model <- n
    a
}))
anova_comparison

## ---- Plotting function updated ----
plot_emm <- function(emm_df, raw_df, group_var, raw_y, colours, y_lab) {
    ggplot(emm_df, aes(x = sampling_age, y = response,
            colour = .data[[group_var]])) +
        geom_point(data = raw_df, aes(y = .data[[raw_y]]),
            alpha = 0.3, size = pt_small,
            position = position_jitterdodge(jitter.width = 0.13,
                dodge.width = 0.5, seed = 1)) +
        geom_pointrange(aes(ymin = lower.CL, ymax = upper.CL),
            linewidth = lw_main, size = pt_mean / 6, fatten = 2,
            position = position_dodge(width = 0.5)) +
        scale_colour_manual(values = colours) +
        scale_y_log10() +
        theme_for_figures(show_labels = TRUE) +
        labs(x = "Sampling age", y = y_lab)
}

## ---- Load emmeans and plots ----
load_plots <- list()
for (n in names(load_models)) {
    is_mom  <- grepl("^mom_", n)
    term    <- if (is_mom) "mom_genotype" else "Genotype"
    colours <- if (is_mom) mom_genotype_colours else genotype_colours
    scaler  <- sub("^(genotype|mom)_", "", n)

    emm_df <- as.data.frame(emmeans(load_models[[n]],
        as.formula(paste0("~ ", term, " | sampling_age")),
        type = "response"))

    load_plots[[n]] <- plot_emm(emm_df,
        subset(load_long, Scaler == scaler),
        group_var = term, raw_y = "cells_per_g",
        colours = colours, y_lab = "Bacterial load (cells/g)")
}
load_plots[["genotype_Imt"]]
load_plots[["mom_Imt"]]
load_plots[["genotype_Allo"]]
load_plots[["mom_Allo"]]
load_plots[["genotype_Both"]]
load_plots[["mom_Both"]]

## ---- Save mom load plots ----
panel_w <- 75   # mm
panel_h <- 75   # mm

for (n in c("mom_Imt", "mom_Allo", "mom_Both")) {
    ggsave(file.path(fig_dir, paste0(n, "_load.pdf")), load_plots[[n]],
        width = panel_w, height = panel_h, units = "mm", dpi = 600)
}

## ---- Absolute abundance phyloseq objects ----
make_abs_phyloseq <- function(abs_mat) {
    ps <- prune_taxa(!taxa_names(Data) %in% spikein_taxa, Data)
    otu_table(ps) <- otu_table(abs_mat, taxa_are_rows = TRUE)
    ps
}

abs_phyloseq <- list(Imt  = make_abs_phyloseq(abs_Imt),
                     Allo = make_abs_phyloseq(abs_Allo),
                     Both = make_abs_phyloseq(abs_Both))

for (s in names(abs_phyloseq)) {
    saveRDS(abs_phyloseq[[s]], file.path(out_dir,
        paste0("Data_absolute_", s, ".rds")))
}