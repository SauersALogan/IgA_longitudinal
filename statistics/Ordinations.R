# =============================================================================
# Script:       XX_ordination.R
# Author:       Logan A. Sauers
# Repository:   github.com/SauersALogan/<repo-name>
# Description:  Absolute-abundance (spike-in adjusted, both spike-ins)
#               community analysis: Euclidean PCoA ordination with
#               PERMANOVA for Genotype and mom_genotype
# Dependencies: abundances.R (creates Data_absolute_Both.rds), plotting.R
# Input:        Data_absolute_Both.rds (from OUT_DIR)
# Output:       A .svg of PCoA coloured by Genotype. [Figure/panel TO FILL]
#               A .svg of PCoA coloured by mom_genotype. [Figure/panel TO FILL]
# =============================================================================
library(phyloseq)
library(vegan)
library(ggplot2)
library(svglite)

# ---- Paths (FIG_DIR / OUT_DIR set in .Renviron) ----
fig_dir <- Sys.getenv("fig_dir")
out_dir <- Sys.getenv("out_dir")
if (fig_dir == "" || out_dir == "") stop("fig_dir or out_dir not set in .Renviron.")

source("statistics/functions/data_structuring.R")
source("statistics/functions/plotting.R")

# ---- Load absolute abundance data ----
Data_rel <- subset_samples(Data, Purpose == "assemblage")

tax <- as.data.frame(tax_table(Data_rel))
spikein_taxa <- rownames(tax)[
    grepl("imtechella|allobacillus", tax$species, ignore.case = TRUE)
]
Data_rel <- prune_taxa(!taxa_names(Data_rel) %in% spikein_taxa, Data_rel)

# =====================================================================
# Distance, PERMANOVA and ordination
# =====================================================================
Data_clr <- microbiome::transform(Data_rel, "clr")
clr.dist <- phyloseq::distance(Data_clr, "euclidean", type = "samples")
ord.data <- as(sample_data(Data_clr), "data.frame")

## ---- Genotype: varies within cage, so permute samples within cages ----
genotype_PERMAnova <- adonis2(clr.dist ~ Genotype * sampling_age,
    data = ord.data, permutations = 9999, strata = ord.data$Cage)
genotype_PERMAnova

## ---- mom_genotype: constant within cage, so permute whole cages ----
clr_mat <- as(otu_table(Data_clr), "matrix")
if (taxa_are_rows(Data_clr)) clr_mat <- t(clr_mat)           # samples x taxa
 
cage_n     <- table(ord.data$Cage)
cage_means <- rowsum(clr_mat, group = ord.data$Cage) / as.vector(cage_n)
 
cage_meta <- unique(ord.data[, c("Cage", "mom_genotype")])
cage_meta <- cage_meta[match(rownames(cage_means), cage_meta$Cage), ]
 
cage.dist <- dist(cage_means, method = "euclidean")
mom_PERMAnova <- adonis2(cage.dist ~ mom_genotype, data = cage_meta,
    permutations = 9999)
mom_PERMAnova

ordination <- ordinate(Data_clr, method = "MDS", distance = clr.dist)

axis_pct <- round(ordination$values$Relative_eig[1:2] * 100, 1)
x_lab <- paste0("PCoA 1 [", axis_pct[1], "%]")
y_lab <- paste0("PCoA 2 [", axis_pct[2], "%]")

# =====================================================================
# Plotting
# =====================================================================
plot_abs_ordination <- function(group_var, colours) {
    p <- plot_ordination(Data_abs, ordination,
        color = group_var, axes = c(1, 2))
    p$layers <- p$layers[-1]

    p +
        geom_point(size = pt_small * 1.5) +
        theme_for_figures(show_labels = TRUE) +
        theme(plot.margin = ggplot2::margin(0, 0, 0, 5, "mm")) +
        scale_color_manual(name = group_var, values = colours) +
        labs(x = x_lab, y = y_lab) +
        scale_x_continuous(expand = expansion(mult = c(0.05, 0.05))) +
        scale_y_continuous(expand = expansion(mult = c(0.05, 0.05)))
}

genotype_ord_plot <- plot_abs_ordination("Genotype", genotype_colours)
genotype_ord_plot

mom_ord_plot <- plot_abs_ordination("mom_genotype", mom_genotype_colours)
mom_ord_plot

# =====================================================================
# Save
# =====================================================================
panel_w <- 75
panel_h <- 75

ggsave(file.path(fig_dir, "ordination_genotype.pdf"), genotype_ord_plot,
    width = panel_w, height = panel_h, units = "mm", dpi = 600)

ggsave(file.path(fig_dir, "ordination_mom_genotype.pdf"), mom_ord_plot,
    width = panel_w, height = panel_h, units = "mm", dpi = 600)