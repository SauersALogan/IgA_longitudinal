# =============================================================================
# Script:       data_structuring.R
# Author:       Logan A. Sauers
# Repository:   github.com/SauersALogan/<repo-name>
# Description:  Loads MetaPhlAn abundances, mouse metadata and SGB taxonomy,
#               coerces metadata types, builds the phyloseq object, and
#               extracts spike-in abundances.
# Input:        First_batch_metaphlan_11-11_2026.csv
#               Mouse_metadata_11-11-2026.csv
#               SGB_taxonomy_11-11-2026.csv
# Output:       Data (phyloseq object)
#               metadata
#               spikein_abundance (samples x spike-in taxa)
# =============================================================================
library(phyloseq)
library(dplyr)

## ---- Global options (inherited by every script that sources this) ----
options(contrasts = c("contr.sum", "contr.poly"))
filter <- dplyr::filter

## ---- Paths (DATA_DIR is set in .Renviron) ----
data_dir <- Sys.getenv("data_dir")
if (data_dir == "") stop("data_dir is not set. Add it to .Renviron and restart R.")

## ---- Load data ----
counts   <- read.csv(file.path(data_dir, "All_batches_metaphlan_18-11-2026.csv"),
    row.names = 1)
metadata <- read.csv(file.path(data_dir, "Mouse_metadata_11-11-2026.csv"),
    row.names = 1)
Taxa     <- read.csv(file.path(data_dir, "SGB_taxonomy_11-11-2026.csv"),
    row.names = 1)

## ---- Factor coercion ----
metadata$Cage             <- as.factor(metadata$Cage)
metadata$Gender           <- as.factor(metadata$Gender)
metadata$Sample_timepoint <- as.factor(metadata$Sample_timepoint)
metadata$Genotype         <- as.factor(metadata$Genotype)
metadata$mom_genotype     <- as.factor(metadata$mom_genotype)
metadata$Batch            <- as.factor(metadata$Batch)
metadata$sampling_age     <- as.factor(metadata$sampling_age)
metadata$abx_timepoint    <- as.factor(metadata$abx_timepoint)
metadata$Mill_run <- as.factor(metadata$Mill_run)

## ---- Phyloseq object ----
OTU <- otu_table(counts, taxa_are_rows = TRUE)
TAX <- tax_table(as.matrix(Taxa))
phyloseq_metadata <- sample_data(metadata)

Data <- phyloseq(OTU, phyloseq_metadata, TAX)

cat("Phyloseq object:", ntaxa(Data), "taxa,", nsamples(Data), "samples\n")
