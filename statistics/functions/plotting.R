# =============================================================================
# Script:       plotting.R
# Author:       Logan A. Sauers
# Repository:   github.com/SauersALogan/<repo-name>
# Description:  Shared figure parameters, theme, label-space reservations,
#               and colour palette. Definitions only; no data or plots.
# Output:       theme_for_figures(), size/linewidth constants, colours
# =============================================================================
library(ggplot2)

## ---- Plotting base lines ----
base_sz  <- 7
pt_to_mm <- 2.46
pt_small <- 0.9
pt_mean  <- 2.4
lw_main  <- 0.6
lw_sm    <- 0.3
lw_axis  <- 0.6
arr_len  <- 0.15

theme_for_figures <- function(show_labels = FALSE) {
  base <- theme_classic(base_size = base_sz) +
    theme(
      axis.text    = element_text(colour = "black", size = base_sz),
      legend.text  = element_text(size = base_sz),
      legend.title = element_text(size = base_sz),
      panel.grid   = element_blank(),
      plot.title   = element_text(face = "bold.italic", hjust = 0.5, size = base_sz),
      legend.position = "bottom",
      axis.line    = element_line(colour = "black", linewidth = lw_axis),
      axis.ticks   = element_line(colour = "black", linewidth = lw_axis),
      plot.margin  = ggplot2::margin(0, 0, 0, 0, "mm"),
      axis.title.y = element_text(
          margin = ggplot2::margin(t = 0, r = 3, b = 0, l = 0)))
  if (show_labels) {
    base
  } else {
    base + theme(
      axis.title.x = element_blank(),   # blank the SPECIFIC elements,
      axis.title.y = element_blank(),   # not the general axis.title
      axis.text    = element_blank()
    )
  }
}

## ---- Genotype colours ----
col_het <- "#DD5566" # Genotype Het (red)
col_ko  <- "#4393C3"     # Genotype KO (blue)
col_mom_het <- "#9970AB"   # mom_genotype Het (purple, ColorBrewer PRGn)
col_mom_wt  <- "#5AAE61"   # mom_genotype WT (green, ColorBrewer PRGn)

genotype_colours <- c("Het" = col_het, "KO" = col_ko)
mom_genotype_colours <- c("Het" = col_mom_het, "KO" = col_mom_wt)