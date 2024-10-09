
# This code is to plot Fig.8 

# WANG. 26-Mar-2024

# Remove all objects created before to prevent clash
rm(list=ls())

library(dplyr)
library(ggseg)
library(ggplot2)
library(tidyr)
library(patchwork)


# Set the working directory to the path where your files are located
setwd("/Users/wang/Desktop/Research_projects/BBSC/Functional/fmri/derivatives/freesurfer-v7.2.0/")
# setwd("//Users/joeywang/Desktop/BBSC/Functional/fmri/derivatives/freesurfer-v7.2.0/")

# generated from 7.1
load("plots_brain_thickness_per_change_qc.RData")
load("plots_brain_thickness_per_change_orig.RData")
load("plots_brain_thickness_cv_qc.RData")
load("plots_brain_thickness_cv_orig.RData")

# generated from 8.1
load("plots_brain_area_per_change_qc.RData")
load("plots_brain_area_per_change_orig.RData")
load("plots_brain_area_cv_qc.RData")
load("plots_brain_area_cv_orig.RData")

# generated from 9.1
load("plots_brain_volume_per_change_qc.RData")
load("plots_brain_volume_per_change_orig.RData")
load("plots_brain_volume_cv_qc.RData")
load("plots_brain_volume_cv_orig.RData")


plots_brain_thickness_cv_orig[[3]] <- plots_brain_thickness_cv_orig[[3]] +
  annotate("text", x = Inf, y = Inf, label = "thickness", vjust = -0.5, hjust = 0.3, family = "Arial", fontface = "bold", size = 8) +
  theme(plot.margin = margin(t = 30, r = 2, b = 10, l = 2, unit = "pt"), legend.position = "none") +
  scale_fill_viridis_c(option = "viridis", direction = 1, limits = c(0, 8), breaks = c(2, 4, 6), labels = c("2", "4", "6"))

plots_brain_thickness_cv_qc[[3]] <- plots_brain_thickness_cv_qc[[3]] + 
  annotate("text", x = Inf, y = Inf, label = "", vjust = 0, hjust = 2.2, family = "Arial", fontface = "bold", size = 8) +
  theme(plot.margin = margin(t = 30, r = 2, b = 10, l = 2, unit = "pt"), legend.position = "none") +
  scale_fill_viridis_c(option = "viridis", direction = 1, limits = c(0, 8), breaks = c(2, 4, 6), labels = c("2", "4", "6"))


plots_brain_area_cv_orig[[3]] <- plots_brain_area_cv_orig[[3]] +
  annotate("text", x = Inf, y = Inf, label = "area", vjust = -0.5, hjust = 0.3, family = "Arial", fontface = "bold", size = 8) +
  theme(plot.margin = margin(t = 30, r = 2, b = 10, l = 2, unit = "pt"), legend.position = "none") +
  scale_fill_viridis_c(option = "viridis", direction = 1, limits = c(0, 8), breaks = c(2, 4, 6), labels = c("2", "4", "6"))

plots_brain_area_cv_qc[[3]] <- plots_brain_area_cv_qc[[3]] + 
  annotate("text", x = Inf, y = Inf, label = "", vjust = 0, hjust = 2.2, family = "Arial", fontface = "bold", size = 8) +
  theme(plot.margin = margin(t = 30, r = 2, b = 10, l = 2, unit = "pt"), legend.position = "none") +
  scale_fill_viridis_c(option = "viridis", direction = 1, limits = c(0, 8), breaks = c(2, 4, 6), labels = c("2", "4", "6"))


plots_brain_volume_cv_orig[[3]] <- plots_brain_volume_cv_orig[[3]] +
  annotate("text", x = Inf, y = Inf, label = "volume", vjust = -0.5, hjust = 0.3, family = "Arial", fontface = "bold", size = 8) +
  theme(plot.margin = margin(t = 30, r = 2, b = 10, l = 2, unit = "pt"), legend.position = "none") +
  scale_fill_viridis_c(option = "viridis", direction = 1, limits = c(0, 8), breaks = c(2, 4, 6), labels = c("2", "4", "6"))


plots_brain_volume_cv_qc[[3]] <- plots_brain_volume_cv_qc[[3]] +
  annotate("text", x = Inf, y = Inf, label = "", vjust = 0, hjust = 2.2, family = "Arial", fontface = "bold", size = 8) +
  theme(plot.margin = margin(t = 30, r = 2, b = 10, l = 2, unit = "pt"),
        legend.title = element_blank(),
        legend.text = element_text(size = 12,family = "Arial", face = "bold", colour = "black")) +
  scale_fill_viridis_c(option = "viridis", direction = 1, limits = c(0, 8), breaks = c(2, 4, 6), labels = c("2", "4", "6"))


plots_brain_thickness_per_change_orig[[3]] <- plots_brain_thickness_per_change_orig[[3]] +
  theme(axis.text.y = element_text(size = 18, family = "Arial",
                                   face = "bold", colour = "black"),
        plot.margin = margin(t = 30, r = 2, b = 10, l = 2, unit = "pt"),
        legend.position = "none")

plots_brain_thickness_per_change_qc[[3]] <- plots_brain_thickness_per_change_qc[[3]] + 
  theme(axis.text.y = element_blank(), 
        plot.margin = margin(t = 30, r = 20, b = 10, l = 2, unit = "pt"),
        legend.position = "none")

plots_brain_area_per_change_orig[[3]] <- plots_brain_area_per_change_orig[[3]] + 
  theme(axis.text.y = element_blank(), 
        plot.margin = margin(t = 30, r = 0, b = 10, l = 10, unit = "pt"),
        legend.position = "none")

plots_brain_area_per_change_qc[[3]] <- plots_brain_area_per_change_qc[[3]] + 
  theme(axis.text.y = element_blank(), 
        plot.margin = margin(t = 30, r = 20, b = 10, l = 0, unit = "pt"),
        legend.position = "none")

plots_brain_volume_per_change_orig[[3]] <- plots_brain_volume_per_change_orig[[3]] + 
  theme(axis.text.y = element_blank(), 
        plot.margin = margin(t = 30, r = 2, b = 10, l = 10, unit = "pt"),
        legend.position = "none")

plots_brain_volume_per_change_qc[[3]] <- plots_brain_volume_per_change_qc[[3]] + 
  theme(axis.text.y = element_blank(), 
        legend.title = element_blank(), 
        legend.text = element_text(size = 16, family = "Arial", face = "bold", colour = "black"), 
        plot.margin = margin(t = 30, r = 2, b = 10, l = 2, unit = "pt"))


Fig8 <- plots_brain_thickness_cv_orig[[3]] + plots_brain_thickness_cv_qc[[3]] + 
  plots_brain_area_cv_orig[[3]] + plots_brain_area_cv_qc[[3]] +
  plots_brain_volume_cv_orig[[3]] + plots_brain_volume_cv_qc[[3]]+
  plots_brain_thickness_per_change_orig[[3]] + plots_brain_thickness_per_change_qc[[3]] +
  plots_brain_area_per_change_orig[[3]] + plots_brain_area_per_change_qc[[3]] +
  plots_brain_volume_per_change_orig[[3]] + plots_brain_volume_per_change_qc[[3]] + 
  plot_layout(ncol = 6)


ggsave("/Users/wang/Library/CloudStorage/OneDrive-UniversityofBergen/Desktop/UiB/Manuscripts/BBSC/4.T1w_Image/Figures/Fig8.png",
       plot = Fig8, width = 16, height = 10, units = 'in', dpi = 300)

