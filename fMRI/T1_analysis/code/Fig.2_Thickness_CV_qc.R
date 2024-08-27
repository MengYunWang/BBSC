
# This code is to plot fig.2

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
load("plots_brain_thickness_cv_qc.RData")

plots_brain_thickness_cv_qc[[1]] <- plots_brain_thickness_cv_qc[[1]] +
  annotate("text", x = Inf, y = Inf, label = "Sub1", vjust = 0, hjust = 2.2, family = "Arial", fontface = "bold", size = 8) +
  theme(plot.margin = margin(t = 30, r = 2, b = 10, l = 2, unit = "pt"), legend.position = "none") +
  scale_fill_viridis_c(option = "viridis", direction = 1, limits = c(0, 8), breaks = c(2, 4, 6), labels = c("2", "4", "6"))


plots_brain_thickness_cv_qc[[2]] <- plots_brain_thickness_cv_qc[[2]] + 
  annotate("text", x = Inf, y = Inf, label = "Sub2", vjust = 0, hjust = 2.2, family = "Arial", fontface = "bold", size = 8) +
  theme(plot.margin = margin(t = 30, r = 2, b = 10, l = 2, unit = "pt"), legend.position = "none") +
  scale_fill_viridis_c(option = "viridis", direction = 1, limits = c(0, 8), breaks = c(2, 4, 6), labels = c("2", "4", "6"))


plots_brain_thickness_cv_qc[[3]] <- plots_brain_thickness_cv_qc[[3]] +
  annotate("text", x = Inf, y = Inf, label = "Sub3", vjust = 0, hjust = 2.2, family = "Arial", fontface = "bold", size = 8) +
  theme(plot.margin = margin(t = 30, r = 2, b = 10, l = 2, unit = "pt"),
        legend.title = element_blank(),
        legend.text = element_text(size = 12,family = "Arial", face = "bold", colour = "black")) +
  scale_fill_viridis_c(option = "viridis", direction = 1, limits = c(0, 8), breaks = c(2, 4, 6), labels = c("2", "4", "6"))


plots_brain_thickness_per_change_qc[[1]] <- plots_brain_thickness_per_change_qc[[1]] + theme(axis.text.y = element_text(size = 18, family = "Arial",
                                                                              face = "bold", colour = "black"),
                                                   plot.margin = margin(t = 30, r = 2, b = 10, l = 2, unit = "pt"),
                                                   legend.position = "none")

plots_brain_thickness_per_change_qc[[2]] <- plots_brain_thickness_per_change_qc[[2]] + theme(axis.text.y = element_blank(), 
                                                   plot.margin = margin(t = 30, r = 2, b = 10, l = 2, unit = "pt"),
                                                   legend.position = "none")


plots_brain_thickness_per_change_qc[[3]] <- plots_brain_thickness_per_change_qc[[3]] + theme(axis.text.y = element_blank(), 
                                                         legend.title = element_blank(), 
                                                         legend.text = element_text(size = 16, family = "Arial", face = "bold", colour = "black"), 
                                                         plot.margin = margin(t = 30, r = 2, b = 10, l = 2, unit = "pt")
                                                         )


Fig2 <- plots_brain_thickness_cv_qc[[1]] + plots_brain_thickness_cv_qc[[2]] + plots_brain_thickness_cv_qc[[3]] +
  plots_brain_thickness_per_change_qc[[1]] + plots_brain_thickness_per_change_qc[[2]] + plots_brain_thickness_per_change_qc[[3]] + plot_layout(ncol = 3)


ggsave("/Users/wang/Library/CloudStorage/OneDrive-UniversityofBergen/Desktop/UiB/Manuscripts/BBSC/4.T1w_Image/Figures/Fig2test.tiff",
       plot = Fig2, width = 13, height = 12, units = 'in', dpi = 300, compression = "lzw")

# ggsave("/Users/joeywang/Library/CloudStorage/OneDrive-UniversityofBergen/Desktop/UiB/Manuscripts/BBSC/4.T1w_Image/Figures/Fig3.tiff", 
#        plot = Fig3, width = 16, height = 12, units = 'in', dpi = 300, compression = "lzw")

######## ----------------------------- ----------------------another way to plot
# rm(list=ls())
# load("brain_thickness_per_change_qc.RData")
# load("brain_thickness_per_change_orig.RData")
# load("brain_thickness_cv_qc.RData")
# load("brain_thickness_cv_orig.RData")
# 
# ##### ------------------do some ajustment for each plot (Fig. 4)
# # Only keep the first y-axis tick lables
# plots_brain_cv[[1]] <- plots_brain_cv[[1]] +
#   annotate("text", x = Inf, y = Inf, label = "Sub1", vjust = 0, hjust = 2.2, family = "Arial", fontface = "bold", size = 8) +
#   theme(plot.margin = margin(t = 30, r = 2, b = 10, l = 2, unit = "pt"), legend.position = "none")
# 
# plots_brain_cv[[2]] <- plots_brain_cv[[2]] + 
#   annotate("text", x = Inf, y = Inf, label = "Sub2", vjust = 0, hjust = 2.2, family = "Arial", fontface = "bold", size = 8) +
#   theme(plot.margin = margin(t = 30, r = 2, b = 10, l = 2, unit = "pt"), legend.position = "none")
# 
# plots_brain_cv[[3]] <- plots_brain_cv[[3]] + 
#   annotate("text", x = Inf, y = Inf, label = "Sub3", vjust = 0, hjust = 2.2, family = "Arial", fontface = "bold", size = 8) +
#   theme(plot.margin = margin(t = 30, r = 2, b = 10, l = 2, unit = "pt"),
#         legend.title = element_blank(),
#         legend.text = element_text(size = 12,family = "Arial", face = "bold", colour = "black"))
# 
# 
# plots_list_box[[1]] <- plots_list_box[[1]] + theme(axis.text.y = element_text(size = 18, family = "Arial",
#                                                                               face = "bold", colour = "black"), 
#                                                    legend.position = "none")
# 
# plots_list_box[[2]] <- plots_list_box[[2]] + theme(axis.text.y = element_blank(), legend.position = "none")
# 
# plots_list_box[[3]] <- plots_list_box[[3]] + theme(axis.text.y = element_blank(),
#                                                    legend.title = element_blank(),
#                                                    legend.text = element_text(size = 16, family = "Arial",
#                                                                               face = "bold", colour = "black"))
# 
# # 
# # plots_brain_mean[[1]] <- plots_brain_mean[[1]] + theme( legend.position = "none") + 
# #   annotate("text", x = Inf, y = Inf, label = "Left", vjust = -1, hjust = 12, family = "Arial", fontface = "bold", size = 8) +
# #   annotate("text", x = Inf, y = Inf, label = "Right", vjust = -1, hjust = 3.5, family = "Arial", fontface = "bold", size = 8) +
# #   annotate("text",  x = Inf, y = Inf, label = "Sub1", angle = 90, vjust = -39, hjust = 1.5,  family = "Arial", fontface = "bold", size = 8)
# # 
# # plots_brain_mean[[2]] <- plots_brain_mean[[2]] + theme(legend.title = element_blank(),
# #                                                        legend.text = element_text(size = 12, family = "Arial", face = "bold", colour = "black"))+
# #   annotate("text",  x = Inf, y = Inf, label = "Sub2", angle = 90, vjust = -39, hjust = 1.5,  family = "Arial", fontface = "bold", size = 8)
# # 
# # plots_brain_mean[[3]] <- plots_brain_mean[[3]] + theme(legend.position = "none")+
# #   annotate("text",  x = Inf, y = Inf, label = "Sub3", angle = 90, vjust = -39, hjust = 1.5,  family = "Arial", fontface = "bold", size = 8)
# 
# 
# plots_brain_cv[[1]] + plots_brain_cv[[2]] + plots_brain_cv[[3]] +
#   plots_list_box[[1]] + plots_list_box[[2]] + plots_list_box[[3]]
# 
# 
# ## for fig. 6
# rm(list=ls())
# 
# load("brain_thickness_per_change_qc.RData")
# load("brain_thickness_per_change_orig.RData")
# load("brain_thickness_cv_qc.RData")
# load("brain_thickness_cv_orig.RData")
# 
# 
# plots_brain_cv[[3]] <- plots_brain_cv[[3]] +
#   scale_fill_viridis_c(option = "viridis", direction = 1, limits = c(0, 6), breaks = c(0, 3, 6), labels = c("0", "3", "6")) +
#   theme(plot.margin = margin(t = 30, r = 2, b = 10, l = 2, unit = "pt"), legend.position = "none")
# 
# plots_brain_cv_qc[[3]] <- plots_brain_cv_qc[[3]] +
#   scale_fill_viridis_c(option = "viridis", direction = 1, limits = c(0, 6), breaks = c(0, 3, 6), labels = c("0", "3", "6")) +
#   theme(plot.margin = margin(t = 30, r = 2, b = 10, l = 2, unit = "pt"),
#         legend.title = element_blank(),
#         legend.text = element_text(size = 12,family = "Arial", face = "bold", colour = "black"))
# 
# plots_list_box[[3]] <- plots_list_box[[3]] + theme(axis.text.y = element_text(size = 18, family = "Arial",
#                                                                               face = "bold", colour = "black"),
#                                                    legend.position = "none")
# 
# plots_list_box_qc[[3]] <- plots_list_box_qc[[3]] + theme(axis.text.y = element_blank(),
#                                                          legend.title = element_blank(),
#                                                          legend.text = element_text(size = 16, family = "Arial",
#                                                                                     face = "bold", colour = "black"))
# 
# plots_brain_cv[[3]] + plots_brain_cv_qc[[3]] +
#   plots_list_box[[3]] + plots_list_box_qc[[3]]






