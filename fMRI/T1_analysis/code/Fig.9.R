# This code is to analyze the effect of FBER and plot Fig.9 and Fig. S1-3

# WANG. 03-April-2024

# Remove all objects created before to prevent clash
rm(list = ls())

library(dplyr)
library(ggplot2)
library(reshape2)
library(patchwork)
library(RColorBrewer) 

# Set the working directory to the path where your files are located
setwd("/Users/wang/Desktop/Research_projects/BBSC/Functional/fmri/derivatives/freesurfer-v7.2.0/")
# setwd("//Users/joeywang/Desktop/BBSC/Functional/fmri/derivatives/freesurfer-v7.2.0/")

#----------------------------------Load the data

# generated from 8.2
percentage_change_area <- read.csv("percentage_change_area.csv", row.names = 1)
percentage_change_area_long <- read.csv("percentage_change_area_long.csv", row.names = 1)

percentage_change_area_long <- percentage_change_area_long %>%
  filter(!((time_seq == 1 | time_seq == 7) & subject == "sub3"))


subj_id <- c("sub1", "sub2", "sub3")


# plot the percentage of change of day
plots_area_day <- 
  lapply (subj_id, function (subj) {
    data2plot <- percentage_change_area_long %>%
      filter(subject == subj) %>%
      filter(region == "temporalpole") %>%
      mutate(day = factor(day, levels = c("morning", "afternoon")))
    
    ggplot(data2plot) +
      geom_violin(
        aes(
          x = day,
          y = area_per_change,
          fill = NULL,
          group = day,
          color = day
        ),
        trim = FALSE,
        position = position_dodge(0.8),
        draw_quantiles = NULL
      ) +
      geom_dotplot(
        aes(
          x = day,
          y = area_per_change,
          fill = hemi,
          group = day,
          color = day
        ),
        binaxis = 'y',
        stackdir = 'down',
        position = position_dodge(0.8),
        dotsize = 0.5
      ) +
      scale_color_manual(values = c("morning" = "#1f77b4", "afternoon" = "#ff7f0e")) +
      scale_fill_manual(values = c("morning" = "#1f77b4", "afternoon" = "#ff7f0e")) +
      theme_minimal() +
      theme(
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(
          size = 14,
          family = "Arial",
          face = "bold",
          colour = "black"),
        axis.text.y = element_text(
          size = 14,
          family = "Arial",
          face = "bold",
          colour = "black"
        ),
        legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank()
      ) +
      labs(title = "" ,
           x = NULL,
           y = NULL)+ 
      ylim(-25, 25)
  })

# plot the percentage change of day after xcluding sessions
percentage_change_area_long <- percentage_change_area_long %>%
  filter(!((time_seq == 16) & subject == "sub1")) %>%
  filter(!((time_seq == 23 | time_seq == 32) & subject == "sub2"))


plots_area_day_qc <- 
  lapply (subj_id, function (subj) {
    data2plot <- percentage_change_area_long %>%
      filter(subject == subj) %>%
      filter(region == "temporalpole") %>%
      mutate(day = factor(day, levels = c("morning", "afternoon")))
    
    ggplot(data2plot) +
      geom_violin(
        aes(
          x = day,
          y = area_per_change,
          fill = NULL,
          group = day,
          color = day
        ),
        trim = FALSE,
        position = position_dodge(0.8),
        draw_quantiles = NULL
      ) +
      geom_dotplot(
        aes(
          x = day,
          y = area_per_change,
          fill = hemi,
          group = day,
          color = day
        ),
        binaxis = 'y',
        stackdir = 'down',
        position = position_dodge(0.8),
        dotsize = 0.5
      ) +
      scale_color_manual(values = c("morning" = "#1f77b4", "afternoon" = "#ff7f0e")) +
      scale_fill_manual(values = c("morning" = "#1f77b4", "afternoon" = "#ff7f0e")) +
      theme_minimal() +
      theme(
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(
          size = 14,
          family = "Arial",
          face = "bold",
          colour = "black"),
        axis.text.y = element_text(
          size = 14,
          family = "Arial",
          face = "bold",
          colour = "black"
        ),
        legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank()
      ) +
      labs(title = "" ,
           x = NULL,
           y = NULL)+ 
      ylim(-25, 25)
  })

####------------------plot Fig.9

Fig9 <- plots_area_day[[1]] +  plots_area_day[[2]] + theme( plot.margin = margin(l = 100, unit = "pt"))+
  plots_area_day_qc[[1]]  + plots_area_day_qc[[2]] + plot_layout(nrow=2, axes = 'collect_x')


ggsave("/Users/wang/Library/CloudStorage/OneDrive-UniversityofBergen/Desktop/UiB/Manuscripts/BBSC/4.T1w_Image/Figures/Fig9.tiff",
       plot = Fig9, width = 9, height = 7, units = 'in', dpi = 300, compression = "lzw")

# ggsave("/Users/joeywang/Library/CloudStorage/OneDrive-UniversityofBergen/Desktop/UiB/Manuscripts/BBSC/4.T1w_Image/Figures/Fig9.tiff",
#        plot = Fig9, width = 10, height = 6.5, units = 'in', dpi = 300, compression = "lzw")

