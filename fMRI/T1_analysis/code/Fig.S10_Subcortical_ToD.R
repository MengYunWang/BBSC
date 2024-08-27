# This code is to plot the time change of subcortex Fig.S10

# WANG. 04-April-2024

# Remove all objects created before to prevent clash
rm(list = ls())

library(dplyr)
library(ggseg)
library(ggplot2)
library(tidyr)
library(patchwork)
library(readr)
library(reshape2)
library(stringr)

# Set the working directory to the path where your files are located
setwd("/Users/wang/Desktop/Research_projects/BBSC/Functional/fmri/derivatives/freesurfer-v7.2.0/")
# setwd("//Users/joeywang/Desktop/BBSC/Functional/fmri/derivatives/freesurfer-v7.2.0/")

#### ------------------------------------------------Part I Load the data

# generated from 10.2
load("percentage_change_subcortex_long.RData")

##-----------------------------------------------

brain_regions <- c( "Thalamus", "Putamen", "Pallidum", "Hippocampus",
                   "Caudate", "Amygdala", "Accumbens-area")

##---generate the subcortex plots
generate_plots_for_subject <- function(subject_id) {
  lapply(brain_regions, function(brain_region) {
    data2plot <- percentage_change_long %>%
      filter(subject == subject_id) %>%
      filter(region == unlist(brain_region))
    
    plot <- ggplot(data2plot) +
      geom_violin(
        aes(
          x = day,
          y = subcortex_per_change,
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
          y = subcortex_per_change,
          group = day,
          color = day
        ),
        binaxis = 'y',
        stackdir = 'down',
        position = position_dodge(0.8),
        dotsize = 0.5
      ) +
      scale_color_manual(values = c(
        "morning" = "#1f77b4",
        "afternoon" = "#ff7f0e"
      )) +
      scale_fill_manual(values = c(
        "morning" = "#1f77b4",
        "afternoon" = "#ff7f0e"
      )) +
      theme_minimal() +
      theme(
        plot.margin = margin(
          t = 5,
          r = 5,
          b = 0,
          l = 20,
          unit = "pt"
        ),
        axis.title.x = element_blank(),
        axis.title.y = element_text(
          size = 16,
          family = "Arial",
          face = "bold",
          colour = "black"
        ),
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
           y = brain_region)
    
    # Add x-axis tile based on the value of brain
    if (brain_region == "Accumbens-area") {
      plot <- plot +  
        ylim(-20, 20) +
        theme(
          axis.text.x = element_text(
            size = 14,
            family = "Arial",
            face = "bold",
            colour = "black",
            angle = 45,
            hjust = 1,
            vjust = 1
          )
        )
    } else {
      plot <- plot + 
        ylim(-8, 8) +
        theme(axis.text.x = element_blank())
    }
    
    
    return(plot)
  })
}

plots_time_subcortex_sub1 <- generate_plots_for_subject("sub1")
plots_time_subcortex_sub2 <- generate_plots_for_subject("sub2")
plots_time_subcortex_sub3 <- generate_plots_for_subject("sub3")

# Function to combine and arrange plots for a given range of indices
combine_plots_for_range <- function(start_index, end_index, plots_list1, plots_list2, plots_list3) {
  combined_plots <- NULL
  
  for (i in start_index:end_index) {
    # Combine the current plot from each subject into a single row
    current_row <- plots_list1[[i]] + 
      plots_list2[[i]] + 
      plots_list3[[i]] +
      plot_layout(ncol = 3, guides = 'collect', axis_titles = "collect")
    
    # Add the current row to the combined plots
    if (is.null(combined_plots)) {
      combined_plots <- current_row
    } else {
      combined_plots <- (combined_plots / current_row) +
        plot_layout(guides = 'collect', axis_titles = "collect")
    }
  }
  
  # After combining all plots, set the global layout
  return(combined_plots + plot_layout(ncol = 1, guides = 'collect'))
}

# Create three segments
subcortex_half1 <- combine_plots_for_range(1, 7, plots_time_subcortex_sub1, plots_time_subcortex_sub2, plots_time_subcortex_sub3)
# subcortex_half2 <- combine_plots_for_range(5, 7, plots_time_subcortex_sub1, plots_time_subcortex_sub2, plots_time_subcortex_sub3)

FigS10 <- (subcortex_half1) + plot_layout(guides = 'collect')

ggsave("/Users/wang/Library/CloudStorage/OneDrive-UniversityofBergen/Desktop/UiB/Manuscripts/BBSC/4.T1w_Image/Figures/FigS10.tiff",
       plot = FigS10, width = 13, height = 18, units = 'in', dpi = 300, compression = "lzw")

# ggsave("/Users/joeywang/Library/CloudStorage/OneDrive-UniversityofBergen/Desktop/UiB/Manuscripts/BBSC/4.T1w_Image/Figures/FigS10.tiff",
#        plot = FigS10, width = 14, height = 20, units = 'in', dpi = 300, compression = "lzw")

