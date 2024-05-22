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

# delete two sessions with exsessive head movements
percentage_change_area_long <- percentage_change_area_long %>%
  filter(!((time_seq == 1 | time_seq == 7) & subject == "sub3"))

# function to plot
plot_day <- function (data2plot, brain) {
  plot <- ggplot(data2plot) +
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
         y = brain)+ 
    ylim(-15, 15)+
    
    if(brain == "posteriorcingulate" | brain == "lateralorbitofrontal"
       | brain == "bankssts") {
      plot <- plot + 
        theme(axis.text.x = element_text(
          size = 14,
          family = "Arial",
          face = "bold",
          colour = "black",
          angle = 45,
          hjust = 1,
          vjust = 1)
        )
    } else {
      plot <- plot +
        theme(axis.text.x = element_blank())
    }
  return(plot)
}
# ---------------------plot the percentage change along the data collection

brain_regions <- colnames(percentage_change_area[, 1:34]) %>%
  sort(decreasing = TRUE) %>%
  as.list()

##---generate the area plots Fig.S1
generate_plots_for_subject <- function(subject_id) {
  lapply(brain_regions, function(brain_region) {
    data2plot <- percentage_change_area_long %>%
      filter(subject == subject_id) %>%
      filter(region == unlist(brain_region))%>%
      mutate(day = factor(day, levels = c("morning", "afternoon")))
    
    plots <- plot_day (data2plot, brain_region)
    return(plots)
  })
}

plots_time_area_sub1 <- generate_plots_for_subject("sub1")
plots_time_area_sub2 <- generate_plots_for_subject("sub2")
plots_time_area_sub3 <- generate_plots_for_subject("sub3")

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
area_half1 <- combine_plots_for_range(1, 11, plots_time_area_sub1, plots_time_area_sub2, plots_time_area_sub3)
area_half2 <- combine_plots_for_range(12, 22, plots_time_area_sub1, plots_time_area_sub2, plots_time_area_sub3)
area_half3 <- combine_plots_for_range(23, 34, plots_time_area_sub1, plots_time_area_sub2, plots_time_area_sub3)

area_time <- (area_half1 | area_half2 | area_half3) + plot_layout(guides = 'collect')

ggsave("/Users/wang/Library/CloudStorage/OneDrive-UniversityofBergen/Desktop/UiB/Manuscripts/BBSC/4.T1w_Image/Figures/FigS2.tiff",
       plot = area_time, width = 27, height = 40, units = 'in', dpi = 300, compression = "lzw")

# ggsave("/Users/joeywang/Library/CloudStorage/OneDrive-UniversityofBergen/Desktop/UiB/Manuscripts/BBSC/4.T1w_Image/Figures/FigS2.tiff",
#        plot = area_time, width = 27, height = 40, units = 'in', dpi = 300, compression = "lzw")

