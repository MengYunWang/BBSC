
# This code is to plot the Fig.4
# a. distribution of the IQM
# b. correlations between IQMs
# c. the percentage change of the mean value of the different pipelines

# WANG. 27-Mar-2024



# Remove all objects created before to prevent clash
rm(list = ls())

library(dplyr)
library(ggseg)
library(ggplot2)
library(tidyr)
library(patchwork)
library(readr)
library(reshape2)

# Set the working directory to the path where your files are located
# setwd("/Users/wang/Desktop/Research_projects/BBSC/Functional/fmri/derivatives/freesurfer-v7.2.0/")
setwd("//Users/joeywang/Desktop/BBSC/Functional/fmri/derivatives/freesurfer-v7.2.0/")

#---------------------------plot the left side of figure, Fig.4AB
load("plots_brain_area_mean_qc.RData") # from 8.1
load("plots_area_ordered.RData") # from 8.2

plots_brain_mean_qc[[1]] <- plots_brain_mean_qc[[1]] + 
  theme(plot.margin = margin(t = 20, r = 5, b = 0, l = 5, unit = "pt"), 
        legend.position = "none") +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = "Sub1",
    angle = 90,
    vjust = -32,
    hjust = 1.2,
    family = "Arial",
    fontface = "bold",
    size = 6
  )

plots_brain_mean_qc[[2]] <- plots_brain_mean_qc[[2]] +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(
      size = 12,
      family = "Arial",
      face = "bold",
      colour = "black"
    )
  ) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = "Sub2",
    angle = 90,
    vjust = -32,
    hjust = 1.2,
    family = "Arial",
    fontface = "bold",
    size = 6
  )

plots_brain_mean_qc[[3]] <-
  plots_brain_mean_qc[[3]] + theme(legend.position = "none") +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = "Sub3",
    angle = 90,
    vjust = -32,
    hjust = 1.2,
    family = "Arial",
    fontface = "bold",
    size = 6
  )

plots_area_mean_ordered <- plots_brain_mean_qc[[1]] + plots_brain_mean_qc[[2]] + plots_brain_mean_qc[[3]] + 
  (plots_area_ordered[[1]] + theme (plot.margin = margin(t = 0, r = 10, b = 10, l = 5, unit = "pt"))) + plot_layout(ncol =1)

#------------------------Plot the area distribution along qc Fig.4C

percentage_change_area_long <- read.csv("percentage_change_area_long.csv", row.names =1)
qc_2_plot <-
  c(
    "cnr",
    "snr_total",
    "wm2max",
    "cjv",
    "efc",
    "inu_med",
    "fber",
    "fwhm_avg",
    "rpve_gm"
  )

subject_colors <-
  c(
    'sub1' = '#D40000',
    'sub2' = '#377EB8',
    'sub3' = '#4DAF4A'
  )

subject_grayscale <- c(
  "sub1" = "#D3D3D3",  # Light gray
  "sub2" = "#A9A9A9", # Medium gray
  "sub3" = "#808080" # Dark gray
)


generate_plot <- function(data, variable_name, subject_colors, use_grayscale = FALSE) {
  p <- ggplot(
    data,
    aes_string(
      x = variable_name,
      y = "area_per_change",
      color = "subject",
      group = "subject"
    )
  ) +
    geom_point(size = 0.5, alpha = 0.5) +
    labs(
      title = variable_name,
      x = "",
      y = ""
    ) +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      plot.title = element_text(
        hjust = 0.5,
        size = 18,
        family = "Arial",
        face = "bold",
        color = "black"
      ),
      axis.ticks.x = element_line(color = "black"),
      axis.text.x = element_text(
        size = 14,
        family = "Arial",
        face = "bold",
        color = "black"
      ),
      axis.line.x = element_line(color = "black"),
      axis.ticks.y = element_blank(),
      axis.text.y = element_text(
        size = 14,
        family = "Arial",
        face = "bold",
        color = "black"
      ),
      axis.line.y = element_line(color = "black"),
      axis.title.y = element_blank()
    ) +
    geom_hline(
      yintercept = c(0),
      linetype = "solid",
      color = "grey"
    ) +
    geom_hline(
      yintercept = c(-10, 10),
      linetype = "dashed",
      color = "green"
    ) +
    geom_hline(
      yintercept = c(-15, 15),
      linetype = "dashed",
      color = "red"
    ) +
    scale_y_continuous(limits = c(-20, 20), breaks = c(-15, 0, 15))
  
  # Choose color scale based on use_grayscale flag
  if(use_grayscale) {
    p <- p + scale_color_manual(values = subject_grayscale)
  } else {
    p <- p + scale_color_manual(values = subject_colors)
  }
  
  return(p)
}

plots_area_qc_list <- lapply(qc_2_plot, function(var) {
  generate_plot(percentage_change_area_long, var, subject_colors, use_grayscale = var %in% c(""))
})


plots_area_qc_list[[1]] <- plots_area_qc_list[[1]] + 
  theme (plot.margin = margin(t = 30, r = 5, b = 0, l = 5, unit = "pt")) +
  scale_x_continuous(breaks = c(1.5, 2.0), labels = c("1.5", "2.0"))

plots_area_qc_list[[2]] <- plots_area_qc_list[[2]] + 
  scale_x_continuous(breaks = c(6.0, 7.0, 8.0), labels = c("6.0", "7.0","8.0"))

plots_area_qc_list[[6]] <- plots_area_qc_list[[6]] + 
  scale_x_continuous(breaks = c(0.4, 0.5), labels = c("0.4", "0.5"))

plots_area_qc_list[[7]] <- plots_area_qc_list[[7]] + 
  scale_x_continuous(trans = 'log10', 
                     breaks = c(300, 3000, 30000), labels = c("300","3000","30000"))

plots_area_qc_list[[8]] <- plots_area_qc_list[[8]] + 
  scale_x_continuous(breaks = c(4.0, 4.2, 4.4), labels = c("4.0", "4.2", "4.4"))

plots_area_qc_list[[9]] <- plots_area_qc_list[[9]] + 
  scale_x_continuous(breaks = c(7.25, 7.75,8.25),labels = c("7.25", "7.75","8.25"))

plots_area_qc <- wrap_plots(plots_area_qc_list, ncol = 3)

#-------------------plot the right bottom, Fig.4DE

percentage_change_area_long <- percentage_change_area_long %>%
  filter(!((time_seq == 1 | time_seq == 7) & subject == "sub3"))

subj_id <- c("sub1", "sub2", "sub3")

plots_area_day <- 
  lapply (subj_id, function (subj) {
    data2plot <- percentage_change_area_long %>%
      filter(subject == subj) %>%
      filter(region == "posteriorcingulate") %>%
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
        plot.margin = margin(t = 20, r = 5, b = 10, l = 5, unit = "pt"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(
          size = 14,
          family = "Arial",
          face = "bold",
          colour = "black",
          angle = 45,
          hjust = 1,
          vjust = 1),
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

#---------------------combine together
design <- "AABBCC
           DDEEFF
           GGHHII
           JJKKLL
           JJKKLL"


Fig4 <- plots_area_mean_ordered | (plots_area_qc + plots_area_day[[1]] + 
                                     plots_area_day[[2]] + plots_area_day[[3]] +
                                     plot_layout(design = design))

# ggsave("/Users/wang/Library/CloudStorage/OneDrive-UniversityofBergen/Desktop/UiB/Manuscripts/BBSC/4.T1w_Image/Figures/Fig4.tiff", 
#        plot = Fig4, width = 15.5, height = 11, units = 'in', dpi = 300, compression = "lzw")
# 
ggsave("/Users/joeywang/Library/CloudStorage/OneDrive-UniversityofBergen/Desktop/UiB/Manuscripts/BBSC/4.T1w_Image/Figures/Fig4.tiff",
       plot = Fig4, width = 15.5, height = 11, units = 'in', dpi = 300, compression = "lzw")



