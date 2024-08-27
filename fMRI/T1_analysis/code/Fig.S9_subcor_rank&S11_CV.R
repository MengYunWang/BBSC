# This code is to ploto the Fig S9 and Fig. S11


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

#### --------------------------------Part I Load the data

# The followings are all generated from 10.1
load("sub_cv_volume_sub_qc.RData")
load("sub_percent_long_qc.RData")

brain_aseg <- as_tibble(aseg)

#### --------------------------------Part II Plot the data
plots_brain_cv <- lapply(sub_cv_volume_sub_qc, function (df) {
  ggplot() +
    geom_brain(atlas = aseg,
               #position = position_brain(hemi ~ side),
               aes(fill = df[brain_aseg$label])) +
    scale_fill_viridis_c(
      option = "magma",
      direction = 1,
      limits = c(0, 10),
      breaks = c(2, 5, 8),
      labels = c("2", "5", "8")
    ) +
    theme_void()
})

plots_list_box <- lapply(names(sub_percent_long_qc), function(name) {
  data <- sub_percent_long_qc[[name]] %>%
    mutate(
      region_order = case_when(
        hemi == "middle" ~ 1,
        # Assign the highest priority to 'middle'
        TRUE ~ as.numeric(factor(hemi, levels = c("lh", "rh"))) + 1  # Ensure 'lh' and 'rh' are ordered after 'middle'
      )) %>%
    arrange(region_order, region) %>%
    mutate(region = factor(region, levels = unique(region))) %>% # Update region to be a factor based on new order
    filter(region != "Total")
  
  # Now, plot with 'region' ordered as desired
  plot <-
    ggplot(data, aes(
      x = percent_change, y = region, fill = hemi
    )) +
    geom_boxplot(position = position_dodge(width = 0.75)) +
    scale_x_continuous(limits = c(-30, 30), breaks = c(-5, 0, 5)) +
    labs(x = "", y = "") +
    theme_minimal() +
    theme(
      strip.background = element_blank(),
      strip.text = element_blank(),
      axis.text.x = element_text(
        size = 18,
        family = "Arial",
        face = "bold",
        colour = "black"
      ),
      axis.text.y = element_text(
        size = 18,
        family = "Arial",
        face = "bold",
        colour = "black"
      ),
      axis.ticks.length = unit(0.2, "cm"),
      axis.ticks.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    geom_vline(
      xintercept = c(-5, 5),
      linetype = "dashed",
      color = "green"
    ) +
    geom_vline(
      xintercept = c(-15, 15),
      linetype = "dashed",
      color = "red"
    )
  return(plot)
})

##### ------------------do some ajustment for each plots
# Only keep the first y-axis tick lables
plots_brain_cv[[1]] <- plots_brain_cv[[1]] +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = "Sub1",
    vjust = 0,
    hjust = 2.2,
    family = "Arial",
    fontface = "bold",
    size = 8
  ) +
  theme(
    plot.margin = margin(
      t = 30,
      r = 2,
      b = 10,
      l = 2,
      unit = "pt"
    ),
    legend.position = "none"
  )

plots_brain_cv[[2]] <- plots_brain_cv[[2]] +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = "Sub2",
    vjust = 0,
    hjust = 2.2,
    family = "Arial",
    fontface = "bold",
    size = 8
  ) +
  theme(
    plot.margin = margin(
      t = 30,
      r = 2,
      b = 10,
      l = 2,
      unit = "pt"
    ),
    legend.position = "none"
  )

plots_brain_cv[[3]] <- plots_brain_cv[[3]] +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = "Sub3",
    vjust = 0,
    hjust = 2.2,
    family = "Arial",
    fontface = "bold",
    size = 8
  ) +
  theme(
    plot.margin = margin(
      t = 30,
      r = 2,
      b = 10,
      l = 2,
      unit = "pt"
    ),
    legend.title = element_blank(),
    legend.text = element_text(
      size = 12,
      family = "Arial",
      face = "bold",
      colour = "black"
    )
  )


plots_list_box[[1]] <-
  plots_list_box[[1]] + theme(
    axis.text.y = element_text(
      size = 18,
      family = "Arial",
      face = "bold",
      colour = "black"
    ),
    legend.position = "none"
  )

plots_list_box[[2]] <-
  plots_list_box[[2]] + theme(axis.text.y = element_blank(), legend.position = "none")

plots_list_box[[3]] <-
  plots_list_box[[3]] + theme(
    axis.text.y = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(
      size = 16,
      family = "Arial",
      face = "bold",
      colour = "black"
    )
  )

FigS11 <-plots_brain_cv[[1]] + plots_brain_cv[[2]] + plots_brain_cv[[3]] +
  plots_list_box[[1]] + plots_list_box[[2]] + plots_list_box[[3]]

ggsave("/Users/wang/Library/CloudStorage/OneDrive-UniversityofBergen/Desktop/UiB/Manuscripts/BBSC/4.T1w_Image/Figures/FigS11.tiff",
       plot = FigS11, width = 13, height = 9, units = 'in', dpi = 300, compression = "lzw")

# ggsave("/Users/joeywang/Library/CloudStorage/OneDrive-UniversityofBergen/Desktop/UiB/Manuscripts/BBSC/4.T1w_Image/Figures/FigS11.tiff",
#        plot = FigS11, width = 13, height = 9, units = 'in', dpi = 300, compression = "lzw")




#-------------------------------------Part III plot the data Fig. S11

load("sub_mean_volume_sub_qc.RData") #from 10.1
load("sub_data_sets_qc.RData") #generated from 10.1

plots_brain_mean <- lapply(sub_mean_volume_sub_qc, function (df) {
  coronal_brain_aseg <- brain_aseg %>%
    filter(side == "coronal", !grepl("\\d", label))
  
  ggplot() +
    geom_brain(atlas = aseg,side = "coronal",
               #position = position_brain(hemi ~ side),
               aes(fill = df[coronal_brain_aseg$label])) +
    scale_fill_viridis_c(
      option = "magma",
      direction = 1,
      limits = c(0, 12000),
      breaks = c(1000, 6000, 11000),
      labels = c("1000", "6000", "11000")
    ) +
    theme_void()
})

plots_volume_ordered  <- lapply(names(data_sets_qc), function(name) {
  df <- data_sets_qc[[name]] %>%
    mutate(subject = name)
  
  data <-
    pivot_longer(
      df,
      cols = !c(subject),
      names_to = "region",
      values_to = "volume"
    ) %>%
    # Create the 'hemi' column based on the 'region' column
    mutate(
      hemi = case_when(
        str_detect(region, "Left") ~ "lh",
        str_detect(region, "Right") ~ "rh",
        TRUE ~ "middle"
      ),
      region = str_replace(region, "Left-", ""),
      region = str_replace(region, "Right-", ""),
      region = str_replace(region, "-Proper", "")
    ) %>%
    filter(hemi != "middle",
           region != "Cerebellum-Cortex",
           region != "Cerebellum-White-Matter",
           region != "Total")
  
  
  average_volume_region <- data %>%
    group_by(region) %>%
    summarise(avg_volume = mean(volume, na.rm = TRUE)) %>%
    arrange(avg_volume) %>%
    mutate(region = factor(region, levels = region))
  
  # update the region
  data_updated <- data %>%
    mutate(region = factor(region, levels = average_volume_region$region))
  
  #
  average_volume <- data_updated %>%
    group_by(region, hemi) %>%
    summarise(avg_volume = mean(volume, na.rm = TRUE)) %>%
    ungroup() %>%
    arrange(region, hemi)
  
  # unique_regions <- unique(average_volume$region)
  # desired_breaks <-
  #   c(head(unique_regions, 3), tail(unique_regions, 3))
  
  plot <-
    ggplot(data_updated, aes(x = region, y = volume, color = hemi)) +
    geom_point(alpha = 0.3) +
    geom_line(data = average_volume,
              aes(x = region, y = avg_volume, group = hemi),
              size = 1) +
    geom_point(data = average_volume,
               aes(x = region, y = avg_volume, group = hemi),
               size = 3) +
    labs(title = "", x = "", y = "") +
    scale_color_manual(values = c("#f66a66", "#14b7bb"),
                       labels = c("lh", "rh")) +
    # scale_x_discrete(breaks = desired_breaks) +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      #
      panel.grid.minor = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(
        size = 16,
        family = "Arial",
        face = "bold",
        colour = "black"
      ),
      axis.ticks.x = element_line(color = "black"),
      #
      axis.text.x = element_text(
        size = 18,
        family = "Arial",
        face = "bold",
        colour = "black",
        angle = 90,
        hjust = 1,
        vjust = 0.5
      ),
      #
      axis.line.x = element_line(color = "black"),
      #
      axis.line.y = element_blank(),
      #
      axis.title.y = element_blank(),
      #
      axis.ticks.y = element_blank(),
      #
      axis.text.y = element_blank() #
    )
  # ylim(700, 26800)
  return(plot)
})

plots_volume_ordered[[1]] <-
  plots_volume_ordered[[1]] + theme(
    legend.position = "none",
    plot.margin = margin(
      t = 30,
      r = 2,
      b = 10,
      l = 2,
      unit = "pt"
    ),
  )
plots_volume_ordered[[2]] <-
  plots_volume_ordered[[2]] + theme(legend.position = "none",plot.margin = margin(
    t = 30,
    r = 2,
    b = 10,
    l = 2,
    unit = "pt"
  ))

plots_volume_ordered[[3]] <-
  plots_volume_ordered[[3]] + theme(plot.margin = margin(
    t = 30,
    r = 2,
    b = 10,
    l = 2,
    unit = "pt"
  ))

plots_brain_mean[[1]] <-
  plots_brain_mean[[1]] + theme(legend.position = "none",plot.margin = margin(
    t = 30,
    r = 2,
    b = 10,
    l = 2,
    unit = "pt"
  )) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = "Sub1",
    vjust = -0.5,
    hjust = 2.2,
    family = "Arial",
    fontface = "bold",
    size = 8
  )

plots_brain_mean[[2]] <- plots_brain_mean[[2]] +
  theme(legend.position = "none",
        plot.margin = margin(
          t = 30,
          r = 2,
          b = 10,
          l = 2,
          unit = "pt"
        )
  ) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = "Sub2",
    vjust = -0.5,
    hjust = 2.2,
    family = "Arial",
    fontface = "bold",
    size = 8
  )

plots_brain_mean[[3]] <-
  plots_brain_mean[[3]] + theme(
    legend.title = element_blank(),
    legend.text = element_text(
      size = 12,
      family = "Arial",
      face = "bold",
      colour = "black"
    ),
    plot.margin = margin(
      t = 30,
      r = 2,
      b = 10,
      l = 2,
      unit = "pt"
    )
  ) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = "Sub3",
    vjust = -0.5,
    hjust = 2.2,
    family = "Arial",
    fontface = "bold",
    size = 8
  )

(plots_brain_mean[[1]] + plots_brain_mean[[2]] + plots_brain_mean[[3]]) / plots_volume_ordered[[1]] 

FigS9 <- (plots_brain_mean[[1]] / plots_volume_ordered[[1]]) |
  (plots_brain_mean[[2]] / plots_volume_ordered[[2]]) |
  (plots_brain_mean[[3]] / plots_volume_ordered[[3]]) 

ggsave("/Users/wang/Library/CloudStorage/OneDrive-UniversityofBergen/Desktop/UiB/Manuscripts/BBSC/4.T1w_Image/Figures/FigS9.tiff",
       plot = FigS9, width = 10, height = 8.5, units = 'in', dpi = 300, compression = "lzw")

# ggsave("/Users/joeywang/Library/CloudStorage/OneDrive-UniversityofBergen/Desktop/UiB/Manuscripts/BBSC/4.T1w_Image/Figures/FigS9.tiff",
#        plot = FigS9, width = 13, height = 9, units = 'in', dpi = 300, compression = "lzw")






