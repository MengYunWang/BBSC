# This code is to ploto the Fig 8 and Fig. S5


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
load("sub_cv_volume_sub.RData")
load("sub_percent_long.RData")

#### --------------------------------Part II Plot the data

coronal_brain_aseg <- as_tibble(aseg) %>%
  filter(side == "coronal", !grepl("\\d", label))

plots_brain_cv <- lapply(sub_cv_volume_sub, function (df) {
  ggplot() +
    geom_brain(atlas = aseg, side = "coronal",
               #position = position_brain(hemi ~ side),
               aes(fill = df[coronal_brain_aseg$label])) +
    scale_fill_viridis_c(
      option = "magma",
      direction = 1,
      limits = c(0, 8),
      breaks = c(2, 4, 6),
      labels = c("2", "4", "6")
    ) +
    theme_void()
})

plots_list_box <- lapply(names(sub_percent_long), function(name) {
  data <- sub_percent_long[[name]] %>%
    mutate(
      region_order = case_when(
        hemi == "middle" ~ 1,
        # Assign the highest priority to 'middle'
        TRUE ~ as.numeric(factor(hemi, levels = c("lh", "rh"))) + 1  # Ensure 'lh' and 'rh' are ordered after 'middle'
      )) %>%
    arrange(region_order, region) %>%
    mutate(region = factor(region, levels = unique(region))) %>% # Update region to be a factor based on new order
    filter(region != "Lateral-Ventricle",
           region != "Cerebellum-Cortex",
           region != "Cerebellum-White-Matter",
           region != "VentralDC",
           hemi != "middle")
  
  # Now, plot with 'region' ordered as desired
  plot <-
    ggplot(data, aes(
      x = percent_change, y = region, fill = hemi
    )) +
    geom_boxplot(position = position_dodge(width = 0.75)) +
    scale_x_continuous(limits = c(-23, 23), breaks = c(-5, 0, 5)) +
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

Fig8 <-plots_brain_cv[[1]] + plots_brain_cv[[2]] + plots_brain_cv[[3]] +
  plots_list_box[[1]] + plots_list_box[[2]] + plots_list_box[[3]]

ggsave("/Users/wang/Library/CloudStorage/OneDrive-UniversityofBergen/Desktop/UiB/Manuscripts/BBSC/4.T1w_Image/Figures/Fig8.tiff",
       plot = Fig8, width = 12, height = 10, units = 'in', dpi = 300, compression = "lzw")

# ggsave("/Users/joeywang/Library/CloudStorage/OneDrive-UniversityofBergen/Desktop/UiB/Manuscripts/BBSC/4.T1w_Image/Figures/Fig8.tiff",
#        plot = Fig8, width = 13, height = 9, units = 'in', dpi = 300, compression = "lzw")


