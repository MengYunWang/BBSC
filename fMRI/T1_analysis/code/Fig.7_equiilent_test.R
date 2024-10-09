# This code is to plot Fig. 7


# WANG. 23-Aug-2024

# Remove all objects created before to prevent clash
rm(list = ls())

library(TOSTER)

# Set the working directory to the path where your files are located
setwd("/Users/wang/Desktop/Research_projects/BBSC/Functional/fmri/derivatives/freesurfer-v7.2.0/")

#from 7.2
load(file="meta_model_day_thickness.RData")
load(file="meta_model_time_thickness.RData")

#from 8.2
load(file="meta_model_day_area.RData")
load(file="meta_model_time_area.RData")

#from 9.2
load(file="meta_model_day_volume.RData")
load(file="meta_model_time_volume.RData")

#from 10.2
load(file="meta_model_day_subcortex.RData")
load(file="meta_model_time_subcortex.RData")


TOST_plot <- function(model_day, model_time, plot_title) {
  TOST_result_day <- TOSTmeta(
    ES = summary(model_day)$beta,
    se = summary(model_day)$se,
    low_eqbound_d = -0.5,
    high_eqbound_d = 0.5,
    alpha = 0.05
  )
  
  TOST_result_time <- TOSTmeta(
    ES = summary(model_time)$beta,
    se = summary(model_time)$se,
    low_eqbound_d = -0.5,
    high_eqbound_d = 0.5,
    alpha = 0.05
  )
  
  
  
  subcortical <- data.frame(
    E = c(TOST_result_day$ES, TOST_result_time$ES),
    LL_CI_TOST = c(TOST_result_day$LL_CI_TOST, TOST_result_time$LL_CI_TOST),
    UL_CI_TOST = c(TOST_result_day$UL_CI_TOST, TOST_result_time$UL_CI_TOST),
    LL_CI_ZTEST = c(TOST_result_day$LL_CI_ZTEST, TOST_result_time$LL_CI_ZTEST),
    UL_CI_ZTEST = c(TOST_result_day$UL_CI_ZTEST, TOST_result_time$UL_CI_ZTEST),
    position = c(0.6, 0.4)
  )
  
  
  ggplot() +
    geom_vline(
      aes(xintercept = 0),
      linetype = "dashed",
      color = "grey",
      size = 1.2
    ) +
    geom_vline(aes(xintercept = -0.5), linetype = "dashed", size = 1.2) +
    geom_vline(aes(xintercept = 0.5), linetype = "dashed", size = 1.2) +
    
    geom_point(
      data = subcortical,
      aes(x = E, y = position),
      shape = 15,
      size = 5,
      color = c("#984EA3", "#FF7F00")
    ) +
    
    geom_segment(
      data = subcortical,
      aes(
        x = LL_CI_TOST,
        xend = UL_CI_TOST,
        y = position,
        yend = position
      ),
      linetype = "solid",
      size = 1.5,
      color = c("#984EA3", "#FF7F00")
    ) +
    geom_segment(
      data = subcortical,
      aes(
        x = LL_CI_ZTEST,
        xend = UL_CI_ZTEST,
        y = position,
        yend = position
      ),
      linetype = "solid",
      size = 0.5,
      color = c("#984EA3", "#FF7F00")
    ) +
    
    coord_cartesian(xlim = c(-1, 1), ylim = c(0, 1)) +
    
    labs(x = "Effect size", title = plot_title) +
    
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      axis.line = element_line(colour = 'black', linewidth = 1),
      # add x, y axis
      strip.background = element_blank(),
      strip.text = element_blank(),
      axis.title.x = element_text(
        # specify x-axis title font
        size = 16,
        family = "Arial",
        face = "bold",
        colour = "black"
      ),
      axis.text.x = element_text(
        size = 14,
        family = "Arial",
        face = "bold",
        colour = "black"
      ),
      plot.title = element_text(
        hjust = 0.5,
        size = 18,
        family = "Arial",
        face = "bold",
        colour = "black"
      )
    )
}

plot_thickness <- TOST_plot(meta_model_day_thickness, meta_model_time_thickness, "")
plot_area <- TOST_plot(meta_model_day_area, meta_model_time_area, "")
plot_volume <- TOST_plot(meta_model_day_volume, meta_model_time_volume, "")
plot_subcortical <- TOST_plot(meta_model_day_subcortex, meta_model_time_subcortex, "")


plot_equi <- (plot_thickness + theme(axis.title.x = element_blank(), plot.margin = margin(1.2, 1.2, 1.2, 1.2, "cm")) |  plot_area + theme(axis.title.x = element_blank(), plot.margin = margin(1.2, 1.2, 1.2, 1.2, "cm"))) /
  (plot_volume |  plot_subcortical)+
  plot_layout(axis_titles = "collect_x") +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 22,
                                family = "Arial",
                                face = "bold",
                                colour = "black"))

ggsave("/Users/wang/Library/CloudStorage/OneDrive-UniversityofBergen/Desktop/UiB/Manuscripts/BBSC/4.T1w_Image/Figures/Fig7_equiv.png",
       plot = plot_equi, width = 9, height = 10, units = 'in', dpi = 300)





