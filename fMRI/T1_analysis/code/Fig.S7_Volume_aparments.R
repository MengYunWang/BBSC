# This code is to plot Fig.S7

# WANG. 20-Feb-2024



# Remove all objects created before to prevent clusing
rm(list=ls())

# Load required packages

library(dplyr)
library(ggseg)
library(ggplot2)
library(tidyr)
library(patchwork)
library(readr)
library(reshape2)

# Set the working directory to the path where your files are located
setwd("/Users/wang/Desktop/Research_projects/BBSC/Functional/fmri/derivatives/freesurfer-v7.2.0/")
# setwd("//Users/joeywang/Desktop/BBSC/Functional/fmri/derivatives/freesurfer-v7.2.0/")


# function to read and normalize the data by dividing the total intracranial volume
read_and_process_file <- function(filepath, cols, icv_col) {
  data <- read.delim(filepath)
  subset_data <- data[, cols]
  normalized_data <- subset_data / subset_data[, icv_col]
  
  return(normalized_data)
}

# define column names
cols <- c("BrainSegVol", "BrainSegVolNotVent", "lhCortexVol", "rhCortexVol", 
          "CortexVol", "lhCerebralWhiteMatterVol", "rhCerebralWhiteMatterVol",
          "CerebralWhiteMatterVol", "CSF", "EstimatedTotalIntraCranialVol")
icv_col <- "EstimatedTotalIntraCranialVol" # name for ICV

# read and normalize the data
sub1_long <- read_and_process_file("sub-1/1.statistics/long_pipeline/brain_seg_volume.txt", cols, icv_col)
sub1_normal <- read_and_process_file("sub-1/1.statistics/normal_pipeline/brain_seg_volume.txt", cols, icv_col)

sub2_long <- read_and_process_file("sub-2/1.statistics/long_pipeline/brain_seg_volume.txt", cols, icv_col)
sub2_normal <- read_and_process_file("sub-2/1.statistics/normal_pipeline/brain_seg_volume.txt", cols, icv_col)

sub3_long <- read_and_process_file("sub-3/1.statistics/long_pipeline/brain_seg_volume.txt", cols, icv_col)
sub3_normal <- read_and_process_file("sub-3/1.statistics/normal_pipeline/brain_seg_volume.txt", cols, icv_col)

datasets_long <- list(sub1_long, sub2_long, sub3_long)
datasets_normal <- list(sub1_normal, sub2_normal, sub3_normal)

morning_sub1 <-
  c(1:4, 6, 8, 9, 10, 12, 15, 17, 18, 20, 23, 24, 27, 30, 31, 32, 34, 35, 37)
afternoon_sub1 <-
  c(5, 7, 11, 13, 14, 16, 19, 21, 22, 25, 26, 28, 29, 33, 36, 38)

morning_sub2 <-
  c(1:4, 6, 7, 10:13, 16, 18, 21, 22, 24, 25, 27, 29, 31, 33, 34, 36, 37, 39)
afternoon_sub2 <-
  c(5, 8, 9, 14, 15, 17, 19, 20, 23, 26, 28, 30, 32, 35, 38, 40)

morning_sub3 <- c(1:6, 9, 11, 13, 15, 17, 18, 20, 22, 24)
afternoon_sub3 <- c(7, 8, 10, 12, 14, 16, 19, 21, 23, 25)

morning_indices_list <- list(morning_sub1, morning_sub2, morning_sub3)

# Define a function that insert several columns day, and time_seq
assign_data <- function(df, morning_indices) {
  df$day <- ifelse(seq_len(nrow(df)) %in% morning_indices, "morning", "afternoon")
  df$time_seq <- seq_len(nrow(df))
  return(df)
}

# Use mapply to change all the data in the datlist
datasets_long <-
  mapply(assign_data,
         datasets_long,
         morning_indices_list,
         SIMPLIFY = FALSE
  )

datasets_normal <-
  mapply(assign_data,
         datasets_normal,
         morning_indices_list,
         SIMPLIFY = FALSE
  )

# --load the qc data and replace 1 2 3 with 01 02 03
qc_data <- read_tsv(file = "group_T1w_v23.tsv") %>%
  mutate(bids_name = sub("ses-([1-9])([^0-9])", "ses-0\\1\\2", bids_name)) %>%
  arrange(bids_name)

qc_indx <-
  c(
    "cnr",
    "snr_total",
    "wm2max",
    "fber",
    "fwhm_avg",
    "rpve_gm",
    "rpve_csf",
    "cjv",
    "efc",
    "inu_med",
    "snr_csf"
  )
qc_2_control <- qc_data[, qc_indx]%>%
  mutate(
    subject = case_when(
      row_number() >= 1 & row_number() <= 38 ~ "sub1",
      row_number() >= 39 & row_number() <= 78 ~ "sub2",
      row_number() >= 79 ~ "sub3"
    )
  )

qc_2_control[c(34,36),"fber"] <- 35000 # an FBER=-1.0 indicates that there is 
#no signal outside the head mask, here very high quality image. 
  
datasets_long[[1]] <- cbind(datasets_long[[1]], qc_2_control[1:38, ])
datasets_normal[[1]] <- cbind(datasets_normal[[1]], qc_2_control[1:38, ])
  
datasets_long[[2]] <- cbind(datasets_long[[2]], qc_2_control[39:78, ])
datasets_normal[[2]] <- cbind(datasets_normal[[2]], qc_2_control[39:78, ])
  
datasets_long[[3]] <- cbind(datasets_long[[3]], qc_2_control[79:103, ])
datasets_normal[[3]] <- cbind(datasets_normal[[3]], qc_2_control[79:103, ]) 


# ##------------- change the time seq to the corresponding month

## -----for sub1
assign_time_seq <- function(data) {
  data <- data %>%
    mutate (
      month = case_when(
        row_number() == 1 ~ 1,
        row_number() >= 2 & row_number() <= 8 ~ 2,
        row_number() >= 9 & row_number() <= 10 ~ 3,
        row_number() >= 11 & row_number() <= 17 ~ 4,
        row_number() >= 18 & row_number() <= 22 ~ 5,
        row_number() >= 23 & row_number() <= 28 ~ 11,
        row_number() >= 29 & row_number() <= 33 ~ 12,
        row_number() >= 34 ~ 2
      )
    )
  return(data)
}

datasets_long[[1]] <- assign_time_seq(datasets_long[[1]])
datasets_normal[[1]] <- assign_time_seq(datasets_normal[[1]])

## ---- for sub2
assign_time_seq <- function(data) {
  data <- data %>%
    mutate (
      month = case_when(
        row_number() == 1 ~ 1,
        row_number() >= 2 & row_number() <= 9 ~ 2,
        row_number() >= 10 & row_number() <= 12 ~ 3,
        row_number() >= 13 & row_number() <= 19 ~ 4,
        row_number() >= 20 & row_number() <= 24 ~ 5,
        row_number() >= 25 & row_number() <= 31 ~ 11,
        row_number() >= 32 & row_number() <= 35 ~ 12,
        row_number() >= 36 ~ 2
      )
    )
  return(data)
}

datasets_long[[2]] <- assign_time_seq(datasets_long[[2]])
datasets_normal[[2]] <- assign_time_seq(datasets_normal[[2]])

## ---- for sub3
assign_time_seq <- function(data) {
  data <- data %>%
    mutate (
      month = case_when(
        row_number() >= 1 & row_number() <= 3 ~ 2,
        row_number() >= 4 & row_number() <= 6 ~ 3,
        row_number() >= 7 & row_number() <= 14 ~ 4,
        row_number() >= 15 & row_number() <= 19 ~ 5,
        row_number() >= 20 & row_number() <= 24 ~ 11,
        row_number() >= 25 ~ 12
      )
    )
  return(data)
}

datasets_long[[3]] <- assign_time_seq(datasets_long[[3]])
datasets_normal[[3]] <- assign_time_seq(datasets_normal[[3]])



convert_long_format <- function(data) {
  data_long <- select(data, -EstimatedTotalIntraCranialVol) %>%
    pivot_longer(
      cols = -c(
        subject,
        time_seq,
        month,
        day,
        cjv,
        efc,
        fwhm_avg,
        wm2max,
        fber,
        inu_med,
        cnr,
        snr_total,
        rpve_gm,
        snr_csf,
        rpve_csf
      ),
      names_to = "region",
      values_to = "volume"
    )
  return(data_long)
}

datasets_long_long <- lapply (datasets_long, convert_long_format)
datasets_normal_long <- lapply (datasets_normal, convert_long_format)

datasets_long_long_all <- rbind(datasets_long_long[[1]],datasets_long_long[[2]],datasets_long_long[[3]])
datasets_normal_long_all <- rbind(datasets_normal_long[[1]],datasets_normal_long[[2]],datasets_normal_long[[3]])


calculate_percentage_change <-
  function(data_set, cols_to_change, cols_to_keep) {
    percentage_change <- data_set[, cols_to_change] %>%
      sweep(
        2,
        colMeans(., na.rm = TRUE),
        FUN = function(x, mean) {
          ((x / mean) - 1) * 100
        }
      )
    cbind(percentage_change, data_set[, cols_to_keep])
  }

columns_to_change <- 1:9
columns_to_keep <- 10:25

percentage_change_longpipe <- do.call(
  rbind,
  mapply(
    calculate_percentage_change,
    datasets_long,
    MoreArgs = list(cols_to_change = columns_to_change, cols_to_keep = columns_to_keep),
    SIMPLIFY = FALSE
  )
)

percentage_change_normalpipe <- do.call(
  rbind,
  mapply(
    calculate_percentage_change,
    datasets_normal,
    MoreArgs = list(cols_to_change = columns_to_change, cols_to_keep = columns_to_keep),
    SIMPLIFY = FALSE
  )
)

percentage_change_longpipe_long <- convert_long_format (percentage_change_longpipe)%>%
  filter(!((time_seq == 1 | time_seq == 7) & subject == "sub3"))

percentage_change_normalpipe_long <- convert_long_format (percentage_change_normalpipe)%>%
  filter(!((time_seq == 1 | time_seq == 7) & subject == "sub3"))

sub_per_change_longpipe_long <-
  list(
    filter(percentage_change_longpipe_long, subject == "sub1"),
    filter(percentage_change_longpipe_long, subject == "sub2"),
    filter(percentage_change_longpipe_long, subject == "sub3")
  )

sub_per_change_normalpipe_long <-
  list(
    filter(percentage_change_normalpipe_long, subject == "sub1"),
    filter(percentage_change_normalpipe_long, subject == "sub2"),
    filter(percentage_change_normalpipe_long, subject == "sub3")
  )

# #### ---------------------------------------------Part IV Linear regression analysis
# library(lme4)
# library(MuMIn)
# library(lmerTest)
# 
# 
# modell_day <- lapply(sub_per_change_longpipe_long, function(df) {
#   
#   df$fber_log <- log(df$fber)
#   
#   normalized_df <- df %>%
#     mutate_at(
#       vars(
#         cjv,
#         efc,
#         fwhm_avg,
#         wm2max,
#         cnr,
#         snr_total,
#         snr_csf,
#         rpve_csf,
#         rpve_gm,
#         inu_med
#       ),
#       list(~ scale(.))
#     )
#   lmer(
#     volume ~ day * region + cjv + efc + fwhm_avg + wm2max + 
#       fber_log+ inu_med + cnr + snr_total + rpve_gm + (1 | month),
#     data = normalized_df
#   )
# })
# 
# summary(modell_day[[2]])
# r.squaredGLMM(modell_day[[2]])
# 
# p_adjusted_small <- lapply(modell_day, function(model) {
#   p_values <- summary(model)$coefficients[, "Pr(>|t|)"]
#   p_adjusted <- p.adjust(p_values, method = "holm")
#   p_adjusted_small <- p_adjusted[p_adjusted < 0.01]
#   print(p_adjusted_small)
# })
# 
# ## across all subj
# modell_day_all <- percentage_change_longpipe_long %>%
#   mutate(fber_log = log(fber)) %>% 
#   mutate_at(
#     vars(
#       cjv,
#       efc,
#       fwhm_avg,
#       wm2max,
#       inu_med,
#       cnr,
#       snr_total,
#       snr_csf,
#       rpve_csf,
#       rpve_gm,
#       fber_log
#     ),
#     list(~ scale(.))
#   ) %>%
#   lmer(
#     volume ~ day * region + cjv + efc + fwhm_avg + wm2max +
#       fber_log + inu_med + cnr + snr_total + rpve_gm + (1 | month),
#     data = .
#   )
# 
# r.squaredGLMM(modell_day_all)
# 
# p_values <- summary(modell_day_all)$coefficients[, "Pr(>|t|)"]
# 
# p_adjusted <- p.adjust(p_values, method = "holm")
# p_adjusted_small <- p_adjusted[p_adjusted < 0.01] %>%
#   print()
# 
# 
# ##-------------------normal pipeline
# 
# modell_day <- lapply(sub_per_change_normalpipe_long, function(df) {
#   
#   df$fber_log <- log(df$fber)
#   
#   normalized_df <- df %>%
#     mutate_at(
#       vars(
#         cjv,
#         efc,
#         fwhm_avg,
#         wm2max,
#         cnr,
#         snr_total,
#         snr_csf,
#         rpve_csf,
#         rpve_gm,
#         inu_med
#       ),
#       list(~ scale(.))
#     )
#   lmer(
#     volume ~ day * region + cjv + efc + fwhm_avg + wm2max + 
#       fber_log+ inu_med + cnr + snr_total + rpve_gm + (1 | month),
#     data = normalized_df
#   )
# })
# 
# summary(modell_day[[2]])
# r.squaredGLMM(modell_day[[2]])
# 
# p_adjusted_small <- lapply(modell_day, function(model) {
#   p_values <- summary(model)$coefficients[, "Pr(>|t|)"]
#   p_adjusted <- p.adjust(p_values, method = "holm")
#   p_adjusted_small <- p_adjusted[p_adjusted < 0.01]
#   print(p_adjusted_small)
# })
# 
# ## across all subj
# modell_day_all <- percentage_change_normalpipe_long %>%
#   mutate(fber_log = log(fber)) %>% 
#   mutate_at(
#     vars(
#       cjv,
#       efc,
#       fwhm_avg,
#       wm2max,
#       inu_med,
#       cnr,
#       snr_total,
#       snr_csf,
#       rpve_csf,
#       rpve_gm,
#       fber_log
#     ),
#     list(~ scale(.))
#   ) %>%
#   lmer(
#     volume ~ day * region + cjv + efc + fwhm_avg + wm2max +
#       fber_log + inu_med + cnr + snr_total + rpve_gm + (1 | month),
#     data = .
#   )
# 
# r.squaredGLMM(modell_day_all)
# 
# p_values <- summary(modell_day_all)$coefficients[, "Pr(>|t|)"]
# 
# p_adjusted <- p.adjust(p_values, method = "holm")
# p_adjusted_small <- p_adjusted[p_adjusted < 0.01] %>%
#   print()
# ## -----------------------------------about time of year of Longitudinal change
# 
# # within each subj
# modell_time <- lapply(sub_per_change_longpipe_long, function(df) {
#   df$fber_log <- log(df$fber)
#     normalized_df <- df %>%
#     mutate_at(
#       vars(
#         cjv,
#         efc,
#         fwhm_avg,
#         wm2max,
#         cnr,
#         snr_total,
#         snr_csf,
#         rpve_gm,
#         rpve_csf,
#         inu_med,
#         fber_log
#       ),
#       list(~ scale(.))
#     )
#   lmer(
#     volume ~ time_seq * region + cjv + efc + fwhm_avg + wm2max +
#       fber_log + inu_med + cnr + snr_total + rpve_gm  + (1 | month),
#     data = normalized_df
#   )
# })
# 
# summary(modell_time[[3]])
# r.squaredGLMM(modell_time[[2]])
# 
# p_adjusted_small <- lapply(modell_time, function(model) {
#   p_values <- summary(model)$coefficients[, "Pr(>|t|)"]
#   p_adjusted <- p.adjust(p_values, method = "holm")
#   p_adjusted_small <- p_adjusted[p_adjusted < 0.01]
#   print(p_adjusted_small)
# })

## across all subj
# modell_time_all <- percentage_change_longpipe_long %>%
#   mutate_at(
#     vars(
#       cjv,
#       efc,
#       fwhm_avg,
#       wm2max,
#       fber,
#       inu_med,
#       cnr,
#       snr_total,
#       rpve_gm
#     ),
#     list(~ scale(.))
#   ) %>%
#   lmer(
#     thickness_per_change ~ time_seq * region + cjv + efc + fwhm_avg + wm2max +
#       fber + inu_med + cnr + snr_total + rpve_gm + (1 | hemi),
#     data = .
#   )
# 
# r.squaredGLMM(modell_time_all)
# 
# p_values <- summary(modell_time_all)$coefficients[, "Pr(>|t|)"]
# 
# p_adjusted <- p.adjust(p_values, method = "BH")
# p_adjusted_small <- p_adjusted[p_adjusted < 0.01] %>%
#   print()

plot_day <- function (data2plot, brain) {
  plot <- ggplot(data2plot) +
    geom_violin(
      aes(
        x = day,
        y = volume,
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
        y = volume,
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
         y = brain)
    
    # Adjust ylim based on the value of brain
    if(brain == "CSF") {
      plot <- plot + 
        ylim(-25, 25) + 
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
      plot <- plot + ylim(-8, 8) +
        theme(axis.text.x = element_blank())
    }
  return(plot)
}

subj_id <- c("sub1", "sub2", "sub3")
brain_segment <- c("BrainSegVol", "BrainSegVolNotVent", "CortexVol", "CerebralWhiteMatterVol", "CSF")

plots_day_long_all <- lapply(brain_segment, function(brain) {
    data2plot <- percentage_change_longpipe_long %>%
      # filter(subject == subj) %>%
      filter(region == brain) %>%
      mutate(day = factor(day, levels = c("morning", "afternoon")))
    plots <- plot_day (data2plot, brain)
    
    return(plots)
  })

plots_day_normal_all <- lapply(brain_segment, function(brain) {
  data2plot <- percentage_change_normalpipe_long %>%
    # filter(subject == subj) %>%
    filter(region == brain) %>%
    mutate(day = factor(day, levels = c("morning", "afternoon")))
  plots <- plot_day (data2plot, brain)
  
  return(plots)
})

plots_day_long_sub <-
  lapply(sub_per_change_longpipe_long, function(sub_per_change) {
    plots_sub <- lapply(brain_segment, function(brain) {
      data2plot <- sub_per_change %>%
        # filter(subject == subj) %>%
        filter(region == brain) %>%
        mutate(day = factor(day, levels = c("morning", "afternoon")))
      plots <- plot_day (data2plot, brain)
      return(plots)
    })
    return(plots_sub)
  })


plots_day_normal_sub <-
  lapply(sub_per_change_normalpipe_long, function(sub_per_change) {
    plots_sub <- lapply(brain_segment, function(brain) {
      data2plot <- sub_per_change %>%
        # filter(subject == subj) %>%
        filter(region == brain) %>%
        mutate(day = factor(day, levels = c("morning", "afternoon")))
      plots <- plot_day (data2plot, brain)
      return(plots)
    })
    return(plots_sub)
  })

FigS7_long <- plots_day_long_sub[[1]][[1]] + plots_day_long_sub[[2]][[1]] + plots_day_long_sub[[3]][[1]] + plots_day_long_all[1] +
  plots_day_long_sub[[1]][[2]] + plots_day_long_sub[[2]][[2]] + plots_day_long_sub[[3]][[2]] + plots_day_long_all[2] +
  plots_day_long_sub[[1]][[3]] + plots_day_long_sub[[2]][[3]] + plots_day_long_sub[[3]][[3]] + plots_day_long_all[3] +
  plots_day_long_sub[[1]][[4]] + plots_day_long_sub[[2]][[4]] + plots_day_long_sub[[3]][[4]] + plots_day_long_all[4] +
  plots_day_long_sub[[1]][[5]] + plots_day_long_sub[[2]][[5]] + plots_day_long_sub[[3]][[5]] + plots_day_long_all[5] +
  plot_layout (ncol=4, axis_titles = "collect_x") +  plot_layout (axis_titles = "collect_y") 

FigS7_normal <- 
  plots_day_normal_sub[[1]][[1]] + plots_day_normal_sub[[2]][[1]] + plots_day_normal_sub[[3]][[1]] + plots_day_normal_all[1] +
  plots_day_normal_sub[[1]][[2]] + plots_day_normal_sub[[2]][[2]] + plots_day_normal_sub[[3]][[2]] + plots_day_normal_all[2] +
  plots_day_normal_sub[[1]][[3]] + plots_day_normal_sub[[2]][[3]] + plots_day_normal_sub[[3]][[3]] + plots_day_normal_all[3] +
  plots_day_normal_sub[[1]][[4]] + plots_day_normal_sub[[2]][[4]] + plots_day_normal_sub[[3]][[4]] + plots_day_normal_all[4] +
  plots_day_normal_sub[[1]][[5]] + plots_day_normal_sub[[2]][[5]] + plots_day_normal_sub[[3]][[5]] + plots_day_normal_all[5] +
  plot_layout (ncol=4, axis_titles = "collect_x") +  plot_layout (axis_titles = "collect_y") 

ggsave("/Users/wang/Library/CloudStorage/OneDrive-UniversityofBergen/Desktop/UiB/Manuscripts/BBSC/4.T1w_Image/Figures/FigS7.png",
       plot = FigS7_long, width = 10, height = 12, units = 'in', dpi = 100)

ggsave("/Users/wang/Library/CloudStorage/OneDrive-UniversityofBergen/Desktop/UiB/Manuscripts/BBSC/4.T1w_Image/Figures/FigS7_hd.png",
       plot = FigS7_long, width = 10, height = 12, units = 'in', dpi = 300)


# ggsave("/Users/wang/Library/CloudStorage/OneDrive-UniversityofBergen/Desktop/UiB/Manuscripts/BBSC/4.T1w_Image/Figures/FigS7_normal
#        .png",
#        plot = FigS7_normal, width = 10, height = 12, units = 'in', dpi = 300, compression = "lzw")


# ggsave("/Users/joeywang/Library/CloudStorage/OneDrive-UniversityofBergen/Desktop/UiB/Manuscripts/BBSC/4.T1w_Image/Figures/FigS7.png",
       # plot = FigS7, width = 10, height = 12, units = 'in', dpi = 300, compression = "lzw")

# ggsave("/Users/joeywang/Library/CloudStorage/OneDrive-UniversityofBergen/Desktop/UiB/Manuscripts/BBSC/4.T1w_Image/Figures/FigS8.png",
# plot = FigS8, width = 10, height = 12, units = 'in', dpi = 300, compression = "lzw")

