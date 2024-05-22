# This code is to 
# 1. calculate the CV of different brain apartments
# 2. add NA to even the session numbers of all subjects
# 3. plot the percentage changes of brain apartments derived from long and normal pipelines

# WANG. 20-Feb-2024

# Remove all objects created before to prevent clustering
rm(list=ls())

# Load required packages

library(berryFunctions)
library(dplyr)

# Set the working directory to the path where your files are located
setwd("/Users/wang/Desktop/Research_projects/BBSC/Functional/fmri/derivatives/freesurfer-v7.2.0/")
# setwd("//Users/joeywang/Desktop/BBSC/Functional/fmri/derivatives/freesurfer-v7.2.0/")

##---------------------------------------Load the data
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

## -------------------------------------Calculate the CV
# function to calculate the coefficient of variation
cv <- function(x) { 100 * (sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)) }

# define calculate cv function
calculate_cv <- function(data) {
  apply(data, 2, cv)
}

for (sub in 1:3) {
  for (pipeline in c("long", "normal")) {
    # calculate CV
    data_name <- paste0("sub", sub, "_", pipeline)
    assign(paste0(data_name, "_cv"), calculate_cv(get(data_name)))
  }
}
brain_cv <-data.frame(label=cols, 
                      sub1_long= sub1_long_cv, sub1_normal= sub1_normal_cv, 
                      sub2_long= sub2_long_cv, sub2_normal= sub2_normal_cv, 
                      sub3_long= sub3_long_cv, sub3_normal= sub3_normal_cv) %>%
  mutate_if(is.numeric, ~round(., 2))

## ------------------------------------ Make the data length the same
# insert nan values to make it even so, 1:24 sessions are from Jan.to May.
sub1_long <- insertRows(sub1_long, 23:24 , new = NA)
sub1_normal <- insertRows(sub1_normal, 23:24 , new = NA)
sub3_long <- insertRows(sub3_long, 20:24 , new = NA)
sub3_normal <- insertRows(sub3_normal, 20:24 , new = NA)

# insert nan values to make the gap, 25:28 sessions are from Jun.to Oct.
sub1_long <- insertRows(sub1_long, 25:28 , new = NA)
sub2_long <- insertRows(sub2_long, 25:28 , new = NA)
sub3_long <- insertRows(sub3_long, 25:28 , new = NA)

sub1_normal <- insertRows(sub1_normal, 25:28 , new = NA)
sub2_normal <- insertRows(sub2_normal, 25:28 , new = NA)
sub3_normal <- insertRows(sub3_normal, 25:28 , new = NA)

# insert nan values to make it even, 29:39 sessions are from Nov.to Dec.
sub3_long[35:39,] <- NA
sub3_normal[35:39,] <- NA

# insert nan values to make it even, 40:41 sessions are from jan. 2022
sub1_long <- insertRows(sub1_long, 40:41 , new = NA)
sub2_long <- insertRows(sub2_long, 40:41 , new = NA)
sub3_long[40:41,] <- NA

sub1_normal <- insertRows(sub1_normal, 40:41 , new = NA)
sub2_normal <- insertRows(sub2_normal, 40:41 , new = NA)
sub3_normal[40:41,] <- NA

# insert nan values to make it even, 42:46 sessions are from Feb. 2022
sub3_long[42:46,] <- NA
sub3_normal[42:46,] <- NA

# calculate mean and limits
sub1_long_mean <- colMeans(sub1_long,na.rm=TRUE)
sub2_long_mean <- colMeans(sub2_long,na.rm=TRUE)
sub3_long_mean <- colMeans(sub3_long,na.rm=TRUE)

sub1_normal_mean <- colMeans(sub1_normal,na.rm=TRUE)
sub2_normal_mean <- colMeans(sub2_normal,na.rm=TRUE)
sub3_normal_mean <- colMeans(sub3_normal,na.rm=TRUE)

## -------------------------------------------Plot percentage change
library(ggplot2)
# set a custom theme for each figures
custom_theme <- theme_bw(base_size = 12) +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(size = 1, color = "black"),
    axis.text = element_text(size = 14, face = "bold", family = "sans", color = "black"),
    axis.title = element_text(size = 16, face = "bold"),
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  )

# initiate the variables
column_indices <- c("BrainSegVolNotVent", "CortexVol", "CerebralWhiteMatterVol")
y_labels <- c("Total Brain", "Grey Matter", "White Matter")


# Plot the long pipeline percentage change
plots_long_cv <- lapply(seq_along(column_indices), function(i) {
  indx <- column_indices[i]
  
  data_sub1 <- data.frame(x = seq(1, 46), y = (sub1_long[[indx]] / sub1_long_mean[[indx]] - 1) * 100)
  data_sub2 <- data.frame(x = seq(1, 46), y = (sub2_long[[indx]] / sub2_long_mean[[indx]] - 1) * 100)
  data_sub3 <- data.frame(x = seq(1, 46), y = (sub3_long[[indx]] / sub3_long_mean[[indx]] - 1) * 100)
  
  ggplot() +
    geom_point(data=data_sub1, aes(x = x, y = y), color = "#D40000") +
    geom_point(data=data_sub2, aes(x = x, y = y), color = "#377EB8") +
    geom_point(data=data_sub3, aes(x = x, y = y), color = "#4DAF4A") +
    geom_hline(yintercept = c(2, -2), linetype = "dashed", color = "grey", size = 1) +
    geom_segment(aes(x = 0.5, xend = 0.5, y = -2, yend = 2), linetype = "solid", color = "#6A3D9A", size = 1) +
    geom_segment(aes(x = 24.5, xend = 24.5, y = -2, yend = 2), linetype = "dashed", color = "#6A3D9A", size = 1) +
    geom_segment(aes(x = 28.5, xend = 28.5, y = -2, yend = 2), linetype = "dashed", color = "#6A3D9A", size = 1) +
    geom_segment(aes(x = 39.5, xend = 39.5, y = -2, yend = 2), linetype = "solid", color = "#6A3D9A", size = 1) +
    geom_segment(aes(x = 41.5, xend = 41.5, y = -2, yend = 2), linetype = "solid", color = "#FFA500", size = 1) +
    geom_segment(aes(x = 46.5, xend = 46.5, y = -2, yend = 2), linetype = "solid", color = "#FFA500", size = 1) +
    xlab("") +
    ylab(y_labels[i]) +
    scale_y_continuous(limits = c(-6, 6), breaks = seq(-5, 5, by = 5)) +
    custom_theme
})

# Plot the normal pipeline percentage change
plots_normal_cv <- lapply(seq_along(column_indices), function(i) {
  indx <- column_indices[i]
  
  data_sub1 <- data.frame(x = seq(1, 46), y = (sub1_normal[[indx]] / sub1_normal_mean[[indx]] - 1) * 100)
  data_sub2 <- data.frame(x = seq(1, 46), y = (sub2_normal[[indx]] / sub2_normal_mean[[indx]] - 1) * 100)
  data_sub3 <- data.frame(x = seq(1, 46), y = (sub3_normal[[indx]] / sub3_normal_mean[[indx]] - 1) * 100)
  
  ggplot() +
    geom_point(data=data_sub1, aes(x = x, y = y), color = "#D40000") +
    geom_point(data=data_sub2, aes(x = x, y = y), color = "#377EB8") +
    geom_point(data=data_sub3, aes(x = x, y = y), color = "#4DAF4A") +
    geom_hline(yintercept = c(2, -2), linetype = "dashed", color = "grey", size = 1) +
    geom_segment(aes(x = 0.5, xend = 0.5, y = -2, yend = 2), linetype = "solid", color = "#6A3D9A", size = 1) +
    geom_segment(aes(x = 24.5, xend = 24.5, y = -2, yend = 2), linetype = "dashed", color = "#6A3D9A", size = 1) +
    geom_segment(aes(x = 28.5, xend = 28.5, y = -2, yend = 2), linetype = "dashed", color = "#6A3D9A", size = 1) +
    geom_segment(aes(x = 39.5, xend = 39.5, y = -2, yend = 2), linetype = "solid", color = "#6A3D9A", size = 1) +
    geom_segment(aes(x = 41.5, xend = 41.5, y = -2, yend = 2), linetype = "solid", color = "#FFA500", size = 1) +
    geom_segment(aes(x = 46.5, xend = 46.5, y = -2, yend = 2), linetype = "solid", color = "#FFA500", size = 1) +
    xlab("") +
    ylab("") +
    scale_y_continuous(limits = c(-6, 6), breaks = seq(-5, 5, by = 5)) +
    custom_theme+
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_line(),
          axis.title.y = element_text(vjust = 0))
})

# add the plots titles
plots_long_cv[[1]] <- plots_long_cv[[1]] + 
  ggtitle("Longitudinal") + 
  theme(plot.title = element_text(size = 18, face = "bold", family = "sans", color = "black", hjust =0.5))

plots_normal_cv[[1]] <- plots_normal_cv[[1]] + 
  ggtitle("Independent") +
  theme(plot.title = element_text(size = 18, face = "bold", family = "sans", color = "black", hjust =0.5))

# arrange the subplots into one plot
library(patchwork)
plots_per_change <- plots_long_cv[[1]] + plots_normal_cv[[1]] + 
  plots_long_cv[[2]] + plots_normal_cv[[2]] + 
  plots_long_cv[[3]] + plots_normal_cv[[3]] + plot_layout(ncol = 2)

# save the plot
 save(plots_per_change, file = "plots_per_change.RData")

