# This code is to
# 1. load the subcortex data
# 2. calculate and plot the CV of the subcortex
# 2. plot the mean, ordered mean; cv and percentage change 

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

##--------------------------------Part I Load the data and add factors
# define the file path
file_paths <- c(
  "sub-1/1.statistics/long_pipeline/brain_seg_volume.txt",
  "sub-2/1.statistics/long_pipeline/brain_seg_volume.txt",
  "sub-3/1.statistics/long_pipeline/brain_seg_volume.txt"
)

fs_data <-
  read_freesurfer_stats("sub-1/sub-1_ses-1_T1w//stats/aseg.stats")

data_sets <- lapply(file_paths, function(path) {
  data <- read.delim(path)[, 2:46]
  colnames(data) = fs_data[["label"]]
  
  colnames(data)[colnames(data) == "Left-Thalamus"] <-
    "Left-Thalamus-Proper"
  colnames(data)[colnames(data) == "Right-Thalamus"] <-
    "Right-Thalamus-Proper"
  colnames(data)[colnames(data) == "3rd-Ventricle"] <-
    "x3rd-ventricle"
  colnames(data)[colnames(data) == "4th-Ventricle"] <-
    "x4th-ventricle"
  colnames(data)[colnames(data) == "5th-Ventricle"] <-
    "x5th-ventricle"
  
  data <- data %>%
    select(-c(33:39)) %>%
    select(-c(
      "Left-vessel",
      "Right-vessel",
      "Left-Inf-Lat-Vent",
      "Right-Inf-Lat-Vent",
      "Optic-Chiasm",
      "Left-choroid-plexus",
      "Right-choroid-plexus"
    )) # delete some columns
  
  data$CC <-
    rowSums(data[, c(
      "CC_Posterior",
      "CC_Mid_Posterior",
      "CC_Central",
      "CC_Mid_Anterior",
      "CC_Anterior"
    )], na.rm = TRUE)
  
  data$'Left-Total' <-
    rowSums(data[, c(
      "Left-Thalamus-Proper",
      "Left-Putamen",
      "Left-Hippocampus",
      "Left-Caudate",
      "Left-Pallidum",
      "Left-Amygdala",
      "Left-Accumbens-area"
    )], na.rm = TRUE)
  
  data$'Right-Total' <-
    rowSums(data[, c(
      "Right-Thalamus-Proper",
      "Right-Putamen",
      "Right-Hippocampus",
      "Right-Caudate",
      "Right-Pallidum",
      "Right-Amygdala",
      "Right-Accumbens-area"
    )], na.rm = TRUE)
  
  data$'total' <-
    rowSums(data[, c(
      "Right-Total",
      "Left-Total"
    )], na.rm = TRUE)
  return(data)
})

names(data_sets) <- c("sub1", "sub2", "sub3")

# exclude the session 1 and 7 in sub3
# data_sets$sub3[c(1,7), ] <- NaN # pay attention to the saving files


# Convert each data frame to long format and calculate percentage change
sub_percent_long <- lapply(names(data_sets), function(name) {
  df <- data_sets[[name]]
  
  # Calculate the percentage change for all columns
  percentage <-
    sweep(df, 2, STATS = colMeans(df, na.rm = TRUE), FUN = "/") %>%
    mutate(across(everything(), ~ (. - 1) * 100)) %>%
    # Add the subject column using the name
    mutate(subject = name)
  
  # Convert the wide data to long format
  percentage_long <-
    pivot_longer(
      percentage,
      cols = !c(subject),
      names_to = "region",
      values_to = "percent_change"
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
    )
  
  return(percentage_long)
})

names(sub_percent_long) <- c("sub1_long", "sub2_long", "sub3_long")

####--------------------------------Part II: Calculate the coefficient of variation for each subset and column

# function to calculate the coefficient of variation
cv <- function(x) {
    100 * (sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE))
  }
calculate_cv <- function(data) {
  apply(data, 2, cv)
}

sub_cv_volume_sub <- lapply(data_sets, function(df) {
  result <- calculate_cv(df)
  round(result, 2)
})

sub_mean_volume_sub <- lapply(data_sets, function (df) {
  sub_mean <- round(colMeans(df, na.rm = TRUE), 2)
})

sub_mean_sd_volume_sub <- lapply(data_sets, function (df) {
  sub_mean <- round(colMeans(df, na.rm = TRUE), 2)
  sub_sd <- round(apply(df, 2, sd, na.rm = TRUE), 2)
  results <- data.frame(sub_mean = sub_mean, sub_sd = sub_sd)  # Explicitly make results a data frame
  
  return(results)
})

# save all the variables
save(data_sets, file="sub_data_sets.RData")
save(sub_percent_long, file = "sub_percent_long.RData")
save(sub_mean_sd_volume_sub, file = "sub_mean_sd_volume_sub.Rdata")
write.csv(sub_cv_volume_sub, "sub_cv_volume_sub.csv", row.names = TRUE)
write.csv(sub_mean_volume_sub, "sub_mean_volume_sub.csv", row.names = TRUE)

# # save all the variables
# data_sets_qc <- data_sets
# save(data_sets_qc, file="sub_data_sets_qc.RData")
# 
# sub_percent_long_qc <- sub_percent_long
# save(sub_percent_long_qc, file = "sub_percent_long_qc.RData")
# 
# sub_cv_volume_sub_qc <- sub_cv_volume_sub
# save(sub_cv_volume_sub_qc, file = "sub_cv_volume_sub_qc.RData")
# write.csv(sub_cv_volume_sub_qc, "sub_cv_volume_sub_qc.csv", row.names = TRUE)
# 
# sub_mean_volume_sub_qc <- sub_mean_volume_sub
# save(sub_mean_volume_sub_qc, file = "sub_mean_volume_sub_qc.RData")
# write.csv(sub_mean_volume_sub_qc, "sub_mean_volume_sub_qc.csv", row.names = TRUE)
# 
# sub_mean_sd_volume_sub_qc <- sub_mean_sd_volume_sub
# save(sub_mean_sd_volume_sub_qc, file = "sub_mean_sd_volume_sub_qc.Rdata")


