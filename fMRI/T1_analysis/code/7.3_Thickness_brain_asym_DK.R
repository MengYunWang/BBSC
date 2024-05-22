# This code is to 
# 1. calculate and plot the CV of the sMRI data
# 2. plot the percentage change of each region with box plot
# 3. combing these two steps together to one figure

# WANG. 26-Feb-2024

# Remove all objects created before to prevent clash
rm(list=ls())

library(dplyr)
library(ggseg)
library(ggplot2)
library(tidyr)
library(patchwork)


# Set the working directory to the path where your files are located
setwd("/Users/wang/Desktop/Research_projects/BBSC/Functional/fmri/derivatives/freesurfer-v7.2.0/")
# setwd("//Users/joeywang/Desktop/BBSC/Functional/fmri/derivatives/freesurfer-v7.2.0/")

#### --------------------------------Part I: Load the data
# define the file path
file_paths <- c(
  "sub-1/1.statistics/long_pipeline/lh_thickness_dk.txt",
  "sub-1/1.statistics/long_pipeline/rh_thickness_dk.txt",
  "sub-2/1.statistics/long_pipeline/lh_thickness_dk.txt",
  "sub-2/1.statistics/long_pipeline/rh_thickness_dk.txt",
  "sub-3/1.statistics/long_pipeline/lh_thickness_dk.txt",
  "sub-3/1.statistics/long_pipeline/rh_thickness_dk.txt"
)

fs_data <- read_freesurfer_stats("sub-1/sub-1_ses-1_T1w//stats/lh.aparc.stats")
data_sets <- lapply(file_paths, function(path) {
  read.delim(path)[, 2:35]
})

data_sets <- lapply(data_sets, function(df) {
  colnames(df) = fs_data[["label"]]
  return(df)
})

names(data_sets) <- c("sub1_lh", "sub1_rh", "sub2_lh", "sub2_rh", "sub3_lh", "sub3_rh")


brain_asym <- lapply(names(data_sets)[seq(1, length(data_sets), 2)], function(subject_prefix) {
  base_name <- gsub("_lh", "", subject_prefix)
  lh <- data_sets[[paste0(base_name, "_lh")]]
  rh <- data_sets[[paste0(base_name, "_rh")]]
  diff <- lh - rh
  mean_lh_rh <- (lh + rh)/2
  normalized_diff <- diff / mean_lh_rh
  normalized_diff
})
names(brain_asym) <- c("sub1", "sub2", "sub3")



