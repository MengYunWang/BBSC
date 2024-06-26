# This code is to 
# 1. Calculate and plot the CV of the sMRI data
# 2. Plot the percentage change of each region with box plot

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

## --------------------------------Load the data
# define the file path
file_paths <- c(
  "sub-1/1.statistics/long_pipeline/lh_volume_dk.txt",
  "sub-1/1.statistics/long_pipeline/rh_volume_dk.txt",
  "sub-2/1.statistics/long_pipeline/lh_volume_dk.txt",
  "sub-2/1.statistics/long_pipeline/rh_volume_dk.txt",
  "sub-3/1.statistics/long_pipeline/lh_volume_dk.txt",
  "sub-3/1.statistics/long_pipeline/rh_volume_dk.txt"
)

data_sets <- lapply(file_paths, function(path) {
  data <- read.delim(path)[, 2:35]
  data$total <- rowSums(data, na.rm = TRUE)
  return(data)
})
names(data_sets) <- c("sub1_lh", "sub1_rh", "sub2_lh", "sub2_rh", "sub3_lh", "sub3_rh")

list2env(data_sets, envir = .GlobalEnv)

sub3_lh[c(1,7), ] <- NA # delete the 1st and 7th
sub3_rh[c(1,7), ] <- NA # delete the 1st and 7th collection because of the qc

##------------------Part I: Calculate the coefficient of variation for each subset and column

# function to calculate the coefficient of variation
cv <- function(x) { 100 * (sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)) }
calculate_cv <- function(data) {
  apply(data, 2, cv)
}

# function to read the freesurfer data and add prefix to label
read_and_label_fs_data <- function(side, path) {
  read_freesurfer_stats(path)  %>%
    add_row(label = "total") %>%
    mutate(across(-label, ~NA)) %>%
    mutate(label = paste0(side, "_", label))
}

# funtion to plot the cv with ggseg
plot_cv_volume <- function(cv_volume_data) {
  ggplot(cv_volume_data) +
    geom_brain(atlas = dk, 
               position = position_brain(hemi ~ side),
               aes(fill = sub_cv)) +
    scale_fill_viridis_c(option = "magma", direction = 1, limits = c(0, 8), breaks = c(2, 4, 6), labels = c("2", "4", "6"))+
    theme_void()
}

# function to plot the mean data with ggseg
plot_mean_volume <- function(mean_volume_data) {
  ggplot(mean_volume_data) +
    geom_brain(atlas = dk, 
               #position = position_brain(hemi ~ side),
               aes(fill = sub_mean)) +
    scale_fill_viridis_c(option = "magma", direction = 1, limits = c(0, 26600), breaks = c(10000, 15000, 20000), labels = c("10000", "15000", "20000")) +
    theme_void()
}

# calculate the cv and mean, and then plot them with ggseg
plots_brain_cv <- list()
plots_brain_mean <- list()
for (sub in 1:3) {
  for (hemi in c("lh", "rh")) {
    # calculate CV
    data_name <- paste0("sub", sub, "_", hemi)
    assign(paste0(data_name, "_cv"), calculate_cv(get(data_name)))
  }
  
  # read labels from a freesurfer data
  fs_data_lh <- read_and_label_fs_data("lh", "sub-1/sub-1_ses-1_T1w//stats/lh.aparc.stats")
  fs_data_rh <- read_and_label_fs_data("rh", "sub-1/sub-1_ses-1_T1w//stats/lh.aparc.stats")
  
  # combine the labels and CVs together
  sub_cv_volume <- data.frame(
    label = rbind(fs_data_lh[, "label"], fs_data_rh[, "label"]),
    sub_cv = c(get(paste0("sub", sub, "_lh_cv")), get(paste0("sub", sub, "_rh_cv")))
  )
  
  sub_mean_volume <- data.frame(
    label = rbind(fs_data_lh[, "label"], fs_data_rh[, "label"]),
    sub_mean = c(colMeans(get(paste0("sub", sub, "_", hemi)), na.rm = TRUE), colMeans(get(paste0("sub", sub, "_", hemi)), na.rm = TRUE))
  )
  
  # plot the brain data with ggseg
  plots_brain_cv[[sub]] <- sub_cv_volume %>%
    slice(-c(35, 70)) %>% # exclude the average row, because it doesnot belong to the DK atlas
    plot_cv_volume()
  
  plots_brain_mean[[sub]] <- sub_mean_volume %>%
    slice(-c(35, 70)) %>% # exclude the average row, because it doesnot belong to the DK atlas
    plot_mean_volume()
}

# save(plots_brain_cv, file = "plots_brain_volume_cv_orig.RData")
# save(plots_brain_mean, file = "plots_brain_volume_mean_orig.RData") # save it for reuse

plots_brain_cv_qc <- plots_brain_cv
plots_brain_mean_qc <- plots_brain_mean
save(plots_brain_cv_qc, file = "plots_brain_volume_cv_qc.RData")
save(plots_brain_mean_qc, file = "plots_brain_volume_mean_qc.RData") # save it for reuse


# ----------------------------------calculate the cv and mean across all sessions
fs_data <- read_freesurfer_stats(paste0("sub-1/sub-1_ses-1_T1w/stats/lh.aparc.stats")) %>%
  add_row(label = "total") %>%
  mutate(across(-label, ~NA))

volume_cv <-data.frame(label=fs_data[["label"]], 
                          sub1_lh= sub1_lh_cv, sub1_rh= sub1_rh_cv, 
                          sub2_lh= sub2_lh_cv, sub2_rh= sub2_rh_cv, 
                          sub3_lh= sub3_lh_cv, sub3_rh= sub3_rh_cv) %>%
  mutate_if(is.numeric, ~round(., 2))

volume_mean <-data.frame(label=fs_data[["label"]], 
                            sub1_lh= colMeans(sub1_lh, na.rm = TRUE), sub1_rh= colMeans(sub1_rh, na.rm = TRUE), 
                            sub2_lh= colMeans(sub2_lh, na.rm = TRUE), sub2_rh= colMeans(sub2_rh, na.rm = TRUE), 
                            sub3_lh= colMeans(sub3_lh, na.rm = TRUE), sub3_rh= colMeans(sub3_rh, na.rm = TRUE)
)%>%
  mutate_if(is.numeric, ~round(., 2))

volume_mean_sd <- data.frame(label=fs_data[["label"]]) %>%
  mutate(
    # Calculating means; adapt if 'sub1_lh', etc., are more complex
    sub1_lh_mean = round(colMeans(sub1_lh, na.rm = TRUE), 2), 
    sub1_rh_mean = round(colMeans(sub1_rh, na.rm = TRUE), 2), 
    sub2_lh_mean = round(colMeans(sub2_lh, na.rm = TRUE), 2), 
    sub2_rh_mean = round(colMeans(sub2_rh, na.rm = TRUE), 2), 
    sub3_lh_mean = round(colMeans(sub3_lh, na.rm = TRUE), 2), 
    sub3_rh_mean = round(colMeans(sub3_rh, na.rm = TRUE), 2),
    # Calculating standard deviations
    sub1_lh_sd = round(apply(sub1_lh, 2, sd, na.rm = TRUE), 2),
    sub1_rh_sd = round(apply(sub1_rh, 2, sd, na.rm = TRUE), 2),
    sub2_lh_sd = round(apply(sub2_lh, 2, sd, na.rm = TRUE), 2),
    sub2_rh_sd = round(apply(sub2_rh, 2, sd, na.rm = TRUE), 2),
    sub3_lh_sd = round(apply(sub3_lh, 2, sd, na.rm = TRUE), 2),
    sub3_rh_sd = round(apply(sub3_rh, 2, sd, na.rm = TRUE), 2)
  )

volume_mean_ranked <- volume_mean %>%
  slice(-c(35, 70)) %>% #exclude the average row
  mutate(across(where(is.numeric), ~as.integer(rank(., ties.method = "first")), .names = "{.col}_rank"))


##------------------------------------ Part II: Plot the boxplot of each region (Fig. 3 in the manuscript)

# function that calculate the percentage change and box plot them
process_subject <- function(subj) {
  fs_data <- read_freesurfer_stats(paste0("sub-1/sub-1_ses-1_T1w/stats/lh.aparc.stats")) %>%
    add_row(label = "total") %>%
    mutate(across(-label, ~NA))
  
  lh_percent <- get(paste0("sub", subj, "_lh")) %>%
    sweep(2, colMeans(., na.rm = TRUE), FUN = function(x, mean) ((x / mean) - 1)*100) %>%
    mutate(hemi = "lh", subject = paste0("sub", subj))
  
  rh_percent <- get(paste0("sub", subj, "_rh")) %>%
    sweep(2, colMeans(., na.rm = TRUE), FUN = function(x, mean) ((x / mean) - 1)*100) %>%
    mutate(hemi = "rh", subject = paste0("sub", subj))
  
  colnames(lh_percent)[1:(ncol(lh_percent)-2)] <- fs_data[["label"]]
  colnames(rh_percent)[1:(ncol(rh_percent)-2)] <- fs_data[["label"]]
  
  percent <- bind_rows(lh_percent, rh_percent)
  
  percent_long <- percent %>% 
    pivot_longer(cols = -c(hemi, subject), names_to = "region", values_to = "percent_change")%>%
    filter(region != "total") %>%
  
  return(percent_long)
}

# process each sub and combing them together
all_subjects_percent_long <- bind_rows(
  process_subject(1),
  process_subject(2),
  process_subject(3)
)
# appendix subID at end of the column 
all_subjects_percent_long$panel_id <- as.integer(factor(all_subjects_percent_long$subject))
all_subjects_percent_long <- all_subjects_percent_long %>% 
  mutate(panel_id = as.integer(factor(subject)))

# split the datafram according to the subID
split_data <- split(all_subjects_percent_long, all_subjects_percent_long$subject)

# plot and save each plot for each sub
plots_list_box <- lapply(split_data, function(data) {
  ggplot(data, aes(x = percent_change, y = region, fill = hemi)) +
    geom_boxplot(position = position_dodge(width = 0.75)) +
    scale_x_continuous(limits = c(-18, 18), breaks = c(-5, 0, 5)) +
    labs(x = "", y = "") +
    # guides(fill = FALSE) +
    theme_minimal() +
    theme(
      strip.background = element_blank(),
      strip.text = element_blank(),
      axis.text.x = element_text(size = 18, family = "Arial", face = "bold", colour = "black"),
      axis.text.y = element_text(size = 18, family = "Arial", face = "bold", colour = "black"),
      axis.ticks.length = unit(0.2, "cm"),
      axis.ticks.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    geom_vline(xintercept = c(-5, 5), linetype = "dashed", color = "green") +
    geom_vline(xintercept = c(-15, 15), linetype = "dashed", color = "red")
})


# save(plots_list_box, file = "plots_brain_volume_per_change_orig.RData")

plots_list_box_qc <- plots_list_box
save(plots_list_box_qc, file = "plots_brain_volume_per_change_qc.RData")

