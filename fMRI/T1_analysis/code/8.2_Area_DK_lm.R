# This code is to
# 1. load the data and add more variables
# 2. calculate and plot the average ordered
# 3. calculate the percentage change of the area
# 4. linear mixed regression to day and time effects
# 5. Plot the brain regions that have day or time effects

# WANG. 04-March-2024

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
setwd("/Users/wang/Desktop/Research_projects/BBSC/Functional/fmri/derivatives/freesurfer-v7.2.0/")
# setwd("//Users/joeywang/Desktop/BBSC/Functional/fmri/derivatives/freesurfer-v7.2.0/")

#### -------------------------------------------------------------Part I Load the data
# define the file path
file_paths <- c(
  "sub-1/1.statistics/long_pipeline/lh_area_dk.txt",
  "sub-1/1.statistics/long_pipeline/rh_area_dk.txt",
  "sub-2/1.statistics/long_pipeline/lh_area_dk.txt",
  "sub-2/1.statistics/long_pipeline/rh_area_dk.txt",
  "sub-3/1.statistics/long_pipeline/lh_area_dk.txt",
  "sub-3/1.statistics/long_pipeline/rh_area_dk.txt"
)

data_sets <- lapply(file_paths, function(path) {
  data <- read.delim(path)[, 2:35]
  data$total <- rowSums(data, na.rm = TRUE)
  return(data)
})
names(data_sets) <-
  c("sub1_lh",
    "sub1_rh",
    "sub2_lh",
    "sub2_rh",
    "sub3_lh",
    "sub3_rh")

data_sets[["sub3_lh"]][c(1,7), ] <- NA # delete the 1st and 7th, pay attention to the files that need to be saved
data_sets[["sub3_rh"]][c(1,7), ] <- NA # delete the 1st and 7th sessions because of the qc

##----------------------add more factors in the columns
# load the labels
fs_data <- read_freesurfer_stats("sub-1/sub-1_ses-1_T1w//stats/lh.aparc.stats") %>%  
  add_row(label = "total") %>%
  mutate(across(-label, ~NA))

# add the morning and after indices
morning_sub1 <- c(1:4, 6, 8, 9, 10, 12, 15, 17, 18, 20, 23, 24, 27, 30, 31, 32, 34, 35, 37)
afternoon_sub1 <- c(5, 7, 11, 13, 14, 16, 19, 21, 22, 25, 26, 28, 29, 33, 36, 38)

morning_sub2 <- c (1:4, 6, 7, 10:13, 16, 18, 21, 22, 24, 25, 27, 29, 31, 33, 34, 36, 37, 39)
afternoon_sub2 <- c(5, 8, 9, 14, 15, 17, 19, 20, 23, 26, 28, 30, 32, 35, 38, 40)

morning_sub3 <- c(1:6, 9, 11, 13, 15, 17, 18, 20, 22, 24)
afternoon_sub3 <- c(7, 8, 10, 12, 14, 16, 19, 21, 23, 25)

morning_indices_list <- list(morning_sub1, morning_sub2, morning_sub3)
morning_indices_list <- rep(morning_indices_list, each = 2)

# add the day length in minutes at the data collection
day_length_sub1 <- c(200, 283, 369, 427, 427, 458, 463, 473, 552, 586, 739, 773, 
                     778, 812, 818, 853, 858, 885, 888, 927, 930, 948, 322, 270, 
                     260, 221, 189, 189, 159, 151, 133, 120, 120, 360, 421, 421, 
                     458, 458)

day_length_sub2 <- c(244, 275, 301, 369, 427, 437, 458, 463, 468, 538, 552, 579, 
                     731, 739, 778, 782, 812, 847, 853, 888, 891, 927, 930, 949, 
                     260, 260, 226, 226, 189, 189, 163, 159, 133, 120, 120, 340, 
                     421, 421, 458, 458)

day_length_sub3 <- c(283, 427, 463, 552, 586, 639, 739, 739, 773, 778, 808, 812, 
                     847, 853, 885, 930, 933, 945, 948, 270, 260, 226, 223, 163, 
                     159)
day_length <- list(day_length_sub1, day_length_sub2, day_length_sub3) %>% 
  rep(each = 2)

# function to insert several columns day, hemi, and time_seq
assign_data <- function(df, morning_indices, day_length, df_name) {
  colnames(df) <- fs_data[["label"]]
  df$day <-
    ifelse(seq_len(nrow(df)) %in% morning_indices, "morning", "afternoon")
  df$day_length <- day_length
  hemi <- ifelse(grepl("lh", df_name), "lh", "rh")
  df$hemi <- hemi
  df$time_seq <- seq_len(nrow(df))
  return(df)
}

# Use mapply to change all the data in the datlist
data_sets <-
  mapply(assign_data, 
         data_sets,
         morning_indices_list,
         day_length,
         names(data_sets),
         SIMPLIFY = FALSE)

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
    "cjv",
    "efc",
    "inu_med"
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

data_sets$sub1_lh <- cbind(data_sets$sub1_lh, qc_2_control[1:38,])
data_sets$sub1_rh <- cbind(data_sets$sub1_rh, qc_2_control[1:38,])

data_sets$sub2_lh <- cbind(data_sets$sub2_lh, qc_2_control[39:78,])
data_sets$sub2_rh <- cbind(data_sets$sub2_rh, qc_2_control[39:78,])

data_sets$sub3_lh <- cbind(data_sets$sub3_lh, qc_2_control[79:103,])
data_sets$sub3_rh <- cbind(data_sets$sub3_rh, qc_2_control[79:103,])

# ##-------------------- change the time seq to the corresponding month

# ## -----for sub1
# assign_time_seq <- function(data) {
#   data$time_seq[1] <- 1
#   data$time_seq[2:8] <- 2
#   data$time_seq[9:10] <- 3
#   data$time_seq[11:17] <- 4
#   data$time_seq[18:22] <- 5
#   data$time_seq[23:28] <- 11
#   data$time_seq[29:33] <- 12
#   data$time_seq[34:38] <- 14 # feb next year
#   return(data)
# }
#
# data_sets$sub1_lh <- assign_time_seq(data_sets$sub1_lh)
# data_sets$sub1_rh <- assign_time_seq(data_sets$sub1_rh)
#
# ## ---- for sub2
# assign_time_seq <- function(data) {
#   data$time_seq[1] <- 1
#   data$time_seq[2:9] <- 2
#   data$time_seq[10:12] <- 3
#   data$time_seq[13:19] <- 4
#   data$time_seq[20:24] <- 5
#   data$time_seq[25:31] <- 11
#   data$time_seq[32:35] <- 12
#   data$time_seq[36:40] <- 14 # feb next year
#   return(data)
# }
#
# data_sets$sub2_lh <- assign_time_seq(data_sets$sub2_lh)
# data_sets$sub2_rh <- assign_time_seq(data_sets$sub2_rh)
#
# ## ---- for sub3
# assign_time_seq <- function(data) {
#   data$time_seq[1:3] <- 2
#   data$time_seq[4:6] <- 3
#   data$time_seq[7:14] <- 4
#   data$time_seq[15:19] <- 5
#   data$time_seq[20:24] <- 11
#   data$time_seq[25] <- 12
#   return(data)
# }
#
# data_sets$sub3_lh <- assign_time_seq(data_sets$sub3_lh)
# data_sets$sub3_rh <- assign_time_seq(data_sets$sub3_rh)


##---------------------------------------------Part II Plot the area from thinnest to thickest
sub1_long <-
  pivot_longer(
    rbind(data_sets$sub1_lh, data_sets$sub1_rh),
    cols = -c(
      subject,
      time_seq,
      hemi,
      day,
      day_length,
      cjv,
      efc,
      fwhm_avg,
      wm2max,
      fber,
      inu_med,
      cnr,
      snr_total,
      rpve_gm
    ),
    names_to = "region",
    values_to = "area"
  )
sub2_long <-
  pivot_longer(
    rbind(data_sets$sub2_lh, data_sets$sub2_rh),
    cols = -c(
      subject,
      time_seq,
      hemi,
      day,
      day_length,
      cjv,
      efc,
      fwhm_avg,
      wm2max,
      fber,
      inu_med,
      cnr,
      snr_total,
      rpve_gm
    ),
    names_to = "region",
    values_to = "area"
  )
sub3_long <-
  pivot_longer(
    rbind(data_sets$sub3_lh, data_sets$sub3_rh),
    cols = -c(
     subject,
     time_seq,
      hemi,
      day,
     day_length,
      cjv,
      efc,
      fwhm_avg,
      wm2max,
      fber,
      inu_med,
      cnr,
      snr_total,
      rpve_gm
    ),
    names_to = "region",
    values_to = "area"
  )
# Combing the long version data together
data_sets_long <- list(sub1_long, sub2_long, sub3_long)
names(data_sets_long) <- c("sub1_long", "sub2_long", "sub3_long")


# define a function
process_data_and_plot <- function(data) {
  # calcualte the average area of each region
  average_area_region <- data %>%
    group_by(region) %>%
    summarise(avg_area = mean(area, na.rm = TRUE)) %>%
    arrange(avg_area) %>%
    mutate(region = factor(region, levels = region))
  
  # update the region
  data_updated <- data %>%
    filter(region != "total") %>%
    mutate(region = factor(region, levels = average_area_region$region))
  
  #
  average_area <- data_updated %>%
    group_by(region, hemi) %>%
    summarise(avg_area = mean(area, na.rm = TRUE)) %>%
    ungroup() %>%
    arrange(region, hemi)
  
  unique_regions <- unique(average_area$region)
  desired_breaks <-
    c(head(unique_regions, 3), tail(unique_regions, 3))
  
  #
  plot <-
    ggplot(data_updated, aes(x = region, y = area, color = hemi)) +
    geom_point(alpha = 0.3) +
    geom_line(
      data = average_area,
      aes(x = region, y = avg_area, group = hemi),
      size = 1
    ) +
    geom_point(
      data = average_area,
      aes(x = region, y = avg_area, group = hemi),
      size = 3
    ) +
    labs(title = "", x = "", y = "") +
    scale_color_manual(
      values = c("#f66a66", "#14b7bb"),
      labels = c("lh", "rh")
    ) +
    scale_x_discrete(breaks = desired_breaks) +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      #
      panel.grid.minor = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(
        size = 12,
        family = "Arial",
        face = "bold",
        colour = "black"
      ),
      axis.ticks.x = element_line(color = "black"),
      #
      axis.text.x = element_text(
        angle = 45, vjust = 1, hjust = 1,
        size = 12,
        family = "Arial",
        face = "bold",
        colour = "black",
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
}

plots_area_ordered <- lapply(data_sets_long, process_data_and_plot)

save(plots_area_ordered, file = "plots_area_ordered.RData")

### ------------------------------------Part III calculate the absolute percentage change of area

calculate_percentage_change <-
  function(data_set, cols_to_change, cols_to_keep) {
    percentage_change <- data_set[, cols_to_change] %>%
      sweep(
        2,
        colMeans(., na.rm = TRUE),
        FUN = function(x, mean)
          ((x / mean) - 1) * 100
      )
    cbind(percentage_change, data_set[, cols_to_keep])
  }

columns_to_change <- 1:35
columns_to_keep <- 36:39

percentage_change_lh <- do.call(
  rbind,
  mapply(
    calculate_percentage_change,
    data_sets[c("sub1_lh", "sub2_lh", "sub3_lh")],
    MoreArgs = list(cols_to_change = columns_to_change, cols_to_keep = columns_to_keep),
    SIMPLIFY = FALSE
  )
) %>%
  cbind(., qc_2_control)

percentage_change_rh <- do.call(
  rbind,
  mapply(
    calculate_percentage_change,
    data_sets[c("sub1_rh", "sub2_rh", "sub3_rh")],
    MoreArgs = list(cols_to_change = columns_to_change, cols_to_keep = columns_to_keep),
    SIMPLIFY = FALSE
  )
) %>%
  cbind(., qc_2_control)

percentage_change <-
  rbind(percentage_change_lh, percentage_change_rh)

# percentage_change_area <- percentage_change
# save (percentage_change_area, file = "percentage_change_area.RData")
write.csv(percentage_change, "percentage_change_area.csv", row.names = TRUE)

##-------------------------calculate the correlation

percentage_change_abs <- abs(select_if(percentage_change, is.numeric))

brain_regions <- 1:35 # brain regions

# calculate the correlation and convert into long form
# qc_time_vars <- 36:44 # if use the absolute percentgae change, qc matrix and time scanning order
# cor_matrix <- cor(percentage_change_abs[, brain_regions], percentage_change_abs[, qc_time_vars], use = "complete.obs") %>%
#   round(., 2)

qc_time_vars <- 40:48 # if use the normal percentage change
cor_matrix <-
  cor(percentage_change[, brain_regions], percentage_change[, qc_time_vars], use = "complete.obs") %>%
  round(., 2)

long_cor_matrix <- melt(cor_matrix)
plot_cor_matrix <-
  ggplot(data = long_cor_matrix, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile(color = "white") +  #
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-0.6, 0.6),
    space = "Lab"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.text.y = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    )
  ) +
  labs(x = '', y = '', fill = '') +
  coord_fixed()

plot_cor_matrix


#### ---------------------------------------------Part IV linear regression analysis

##---------------converted to long form

percentage_change_long <- pivot_longer(
  percentage_change,
  cols = -c(
    subject,
    hemi,
    day,
    day_length,
    time_seq,
    cjv,
    efc,
    fwhm_avg,
    wm2max,
    fber,
    inu_med,
    cnr,
    snr_total,
    rpve_gm
  ),
  names_to = "region",
  values_to = "area_per_change"
) %>%
  mutate(
    day = as.factor(day),
    hemi = as.factor(hemi),
    region = as.factor(region),
    subject = as.factor(subject),
  )
# percentage_change_area_long <- percentage_change_long
# save(percentage_change_area_long, file="percentage_change_area_long.RData")

write.csv(percentage_change_long, "percentage_change_area_long.csv", row.names = TRUE)

##----------------to exclude the outliers
# percentage_change_long <- percentage_change_long %>%
#   filter(!(time_seq >= 34 & time_seq <= 38 & subject == "sub1")) %>%
#   filter(!(time_seq >= 36 & time_seq <= 40 & subject == "sub2")) %>%
#   filter(!((time_seq == 1 | time_seq == 7) & subject == "sub3"))

percentage_change_long <- percentage_change_long %>%
  # filter(!((time_seq == 16) & subject == "sub1")) %>%
  # filter(!((time_seq == 23 | time_seq == 32) & subject == "sub2")) %>%
  filter(!((time_seq == 1 | time_seq == 7) & subject == "sub3"))


library(lme4)
library(MuMIn)
library(lmerTest)

# modell_qc <- percentage_change_long %>%
#   mutate_at(
#     vars(
#       "cnr",
#       "snr_total",
#       "fber",
#       "wm2max",
#       "rpve_gm",
#       "fwhm_avg",
#       "cjv",
#       "efc",
#       "inu_med"
#     ),
#     list( ~ scale(.))
#   ) %>%
#   lmer(
#     area_per_change ~  cjv +  efc + fwhm_avg + wm2max + fber + inu_med +
#       cnr + snr_total + rpve_gm + (1 | region),
#     data = .
#   )
# 
# r.squaredGLMM(modell_qc)
# 
# p_values <- summary(modell_qc)$coefficients[, "Pr(>|t|)"]
# 
# p_adjusted <- p.adjust(p_values, method = "BH")
# p_adjusted_small <- p_adjusted[p_adjusted < 0.01] %>%
#   print()

# # summary(modell_qc)
# p_adjusted_small <- lapply(modell_qc, function(model) {
#   p_values <- summary(model)$coefficients[, "Pr(>|t|)"]
#   p_adjusted <- p.adjust(p_values, method = "BH")
#   p_adjusted_small <- p_adjusted[p_adjusted < 0.01]
#   print(p_adjusted_small)
# })


## list the data into its own
sub_per_change_long <-
  list(
    filter(percentage_change_long, subject == "sub1"),
    filter(percentage_change_long, subject == "sub2"),
    filter(percentage_change_long, subject == "sub3")
  )

## ---------------------------------------about the time of day effect

# with each subj
modell_day <- lapply(sub_per_change_long, function(df) {
  
  df$fber_log <- log(df$fber)
  
  normalized_df <- df %>%
    mutate_at(
      vars(
        cjv,
        efc,
        fwhm_avg,
        wm2max,
        cnr,
        snr_total,
        rpve_gm,
        inu_med
      ),
      list(~ scale(.))
    )
  lmer(
    area_per_change ~ day * region + cjv + efc + fwhm_avg + wm2max + 
      fber_log+ inu_med + cnr + snr_total + rpve_gm + (1 | hemi),
    data = normalized_df
  )
})

summary(modell_day[[1]])
r.squaredGLMM(modell_day[[1]])

p_adjusted_small <- lapply(modell_day, function(model) {
  p_values <- summary(model)$coefficients[, "Pr(>|t|)"]
  p_adjusted <- p.adjust(p_values, method = "holm")
  p_adjusted_small <- p_adjusted[p_adjusted < 0.01]
  print(p_adjusted_small)
})

## across all subj
# modell_day_all <- percentage_change_long %>%
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
#       rpve_gm,
#       fber_log
#     ),
#     list(~ scale(.))
#   ) %>%
#   lmer(
#     area_per_change ~ day * region + cjv + efc + fwhm_avg + wm2max +
#       fber_log + inu_med + cnr + snr_total + rpve_gm + (1 | hemi),
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

## -----------------------------------about day length

# within each subj
modell_time <- lapply(sub_per_change_long, function(df) {
  df$fber_log <- log(df$fber)
  normalized_df <- df %>%
    mutate_at(
      vars(
        day_length,
        cjv,
        efc,
        fwhm_avg,
        wm2max,
        cnr,
        snr_total,
        rpve_gm,
        inu_med,
        fber_log
      ),
      list(~ scale(.))
    )
  lmer(
    area_per_change ~ day_length * region + cjv + efc + fwhm_avg + wm2max +
      fber_log + inu_med + cnr + snr_total + rpve_gm + (1 | hemi),
    data = normalized_df
  )
})

summary(modell_time[[1]])
r.squaredGLMM(modell_time[[1]])

p_adjusted_small <- lapply(modell_time, function(model) {
  p_values <- summary(model)$coefficients[, "Pr(>|t|)"]
  p_adjusted <- p.adjust(p_values, method = "holm")
  p_adjusted_small <- p_adjusted[p_adjusted < 0.01]
  print(p_adjusted_small)
})

## across all subj
# modell_time_all <- percentage_change_long %>%
#   mutate_at(
#     vars(
#       day_length,
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
#     area_per_change ~ day_length * region + cjv + efc + fwhm_avg + wm2max +
#       fber + inu_med + cnr + snr_total + rpve_gm + (1 | subject),
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
