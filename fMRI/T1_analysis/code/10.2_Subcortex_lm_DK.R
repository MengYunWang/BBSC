# This code is to
# 1. load data and add more factors
# 2. calculate the percentage change
# 3. linear regression, meta, and equivalent

# WANG. 04-April-2024
#updated 26-Aug-2024

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

#### -------------------------------------------------------Part I Load the data and add factors----------------------------
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
    dplyr::select(-c(33:39)) %>%
    dplyr::select(-c(
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

# data_sets[["sub3"]][c(1,7), ] <- NA # delete the 1st and 7th, pay attention to the files that need to be saved

# add more factors into the data
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

morning_indices_list <-
  list(morning_sub1, morning_sub2, morning_sub3)

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
day_length <- list(day_length_sub1, day_length_sub2, day_length_sub3) 

# Define a function that insert several columns day, and time_seq
assign_data <- function(df, morning_indices, day_length, df_name) {
  df$day <-
    ifelse(seq_len(nrow(df)) %in% morning_indices, "morning", "afternoon")
  df$day_length <- day_length
  df$time_seq <- seq_len(nrow(df))
  return(df)
}

# Use mapply to change all the data in the datalist
data_sets <-
  mapply(assign_data,
         data_sets,
         morning_indices_list,
         day_length,
         names(data_sets),
         SIMPLIFY = FALSE
  )

# --load the qc data and replace 1 2 3 with 01 02 03
qc_data <- read_tsv(file = "group_T1w_v23.tsv") %>%
  mutate(bids_name = sub("ses-([1-9])([^0-9])", "ses-0\\1\\2", bids_name)) %>%
  arrange(bids_name)

qc_indx <-
  c(
    "cjv",
    "efc",
    "fwhm_avg",
    "cnr",
    "snr_gm",
    "snr_csf",
    "wm2max",
    "rpve_gm",
    "rpve_csf",
    "fber",
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

data_sets$sub1 <- cbind(data_sets$sub1, qc_2_control[1:38, ])

data_sets$sub2 <- cbind(data_sets$sub2, qc_2_control[39:78, ])

data_sets$sub3 <- cbind(data_sets$sub3, qc_2_control[79:103, ])



data_sets_long <- lapply (data_sets, function(df) {
  pivot_longer(
    df,
    cols = -c(
      subject,
      time_seq,
      day,
      day_length,
      cjv,
      efc,
      fwhm_avg,
      cnr,
      snr_gm,
      snr_csf,
      wm2max,
      fber,
      inu_med,
      rpve_gm,
      rpve_csf
    ),
    names_to = "region",
    values_to = "subcortex"
  )
})


#----------------------------------------------Part II calculate the percentage change and converted into long form-----------------

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

percentage_change_sub1 <- calculate_percentage_change(data_sets$sub1, 1:35, 36:50)
percentage_change_sub2 <- calculate_percentage_change(data_sets$sub2, 1:35, 36:50)
percentage_change_sub3 <- calculate_percentage_change(data_sets$sub3, 1:35, 36:50)
percentage_change <- rbind(percentage_change_sub1, percentage_change_sub2, percentage_change_sub3)

#-------
sub_per_change_long <- lapply(names(data_sets), function(name) {
  df <- data_sets[[name]]
  
  # Calculate the percentage change for all columns except the ones to exclude
  percentage <- calculate_percentage_change(df, 1:35, 36:50)
  
  # Convert the wide data to long format
  percentage_long <-
    pivot_longer(
      percentage,
      cols = !c(subject,
                day,
                day_length,
                time_seq,
                cjv,
                efc,
                fwhm_avg,
                fber,
                inu_med,
                cnr,
                wm2max,
                snr_gm,
                rpve_gm,
                snr_csf,
                rpve_csf),
      names_to = "region",
      values_to = "subcortex_per_change"
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
    mutate(
      day = as.factor(day),
      hemi = as.factor(hemi),
      region = as.factor(region),
      subject = as.factor(subject),
    )
  
  return(percentage_long)
})

names(sub_per_change_long) <- c("sub1_long", "sub2_long", "sub3_long")

# combine them together into one
percentage_change_long <- rbind(sub_per_change_long$sub1_long, sub_per_change_long$sub2_long, sub_per_change_long$sub3_long)
save(percentage_change_long, file = "percentage_change_subcortex_long.RData")

## -------------------------------------------------------------Part III linear regression------------------------

## ------------------------------------about the time of day effect

# 1. lm modeling the time of day effect within subjects
model_day_subcortex <- lapply(data_sets_long, function(df) {
  
  df <- filter(df, region == "total") # only look at the total subcortex
  
  # Standardize the variables
  df$subcortex_normalized <- scale(df$subcortex)
  
  model <- lm(subcortex_normalized ~ day, data = df)
  
  # extract the beta se and p values
  beta_values <- summary(model)$coefficients[, "Estimate"]
  standard_errors <- summary(model)$coefficients[, "Std. Error"]
  p_value <- summary(model)$coefficients[, "Pr(>|t|)"]
  
  return(list("model" = model, "p_value" = p_value,
              "beta" = beta_values, "SE" = standard_errors
  )
  )
})

# 2. do the meta analysis based on three sub using metafor
library(metafor)
data2meta_day <- data.frame(
  beta = sapply(model_day_subcortex, function(x) x$beta["daymorning"]),
  se_beta = sapply(model_day_subcortex, function(x) x$SE["daymorning"])
)

meta_model_day_subcortex <- rma(yi=beta, sei=se_beta, data = data2meta_day)
summary(meta_model_day_subcortex)
save(meta_model_day_subcortex, file="meta_model_day_subcortex.RData")

# 3. do the equivalence testing based on the meta results
library(TOSTER)
TOSTmeta(ES = summary(meta_model_day_subcortex)$beta, se = summary(meta_model_day_subcortex)$se,
         low_eqbound_d=-0.5, high_eqbound_d=0.5, alpha=0.05)


## -----------------------------------about day_length effect

# 1. lm modelling within each subj
model_time_subcortex <- lapply(data_sets_long, function(df) {
  
  df <- filter(df, region == "total")
  
  # Standardize the variables
  df$subcortex_normalized <- scale(df$subcortex)
  df$day_length <- scale(df$day_length)
  
  model <- lm(subcortex_normalized ~ day_length, data = df)
  
  beta_values <- summary(model)$coefficients[, "Estimate"]
  standard_errors <- summary(model)$coefficients[, "Std. Error"]
  p_value <- summary(model)$coefficients[, "Pr(>|t|)"]
  
  
  return(list("model" = model, "p_value" = p_value,
              "beta" = beta_values, "SE" = standard_errors
  )
  )
})

# 2. do the meta analysis based on three sub using metafor
data2meta_time <- data.frame(
  beta = sapply(model_time_subcortex, function(x) x$beta["day_length"]),
  se_beta = sapply(model_time_subcortex, function(x) x$SE["day_length"])
)

meta_model_time_subcortex <- rma(yi=beta, sei=se_beta, data = data2meta_time)
summary(meta_model_time_subcortex)
save(meta_model_time_subcortex, file="meta_model_time_subcortex.RData")

# 3. equivalence testing based on the meta results
TOST_result_time <- TOSTmeta(ES = summary(meta_model_time_subcortex)$beta, se = summary(meta_model_time_subcortex)$se,
         low_eqbound_d=-0.5, high_eqbound_d=0.5, alpha=0.05)

  

  
  
  
###### obsolete code from before
  # #--------------------------------
  # brain_regions <- 1:35 # brain regions
  # qc_time_vars <- 39:49 # other variables
  # 
  # cor_matrix <-
  #   cor(percentage_change[, brain_regions], percentage_change[, qc_time_vars], use = "complete.obs") %>%
  #   round(., 2)
  # 
  # # qc_cor_matrix <-
  # #   cor(percentage_change[, qc_time_vars], percentage_change[, qc_time_vars], use = "complete.obs") %>%
  # #   round(., 2)
  # 
  # long_cor_matrix <- melt(cor_matrix)
  # plot_cor_matrix <-
  #   ggplot(data = long_cor_matrix, aes(x = Var2, y = Var1, fill = value)) +
  #   geom_tile(color = "white") +  #
  #   scale_fill_gradient2(
  #     low = "blue",
  #     high = "red",
  #     mid = "white",
  #     midpoint = 0,
  #     limit = c(-1, 1),
  #     space = "Lab"
  #   ) +
  #   theme_minimal() +
  #   theme(
  #     axis.text.x = element_text(
  #       angle = 45,
  #       hjust = 1,
  #       vjust = 1
  #     ),
  #     axis.text.y = element_text(
  #       angle = 45,
  #       hjust = 1,
  #       vjust = 1
  #     )
  #   ) +
  #   labs(x = '', y = '', fill = '') +
  #   coord_fixed()
  # 
  # plot_cor_matrix
  # 
  # library(lme4)
  # library(MuMIn)
  # library(lmerTest)
  # 
  # 
  # qc_2_plot <-
  #   c(
  #     "cnr",
  #     "snr_gm",
  #     "fber",
  #     "cjv",
  #     "efc",
  #     "inu_med",
  #     "wm2max",
  #     "fwhm_avg",
  #     "rpve_gm"
  #   )
  # 
  # subject_colors <-
  #   c(
  #     "sub1" = "#D40000",
  #     "sub2" = "#377EB8",
  #     "sub3" = "#4DAF4A"
  #   )
  # 
  # subject_grayscale <- c(
  #   "sub1" = "#D3D3D3",  # Light gray
  #   "sub2" = "#A9A9A9", # Medium gray
  #   "sub3" = "#808080" # Dark gray
  # )
  # 
  # 
  # generate_plot <- function(data, variable_name, subject_colors, use_grayscale = FALSE) {
  #   p <- ggplot(
  #     data,
  #     aes_string(
  #       x = variable_name,
  #       y = "percent_change",
  #       color = "subject",
  #       group = "subject"
  #     )
  #   ) +
  #     geom_point(size = 0.5, alpha = 0.5) +
  #     labs(
  #       title = variable_name,
  #       x = "",
  #       y = ""
  #     ) +
  #     theme_minimal() +
  #     theme(
  #       panel.grid.major = element_blank(),
  #       panel.grid.minor = element_blank(),
  #       legend.position = "none",
  #       plot.title = element_text(
  #         hjust = 0.5,
  #         size = 16,
  #         family = "Arial",
  #         face = "bold",
  #         color = "black"
  #       ),
  #       axis.ticks.x = element_line(color = "black"),
  #       axis.text.x = element_text(
  #         size = 14,
  #         family = "Arial",
  #         face = "bold",
  #         color = "black"
  #       ),
  #       axis.line.x = element_line(color = "black"),
  #       axis.ticks.y = element_blank(),
  #       axis.text.y = element_text(
  #         size = 14,
  #         family = "Arial",
  #         face = "bold",
  #         color = "black"
  #       ),
  #       axis.line.y = element_line(color = "black"),
  #       axis.title.y = element_blank()
  #     ) +
  #     geom_hline(
  #       yintercept = c(0),
  #       linetype = "solid",
  #       color = "grey"
  #     ) +
  #     geom_hline(
  #       yintercept = c(-10, 10),
  #       linetype = "dashed",
  #       color = "green"
  #     ) +
  #     geom_hline(
  #       yintercept = c(-15, 15),
  #       linetype = "dashed",
  #       color = "red"
  #     ) +
  #     scale_y_continuous(limits = c(-20, 20), breaks = c(-15, 0, 15))
  #   
  #   # Choose color scale based on use_grayscale flag
  #   if(use_grayscale) {
  #     p <- p + scale_color_manual(values = subject_grayscale)
  #   } else {
  #     p <- p + scale_color_manual(values = subject_colors)
  #   }
  #   
  #   return(p)
  # }
  # 
  # 
  # plots_per_hemisphere <- lapply(qc_2_plot, function(var) {
  #   generate_plot(percentage_change_long, var, subject_colors, use_grayscale = var %in% c(" "))
  # })
  # 
  # 
  # plots_per_hemisphere[[1]] <- plots_per_hemisphere[[1]] + 
  #   theme (plot.margin = margin(t = 30, r = 5, b = 0, l = 5, unit = "pt")) +
  #   scale_x_continuous(breaks = c(1.5, 2.0), labels = c("1.5", "2.0"))
  # 
  # plots_per_hemisphere[[2]] <- plots_per_hemisphere[[2]] + 
  #   scale_x_continuous(breaks = c(6.0, 7.0), labels = c("6.0", "7.0"))
  # 
  # plots_per_hemisphere[[3]] <- plots_per_hemisphere[[3]] + 
  #   scale_x_continuous(trans = 'log10', 
  #                      breaks = c(300, 3000, 30000), labels = c("300","3000","30000"))
  # 
  # plots_per_hemisphere[[6]] <- plots_per_hemisphere[[6]] + 
  #   scale_x_continuous(breaks = c(0.4, 0.5), labels = c("0.4", "0.5"))
  # # 
  # # plots_per_hemisphere[[7]] <- plots_per_hemisphere[[7]] + 
  # #   scale_x_continuous(breaks = c(0.001, 0.003), labels = c("0.001", "0.003"))
  # 
  # plots_per_hemisphere[[8]] <- plots_per_hemisphere[[8]] + 
  #   scale_x_continuous(breaks = c(4.0, 4.2, 4.4), labels = c("4.0", "4.2", "4.4"))
  # 
  # plots_per_hemisphere[[9]] <- plots_per_hemisphere[[9]] + 
  #   scale_x_continuous(breaks = c(7.25, 7.75,8.25),labels = c("7.25", "7.75","8.25"))
  # 
  # 
  # plots_subcortex <- wrap_plots(plots_per_hemisphere, ncol = 3)
  # 
  # # save(plots_subcortex, file = "plots_subcortex.RData")
  # 
  # # with each subj
  # modell_day <- lapply(sub_per_change_long, function(df) {
  #   normalized_df <- df %>%
  #     mutate_at(
  #       vars(
  #         cjv,
  #         efc,
  #         fwhm_avg,
  #         fber,
  #         cnr,
  #         snr_gm,
  #         rpve_gm,
  #         wm2max,
  #         inu_med
  #       ),
  #       list(~ scale(.))
  #     )
  #   lmer(
  #     percent_change ~ day * region + cjv + efc + fwhm_avg + fber + wm2max + inu_med +
  #       cnr + snr_gm + rpve_gm + (1 | hemi),
  #     data = normalized_df
  #   )
  # })
  # 
  # summary(modell_day[[1]])
  # r.squaredGLMM(modell_day[[2]])
  # 
  # p_adjusted_small <- lapply(modell_day, function(model) {
  #   p_values <- summary(model)$coefficients[, "Pr(>|t|)"]
  #   p_adjusted <- p.adjust(p_values, method = "BH")
  #   p_adjusted_small <- p_adjusted[p_adjusted < 0.01]
  #   print(p_adjusted_small)
  # })
  # 
  # ## across all subj
  # modell_day_all <- percentage_change_long %>%
  #   mutate_at(
  #     vars(
  #       cjv,
  #       efc,
  #       fwhm_avg,
  #       fber,
  #       wm2max,
  #       inu_med,
  #       cnr,
  #       snr_gm,
  #       rpve_gm
  #     ),
  #     list(~ scale(.))
  #   ) %>%
  #   lmer(
  #     percent_change ~ day * region + cjv + efc + fwhm_avg + fber + wm2max + inu_med +
  #       cnr + snr_gm + rpve_gm + (1 | subject),
  #     data = .
  #   )
  # 
  # r.squaredGLMM(modell_day_all)
  # 
  # p_values <- summary(modell_day_all)$coefficients[, "Pr(>|t|)"]
  # 
  # p_adjusted <- p.adjust(p_values, method = "BH")
  # p_adjusted_small <- p_adjusted[p_adjusted < 0.01] %>%
  #   print()
  # 
  # ## -----------------------------------about day_length
  # 
  # # within each subj
  # modell_time <- lapply(sub_per_change_long, function(df) {
  #   normalized_df <- df %>%
  #     mutate_at(
  #       vars(
  #         day_length,
  #         cjv,
  #         efc,
  #         fwhm_avg,
  #         fber,
  #         cnr,
  #         snr_gm,
  #         rpve_gm,
  #         wm2max,
  #         inu_med
  #       ),
  #       list(~ scale(.))
  #     )
  #   lmer(
  #     percent_change ~  day_length * region + cjv + efc + fwhm_avg + fber +
  #       wm2max + inu_med + cnr + snr_gm + rpve_gm + (1 | hemi),
  #     data = normalized_df
  #   )
  # })
  # 
  # summary(modell_time[[2]])
  # r.squaredGLMM(modell_time[[2]])
  # 
  # p_adjusted_small <- lapply(modell_time, function(model) {
  #   p_values <- summary(model)$coefficients[, "Pr(>|t|)"]
  #   p_adjusted <- p.adjust(p_values, method = "BH")
  #   p_adjusted_small <- p_adjusted[p_adjusted < 0.01]
  #   print(p_adjusted_small)
  # })
  # 
  # ## across all subj
  # modell_time_all <- percentage_change_long %>%
  #   mutate_at(
  #     vars(
  #       day_length,
  #       cjv,
  #       efc,
  #       fwhm_avg,
  #       fber,
  #       wm2max,
  #       inu_med,
  #       cnr,
  #       snr_gm,
  #       rpve_gm
  #     ),
  #     list(~ scale(.))
  #   ) %>%
  #   lmer(
  #     percent_change ~ day_length * region + cjv + efc + fwhm_avg + fber +
  #       wm2max + inu_med + cnr + snr_gm + rpve_gm + (1 | subject),
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
  # ##------------------ Look at the CSF values and its own quality
  # 
  # cor.test(percentage_change[, "CSF"], percentage_change[, "snr_csf"]) 
  # 
  # cor.test(percentage_change[, "CSF"], percentage_change[, "rpve_csf"]) 
  # 
  # 
  # 
  # ggplot(percentage_change, aes(x = snr_csf, y = CSF)) +
  #   geom_point(color = "blue") +  # Adds scatter plot points
  #   theme_minimal() +  # Uses a minimal theme for the plot
  #   labs(title = "Scatter Plot of CSF vs SNR_CSF",
  #        x = "SNR_CSF",
  #        y = "CSF")