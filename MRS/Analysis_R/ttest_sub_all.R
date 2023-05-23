
# This code is to calculate paired ttests of the MRS data

# WANG. 25-April-2023

# Remove all objects created before to prevent clusing
rm(list=ls())

# Set the working directory to the path where your files are located
setwd("/Users/wang/Desktop/Research_projects/BBSC/MRS/Output_v2_4/QuantifyResults")


# create a vector of row indices for morning and afternoon
morning_sub1 <- c(1:4, 6, 8, 9, 10, 12, 15, 17, 18, 20, 23, 24, 27, 30, 31, 32, 34, 35, 37)
afternoon_sub1 <- c(5, 7, 11, 13, 14, 16, 19, 21, 22, 25, 26, 28, 29, 33, 36, 38)

morning_sub2 <- c (2:4, 6, 7, 10:13, 16, 18, 21, 22, 24, 25, 27, 29, 31, 33, 34, 36, 37, 39) + 37
afternoon_sub2 <- c(5, 8, 9, 14, 15, 17, 19, 20, 23, 26, 28, 30, 32, 35, 38, 40) + 37

morning_sub3 <- c(1:6, 9, 11, 13, 15, 17, 18, 20, 22, 24) + 77
afternoon_sub3 <- c(7, 8, 10, 12, 14, 16, 19, 21, 23, 25) + 77

morning <- c(morning_sub1, morning_sub2, morning_sub3)
afternoon <- c(afternoon_sub1, afternoon_sub2, afternoon_sub3)

# the following metabolites were analysized
mrs_of_interest <- c("Cr", "PCr", "Cr_PCr", "Glu", "Glu_Gln", "NAA", "NAA_NAAG", "PCh_GPC", "mI")


# read the data
MRS_tCr <- read.delim("A_tCr_Voxel_1_Basis_1.tsv")
MRS_water <- read.delim("A_rawWaterScaled_Voxel_1_Basis_1.tsv")
MRS_water_csf <- read.delim("A_CSFWaterScaled_Voxel_1_Basis_1.tsv")
MRS_water_tc <- read.delim("A_TissCorrWaterScaled_Voxel_1_Basis_1.tsv")

# subset the data based on the row indices
MRS_tCr_morning <- MRS_tCr[morning, mrs_of_interest]
MRS_tCr_afternoon <- MRS_tCr[afternoon,mrs_of_interest]

MRS_water_morning <- MRS_water[morning, mrs_of_interest]
MRS_water_afternoon <- MRS_water[afternoon,mrs_of_interest]

MRS_water_csf_morning <- MRS_water_csf[morning, mrs_of_interest]
MRS_water_csf_afternoon <- MRS_water_csf[afternoon,mrs_of_interest]

MRS_water_tc_morning <- MRS_water_tc[morning, mrs_of_interest]
MRS_water_tc_afternoon <- MRS_water_tc[afternoon,mrs_of_interest]

# calculate the mean and std for the morning and afternoon 
mean_morning <- apply(MRS_water_tc_morning, 2, mean)
mean_afternoon <- apply(MRS_water_tc_afternoon, 2, mean)

std_morning <- apply(MRS_water_tc_morning, 2, sd)
std_afternoon <- apply(MRS_water_tc_afternoon, 2, sd)
stats_water_tc <- data.frame(mean_moring = mean_morning,
                            std_morning = std_morning,
                            mean_afternoon = mean_afternoon,
                            std_afternoon = std_afternoon)

# initialize a list to store the t-test results for each column
ttest_results <- list()

# loop over each column of the datasets and perform a t-test
for (i in 1:ncol(MRS_water_tc_morning)) {
  colname <- paste0("Column", i)
  ttest_results[[colname]] <- t.test(MRS_water_tc_morning[, i], MRS_water_tc_afternoon[, i])
}


# correct for multiple comparisons using FDR method
p_values <- sapply(ttest_results, function(x) x$p.value)
adjusted_p_values <- p.adjust(p_values, method = "fdr")

# print the corrected p-values
print(adjusted_p_values)

#------------------------------------------------------------------------------- not finished yet
## the folloiwng is trying to plot the data with rain cloud style
library(ggplot2)
library(ggridges)
library(gridExtra)

# Define the colors to use in the plot
colors <- c("#FF7F00", "#E41A1C", "#4DAF4A", "#984EA3", "#377EB8", "#FFFF33", "#A65628", "#F781BF", "#999999")

# Create a list of variable names
var_names <- c("Cr", "PCr", "Cr_PCr", "Glu", "Glu_Gln", "NAA", "NAA_NAAG", "PCh_GPC", "mI")
#mrs_of_interest <- c("Cr", "PCr", "Cr_PCr", "Glu", "Glu_Gln", "NAA", "NAA_NAAG", "PCh_GPC", "mI")

# Create a list to store the plots
p_list <- list()

# Create a raincloud plot for each variable
for (i in seq_along (var_names)) {
  p_list[[i]] <- ggplot(MRS_water_tc_morning, aes(x = eval(parse(text = var_names[i])), y = factor(var_names[i]), fill = factor(var_names[i]))) +
    geom_density_ridges_gradient(scale = 2, rel_min_height = 0.01, gradient_lwd = 0.2, gradient_length = unit(0.4, "inches"), gradient_color = colors[i]) +
    geom_point(aes(y = factor(var_names[i]), color = factor(var_names[i])), alpha = 0.5, size = 1) +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    theme_ridges() +
    labs(x = "", y = "") +
    ggtitle(var_names[i]) +
    theme(plot.title = element_text(hjust = 0.5))
}
# Combine all the plots into one using the cowplot package

grid.arrange(
  p_list[[1]], p_list[[2]], p_list[[3]],
  p_list[[4]], p_list[[5]], p_list[[6]],
  p_list[[7]], p_list[[8]], p_list[[9]],
  nrow = 3, ncol = 3
)

# Optional: save the plot as a PDF file
ggsave("raincloud_plot.pdf", width = 8, height = 11, units = "in")






