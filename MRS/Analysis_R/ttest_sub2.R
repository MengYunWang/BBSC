
# This code is to calculate time of day effect on the MRS data for sub2

# WANG. 25-April-2023

# Remove all objects created before to prevent clusing


# Set the working directory to the path where your files are located
setwd("/Users/wang/Desktop/Research_projects/BBSC/MRS/Output_v2_4/QuantifyResults")


# create a vector of row indices for morning and afternoon
morning <- c (2:4, 6, 7, 10:13, 16, 18, 21, 22, 24, 25, 27, 29, 31, 33, 34, 36, 37, 39) + 37
afternoon <- c(5, 8, 9, 14, 15, 17, 19, 20, 23, 26, 28, 30, 32, 35, 38, 40) + 37
# the following metabolites were analysized
mrs_of_interest <- c("NAA_NAAG", "NAA", "Cr_PCr", "PCr", "Cr", "Glu_Gln", "Glu", "mI", "PCh_GPC")


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
stats_water_tc_sub2 <- data.frame(mean_morning = mean_morning,
                             std_morning = std_morning,
                             mean_afternoon = mean_afternoon,
                             std_afternoon = std_afternoon)

# initialize a list to store the t-test results for each column
ttest_results_sub2 <- list()

# loop over each column of the datasets and perform a t-test
for (i in 1:ncol(MRS_water_tc_morning)) {
  colname <- paste0(mrs_of_interest[i])
  ttest_results_sub2[[colname]] <- t.test(MRS_water_tc_morning[, i], MRS_water_tc_afternoon[, i])
}

# extract the p-value
p_values_sub2 <- sapply(ttest_results_sub2, function(x) x$p.value)

# adjusted_p_values <- p.adjust(p_values_sub2, method = "fdr")

# print the corrected p-values
# print(adjusted_p_values)




