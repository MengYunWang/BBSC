
# This code is to calculate time of day effect on the MRS data for sub3

# WANG. 25-April-2023

# Remove all objects created before to prevent clusing


# Set the working directory to the path where your files are located
setwd("/Users/wang/Desktop/Research_projects/BBSC/MRS/Output_v2_4/QuantifyResults")


# create a vector of row indices for morning and afternoon
morning <- c(1:6, 9, 11, 13, 15, 17, 18, 20, 22, 24) + 77
afternoon <- c(7, 8, 10, 12, 14, 16, 19, 21, 23, 25) + 77
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
stats_water_tc_sub3 <- data.frame(mean_morning = mean_morning,
                             std_morning = std_morning,
                             mean_afternoon = mean_afternoon,
                             std_afternoon = std_afternoon)

# initialize a list to store the t-test results for each column
ttest_results_sub3 <- list()

# loop over each column of the datasets and perform a t-test
for (i in 1:ncol(MRS_water_tc_morning)) {
  colname <- paste0(mrs_of_interest[i])
  ttest_results_sub3[[colname]] <- t.test(MRS_water_tc_morning[, i], MRS_water_tc_afternoon[, i])
}

# correct for multiple comparisons using FDR method
p_values_sub3 <- sapply(ttest_results_sub3, function(x) x$p.value)

#adjusted_p_values <- p.adjust(p_values_sub3, method = "fdr")

p_values <- c (p_values_sub1, p_values_sub2, p_values_sub3)
adjusted_p_values <- p.adjust(p_values, method = "fdr")
# print the corrected p-values
print(adjusted_p_values)




