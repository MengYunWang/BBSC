
# This code is to calculate the CV of the MRS data

# WANG. 24-April-2023

# Remove all objects created before to prevent clusing
rm(list=ls())

# Set the working directory to the path where your files are located
setwd("/Users/wang/Desktop/Research_projects/BBSC/MRS/Output_v2_4/QuantifyResults")

# Read in the tsv files
file1 <- read.delim("A_tCr_Voxel_1_Basis_1.tsv")
file2 <- read.delim("A_rawWaterScaled_Voxel_1_Basis_1.tsv")
file3 <- read.delim("A_CSFWaterScaled_Voxel_1_Basis_1.tsv")
file4 <- read.delim("A_TissCorrWaterScaled_Voxel_1_Basis_1.tsv")

# Create a function to calculate the coefficient of variation
cv <- function(x) { 100 * (sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)) }

# Subset the first 1:38 rows and the specified columns
cols <- c("NAA_NAAG", "NAA", "Cr_PCr", "PCr", "Cr", "Glu_Gln", "Glu", "mI", "PCh_GPC")
subset1 <- file1[1:38, cols]
subset2 <- file2[1:38, cols]
subset3 <- file3[1:38, cols]
subset4 <- file4[1:38, cols]

# Calculate the coefficient of variation for each subset and column
cv1 <- apply(subset1, 2, cv)
cv2 <- apply(subset2, 2, cv)
cv3 <- apply(subset3, 2, cv)
cv4 <- apply(subset4, 2, cv)

# Store the results in a data frame
CV_sub1 <- data.frame(column = c("NAA_NAAG", "NAA", "Cr_PCr", "PCr", "Cr", "Glu_Gln", "Glu", "mI", "PCh_GPC"), tCr = cv1, water_raw = cv2, water_CSF = cv3, water_tissue_corrected = cv4)

# Specify the directory to save the file in
output_dir <- "/Users/wang/Desktop/Research_projects/BBSC/MRS/analysis/"

# Save the results as sub1_CV.tsv in the output directory
write.table(CV_sub1, file = paste(output_dir, "sub1_CV.tsv", sep = "/"), sep = "\t", row.names = FALSE)

# Print the result
print(CV_sub1)

