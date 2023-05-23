

# Set the working directory to the path where your files are located
setwd("/Users/wang/Desktop/Research_projects/BBSC/MRS/Output_v2_4/QuantifyResults")

# Read in the tsv files and define the intersted columns
file_names <- c("A_tCr_Voxel_1_Basis_1.tsv", "A_rawWaterScaled_Voxel_1_Basis_1.tsv", 
                "A_CSFWaterScaled_Voxel_1_Basis_1.tsv", "A_TissCorrWaterScaled_Voxel_1_Basis_1.tsv")
cols <- c("Cr", "PCr", "Cr_PCr", "Glu", "Glu_Gln", "NAA", "NAA_NAAG", "PCh_GPC", "mI")

# sub1's data is from 1:38; sub2 is 39:77; sub3 is 78:102
rows_of_interest <- list(1:38, 39:77, 78:102)

# Create a function to calculate the coefficient of variation
cv <- function(x) { 100 * (sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)) }

# Create the a blank results list to stor the final results
results <- list()

# Iterate over the subset of files
for (i in 1:3) {
  
  for (j in 1:length(file_names)) {
  # Subset the specified rows and the specified columns
  subset_file <- read.delim(file_names[j])
  subset <- subset_file[rows_of_interest[[i]], cols]
  
  # Calculate the coefficient of variation for each subset and column
  cv_subset <- apply(subset, 2, cv)
  
  # Store the results in a list
  results[[i]] <- data.frame(column = c("Cr", "PCr", "tCr", "Glu", "Glx", "NAA", "tNAA", "tCh", "mI"), tCr = cv_subset)
}

# Combine the results into a single data frame
CV_sub <- do.call(cbind, results)

# Specify the directory to save the file in
output_dir <- "/Users/wang/Desktop/Research_projects/BBSC/MRS/analysis/"

# Save the results as sub1_CV.tsv in the output directory
#write.table(CV_sub, file = paste(output_dir, "sub1_CV.tsv", sep = "/"), sep = "\t", row.names = FALSE)

# Print the result
print(CV_sub)
