# This code is to explore the relationship between three phenotypes and 
# plot Fig. 6 and S12

# WANG. 03-April-2024

# Remove all objects created before to prevent clash
rm(list = ls())

library(dplyr)
library(ggplot2)
library(reshape2)
library(patchwork)
library(RColorBrewer) 

# Set the working directory to the path where your files are located
setwd("/Users/wang/Desktop/Research_projects/BBSC/Functional/fmri/derivatives/freesurfer-v7.2.0/")
# setwd("//Users/joeywang/Desktop/BBSC/Functional/fmri/derivatives/freesurfer-v7.2.0/")

#----------------------------------------------------------Load the data

# generated from 7.2
percentage_change_thickness <- read.csv("percentage_change_thickness_qc.csv", row.names = 1)
percentage_change_thickness_long <- read.csv("percentage_change_thickness_long_qc.csv", row.names = 1)

# generated from 8.2
percentage_change_area <- read.csv("percentage_change_area_qc.csv", row.names = 1)
percentage_change_area_long <- read.csv("percentage_change_area_long_qc.csv", row.names = 1)

# generated from 9.2
percentage_change_volume <- read.csv("percentage_change_volume_qc.csv", row.names = 1)
percentage_change_volume_long <- read.csv("percentage_change_volume_long_qc.csv", row.names = 1)


#------------------------------------------Plot correlation matrix within and between phenotypes
percentage_change <- list (percentage_change_thickness, percentage_change_area, percentage_change_volume)

#-------to exclude the outliers
percentage_change <- lapply(percentage_change, function (data) {
  data_qc <- data %>%
    filter(!((time_seq == 1 | time_seq == 7) & subject == "sub3"))
  } )


# brain_regions <- 1:34 # brain regions
# qc_time_vars <- 38:46 # co-variance


# Function to create a heatmap using ComplexHeatmap and return it
create_ggplot_heatmap <- function(corr_matrix, title = "") {
  
  # Perform hierarchical clustering
  hc <- hclust(dist(corr_matrix))
  hc_order <- hc$order
  
  # Reorder the correlation matrix based on the hierarchical clustering
  corr_matrix_ordered <- corr_matrix[hc_order, hc_order]
  
  # Melt the correlation matrix into a long format
  long_corr_matrix <- melt(corr_matrix_ordered)
  
  colors <- rev(colorRampPalette(brewer.pal(n = 11, name = "RdBu"))(100))
  
  # Create the heatmap
  heatmap_plot <- ggplot(long_corr_matrix, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile() +
    scale_fill_gradientn(colors = colors, limits = c(-1, 1), breaks = c(-0.5, 0, 0.5), labels = c("-0.5", "0", "0.5"), name = "") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14,
                                     family = "Arial",
                                     face = "bold",
                                     color = "black"),
          axis.text.y = element_text(size = 14,
                                     family = "Arial",
                                     face = "bold",
                                     color = "black"),
          axis.title = element_blank(),
          legend.text = element_text(size = 20, family = "Arial", face = "bold", colour = "black")) +
    ggtitle(title)+
    coord_fixed()
  
  return(heatmap_plot)
}

# Calculate the correlation matrix within each phenotypes
corr_within <- lapply(percentage_change, function(data) {
  data_sub1 <- data[grep("sub1", rownames(data)), ]
  data_sub2 <- data[grep("sub2", rownames(data)), ]
  data_sub3 <- data[grep("sub3", rownames(data)), ]
  
  # Calculate  correlation for sub1
  corr_sub1 <- cor(data_sub1[, 1:34])
  
  # Calculate  correlation for sub2
  corr_sub2 <- cor(data_sub2[, 1:34])
  
  # Calculate  correlation for sub2
  corr_sub3 <- cor(data_sub3[, 1:34])
  
  # Calculate  correlation for whole data
  corr_whole <- cor(data[, 1:34])
  
  corr_within <- list(corr_sub1, corr_sub2, corr_sub3, corr_whole)
  
  return(corr_within)
})  


# Calcualte the correlaiton matrix between each phenotypes
corr_thick_area <-  cor(na.omit(percentage_change_thickness[, 1:34]), na.omit(percentage_change_area[, 1:34]))
corr_thick_volume <-  cor(na.omit(percentage_change_thickness[, 1:34]), na.omit(percentage_change_volume[, 1:34]))
corr_area_volume <-  cor(na.omit(percentage_change_area[, 1:34]), na.omit(percentage_change_volume[, 1:34]))

corr_between_diagnal <- cbind (diag(corr_thick_volume), diag(corr_thick_area), diag(corr_area_volume))
##------------------------------------------Summary of correlaiton matrices

corr_list <- list(corr_within[[1]][[4]], corr_within[[2]][[4]],corr_within[[3]][[4]], 
  corr_thick_volume, corr_thick_area, corr_area_volume)

corr_mean <- lapply(corr_list, function(corr_matrix) {
  # find the lower triangle of the correlation matrix and excluding the diagnal values
  corr_matrix_lower_triangle <- corr_matrix[lower.tri(corr_matrix)]
  
  # Apply Fisher's Z-transformation
  Z_transformed_correlations = atanh(corr_matrix_lower_triangle)
  
  # Calculate the mean and standard deviation of the Z-transformed correlations
  mean_Z = mean(Z_transformed_correlations, na.rm = TRUE)
  std_dev_Z = sd(Z_transformed_correlations, na.rm = TRUE)
  
  # Convert the mean Z-score back to a correlation coefficient
  mean_correlation = tanh(mean_Z)
  
  return(mean_correlation)
  
})
##---------------------------------------Distribution of correlation 
# Labels for each matrix
labels <- c("thickness", "area", "volume", "thickness-volume", "thickness-area", "area-volume")

# Function to prepare data: Extract lower triangle, add labels
prepare_data <- function(corr_matrix, label) {
  corr_data <- corr_matrix[lower.tri(corr_matrix)]
  data.frame(value = corr_data, label = label)
}

# Apply the function to each matrix with appropriate labels
corr_data_frames <- lapply(1:length(corr_list), function(i) {
  prepare_data(corr_list[[i]], labels[i])
})

# Combine all data frames into one
combined_corr_data <- do.call(rbind, corr_data_frames)

plot_corr_distribution <-
  ggplot(combined_corr_data, aes(x = value, fill = label)) +
  geom_density(alpha = 0.6) +  # Adjust transparency
  scale_fill_manual(values = c("thickness" = "#FF6666", 
                               "area" = "#FFCC66", 
                               "volume" = "#66CC66", 
                               "thickness-volume" = "#66CCCC", 
                               "thickness-area" = "#6666CC", 
                               "area-volume" = "#CC66CC")) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t = 20, r = 0, b = 20, l = 30, unit = "pt"),
    legend.title = element_blank(),
    legend.position = "right",
    legend.text = element_text(
      size = 12,
      family = "Arial",
      face = "bold",
      colour = "black"
    ),
    axis.ticks.x = element_line(color = "black"),
    axis.text.x = element_text(
      size = 16,
      family = "Arial",
      face = "bold",
      colour = "black",
    ),
    axis.title.x = element_text(
      size = 16,
      family = "Arial",
      face = "bold",
      colour = "black",
    ),
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_blank(),
    axis.title.y = element_text(
      size = 16,
      family = "Arial",
      face = "bold",
      colour = "black",
    ),
    axis.ticks.y = element_line(color = "black"),
    axis.text.y = element_text(
      size = 16,
      family = "Arial",
      face = "bold",
      colour = "black",
    )
    )+
  labs(title = "",
       x = "Correlation Coefficient", 
       y = "Density", fill = "Phenotype") 


FigS12 <- plot_corr_distribution

ggsave("/Users/wang/Library/CloudStorage/OneDrive-UniversityofBergen/Desktop/UiB/Manuscripts/BBSC/4.T1w_Image/Figures/FigS12.png",
       plot = FigS12, width = 9, height = 6.5, units = 'in', dpi = 300)

# ggsave("/Users/joeywang/Library/CloudStorage/OneDrive-UniversityofBergen/Desktop/UiB/Manuscripts/BBSC/4.T1w_Image/Figures/FigS12.png",
#        plot = FigS12, width = 9, height = 6.5, units = 'in', dpi = 300)

##-----------------------community detection



##------------------------------------plot
# Plot heatmaps within each phenotyes
heatmap_within_thickness <- create_ggplot_heatmap(corr_within[[1]][[4]], "")
heatmap_within_area <- create_ggplot_heatmap(corr_within[[2]][[4]], "")
heatmap_within_volume <- create_ggplot_heatmap(corr_within[[3]][[4]], "")

heatmaps_within <- list(heatmap_within_thickness, heatmap_within_area, heatmap_within_volume)

# Plot the correlaiton matrix between each phenotypes
heatmap_thick_area <- create_ggplot_heatmap(corr_thick_area, "")
heatmap_thick_volume <- create_ggplot_heatmap(corr_thick_volume, "")
heatmap_area_volume <- create_ggplot_heatmap(corr_area_volume, "")

heatmaps_between <- list(heatmap_thick_area, heatmap_thick_volume, heatmap_area_volume)


# To plot the combined heatmaps for the first dataset
# design <-
#   "AD
#    BE
#    CF"

Fig6 <- heatmaps_within[[1]] + heatmaps_within[[2]] + 
  theme(plot.margin = margin(t = 60, r = 30, b = 30, l = 30, unit = "pt")) + heatmaps_within[[3]] +
  heatmaps_between[[2]] + heatmaps_between[[1]] +
  theme(plot.margin = margin(t = 30, r = 30, b = 30, l = 30, unit = "pt")) + heatmaps_between[[3]] +
  plot_layout(ncol = 3, guides = 'collect') +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 50,
                                family = "Arial",
                                face = "bold",
                                colour = "black"))


ggsave("/Users/wang/Library/CloudStorage/OneDrive-UniversityofBergen/Desktop/UiB/Manuscripts/BBSC/4.T1w_Image/Figures/Fig6.png",
plot = Fig6, width = 28, height = 18, units = 'in', dpi = 300)

# ggsave("/Users/joeywang/Library/CloudStorage/OneDrive-UniversityofBergen/Desktop/UiB/Manuscripts/BBSC/4.T1w_Image/Figures/Fig6.png",
#        plot = Fig6, width = 27, height = 18, units = 'in', dpi = 300)


