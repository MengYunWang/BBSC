# This code is to evaluate the relation between cv and brain size

# WANG. 18-Sep.-2024

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

# Function to read, clean, and scale data
read_clean_scale <- function(file) {
  data <- read.csv(file)[,-1]
  data <- data[-nrow(data),]
  scaled <- scale(data[,2:ncol(data)])
  return(scaled)
}

# Function to generate long version
pivot_long <- function(data, value_name) {
  data %>% 
    data.frame() %>%
    pivot_longer(cols = everything(), names_to = "subj", values_to = value_name)
}

files <- c("thickness_cv.csv", "area_cv.csv", "volume_cv.csv",
           "thickness_mean.csv", "area_mean.csv", "volume_mean.csv")

names <- c("thickness_cv", "area_cv", "volume_cv", 
           "thickness_scaled", "area_scaled", "volume_scaled")

#-----------read and convert data into long version
data_list <- lapply(files, FUN = read_clean_scale)
names(data_list) <- names

long_data <- lapply(seq_along(data_list), function(i) {
  pivot_long(data_list[[i]], names[i])
}) 

#-------combine all the data together
data2analysis <- Reduce(function(dtf1, dtf2) {cbind(dtf1, dtf2[,-1])}, long_data) %>%
  separate(subj, into = c("sub", "hemi"), sep = "_")


#--test and plot
cor_test <- cor.test(data2analysis$volume_scaled, data2analysis$thickness_cv)
print(cor_test)

cor_test <- cor.test(data2analysis$volume_scaled, data2analysis$area_cv)
print(cor_test)

cor_test <- cor.test(data2analysis$volume_scaled, data2analysis$volume_cv)
print(cor_test)


color_mapping <- c("sub1" = "#D40000", "sub2" = "#377EB8", "sub3" = "#4DAF4A")

plot_cv_brain_size <- function(data, cv, name){
  ggplot(data, aes(x = volume_scaled, y = cv, color = sub)) +
  geom_point(size = 2) +
  geom_smooth(method=lm, se=TRUE, color="black", fill="grey") +
  labs(title = name, x = "Scaled Brain Size", y = "Scaled CV") +
  scale_color_manual(values=color_mapping) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 18, family = "Arial", face = "bold", colour = "black", hjust = 0.5),
    axis.title.x = element_text(size = 18, family = "Arial", face = "bold", colour = "black", vjust = -0.5),
    axis.title.y = element_text(size = 18, family = "Arial", face = "bold", colour = "black", vjust = 1.5),
    axis.text.x = element_text(size = 18, family = "Arial", face = "bold", colour = "black"),
    axis.text.y = element_text(size = 18, family = "Arial", face = "bold", colour = "black"),
    axis.ticks.length = unit(0.2, "cm"),
    axis.line = element_line(colour = 'black', linewidth = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank(),
    axis.text.x.top = element_text(angle = 90),
    axis.text.y.right = element_text(angle = 90),
    axis.line.x.top = element_line(color = "black"),
    axis.line.y.right = element_line(color = "black"),
    legend.position = "none"
  )
}

plot_thickness_brain_size <- plot_cv_brain_size(data2analysis, data2analysis$thickness_cv, "Thickness")
plot_area_brain_size <- plot_cv_brain_size(data2analysis, data2analysis$area_cv, "Surface Area")
plot_volume_brain_size <- plot_cv_brain_size(data2analysis, data2analysis$volume_cv, "Cortical Volume")


read_clean_scale <- function(file) {
  data <- read.csv(file)
  data <- data[grepl("Thalamus-Proper|Putamen|Hippocampus|Caudate|Pallidum|Amygdala|Accumbens-area", 
                                                     data$X), ]
  scaled <- scale(data[,2:ncol(data)])
  return(scaled)
}

sub_files <- c("sub_cv_volume_sub_qc.csv", "sub_mean_volume_sub_qc.csv")

names <- c("sub_cv", "sub_mean")

data_list <- lapply(sub_files, FUN = read_clean_scale)
names(data_list) <- names

sub_long_data <- lapply(seq_along(data_list), function(i) {
  pivot_long(data_list[[i]], names[i])
}) 
sub_data2analysis <- Reduce(function(dtf1, dtf2) {cbind(dtf1, dtf2[,-1])}, sub_long_data)

cor_test <- cor.test(sub_data2analysis$sub_mean, sub_data2analysis$sub_cv)
print(cor_test)


plot_subvolume_brain_size <-   ggplot(sub_data2analysis, aes(x = sub_mean, y = sub_cv, color = subj)) +
  geom_point(size = 2) +
  geom_smooth(method=lm, se=TRUE, color="black", fill="grey") +
  labs(title = "Subcortical Volume", x = "Scaled Brain Size", y = "Scaled CV") +
  scale_color_manual(values=color_mapping) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 18, family = "Arial", face = "bold", colour = "black", hjust = 0.5),
    axis.title.x = element_text(size = 18, family = "Arial", face = "bold", colour = "black", vjust = -0.5),
    axis.title.y = element_text(size = 18, family = "Arial", face = "bold", colour = "black", vjust = 1.5),
    axis.text.x = element_text(size = 18, family = "Arial", face = "bold", colour = "black"),
    axis.text.y = element_text(size = 18, family = "Arial", face = "bold", colour = "black"),
    axis.ticks.length = unit(0.2, "cm"),
    axis.line = element_line(colour = 'black', linewidth = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank(),
    axis.text.x.top = element_text(angle = 90),
    axis.text.y.right = element_text(angle = 90),
    axis.line.x.top = element_line(color = "black"),
    axis.line.y.right = element_line(color = "black"),
    legend.position = "none"
  )


Fig9 <- plot_thickness_brain_size + theme(axis.title.x = element_blank(), plot.margin = margin(1.2, 1.2, 1.2, 1.2, "cm"))+
  plot_area_brain_size + theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  plot_volume_brain_size +
  plot_subvolume_brain_size + theme(axis.title.y = element_blank())
  
ggsave("/Users/wang/Library/CloudStorage/OneDrive-UniversityofBergen/Desktop/UiB/Manuscripts/BBSC/4.T1w_Image/Figures/Fig9.tiff", 
       plot = Fig9, width = 10, height = 10, units = 'in', dpi = 300, compression = "lzw")

