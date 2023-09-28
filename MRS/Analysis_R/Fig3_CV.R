# This code is to plot the Fig. 3

# WANG. 28-April-2023
# Updated 28-Aug-2023
# Set the working directory to the path where your files are located
setwd("/Users/wang/Desktop/Research_projects/BBSC/MRS/Output_v2_4/QuantifyResults")


# Remove all objects created before to prevent clusing
rm()

# Load required packages
library(ggplot2)
library(gridExtra)
library(grid)
library(berryFunctions)

# Define file names and column names
file_names <- "A_TissCorrWaterScaled_Voxel_1_Basis_1.tsv"

mrs_of_interest <- c("Cr", "PCr", "Cr_PCr", "Glu", "Glu_Gln", "NAA", "NAA_NAAG", "PCh_GPC", "mI")
mrs_of_name <- c("Cr", "PCr", "tCr", "Glu", "Glx", "NAA", "tNAA", "tCho", "mI")

# Read in the data
MRS_data <- read.table(file_names, header = TRUE, sep = "\t")

# extract the data for each sub
MRS_data_sub1 <- MRS_data[1:38, mrs_of_interest]
MRS_data_sub2 <- MRS_data[39:77, mrs_of_interest]
MRS_data_sub3 <- MRS_data[78:102, mrs_of_interest]


# insert nan values to make it even so, 1:23 sessions are from Jan.to May.
MRS_data_sub1 <- insertRows(MRS_data_sub1, 23 , new = NA)
MRS_data_sub3 <- insertRows(MRS_data_sub3, 20:23 , new = NA)

# insert nan values to make the gap, 24:27 sessions are from Jun.to Oct.
MRS_data_sub1 <- insertRows(MRS_data_sub1, 24:27 , new = NA)
MRS_data_sub2 <- insertRows(MRS_data_sub2, 24:27 , new = NA)
MRS_data_sub3 <- insertRows(MRS_data_sub3, 24:27 , new = NA)

# insert nan values to make it even, 28:38 sessions are from Nov.to Dec.
MRS_data_sub3[34:38,] <- NA

# insert nan values to make it even, 39:40 sessions are from jan. 2022
MRS_data_sub1 <- insertRows(MRS_data_sub1, 39:40 , new = NA)
MRS_data_sub2 <- insertRows(MRS_data_sub2, 39:40 , new = NA)
MRS_data_sub3[39:40,] <- NA

# insert nan values to make it even, 41:45 sessions are from Feb. 2022
MRS_data_sub3[41:45,] <- NA

# # Add NA values to sub1
# n_missing_rows <- nrow(MRS_data_sub2) - nrow(MRS_data_sub1)
# temp <- data.frame(matrix(NA, ncol = 9, nrow = n_missing_rows))
# names(temp) <- mrs_of_interest
# MRS_data_sub1 <- rbind(MRS_data_sub1, temp)
# 
# # Add NA values to sub3
# n_missing_rows <- nrow(MRS_data_sub2) - nrow(MRS_data_sub3)
# temp <- data.frame(matrix(NA, ncol = 9, nrow = n_missing_rows))
# names(temp) <- mrs_of_interest
# MRS_data_sub3 <- rbind(MRS_data_sub3, temp)

# calculate mean and limits
MRS_mean_sub1 <- colMeans(MRS_data_sub1,na.rm=TRUE)
MRS_mean_sub2 <- colMeans(MRS_data_sub2,na.rm=TRUE)
MRS_mean_sub3 <- colMeans(MRS_data_sub3,na.rm=TRUE)


# Set base font size
base_size <- 12

p1_cr <- ggplot() +
  geom_point(data=MRS_data_sub1, aes(x = seq(1, 45), y = (MRS_data_sub1[,mrs_of_interest[1]]/MRS_mean_sub1[mrs_of_interest[1]] - 1)*100, color = "#D40000")) +
  geom_point(data=MRS_data_sub2, aes(x = seq(1, 45), y = (MRS_data_sub2[,mrs_of_interest[1]]/MRS_mean_sub2[mrs_of_interest[1]] - 1)*100, color = "#377EB8")) +
  geom_point(data=MRS_data_sub3, aes(x = seq(1, 45), y = (MRS_data_sub3[,mrs_of_interest[1]]/MRS_mean_sub3[mrs_of_interest[1]] - 1)*100, color = "#4DAF4A")) +
  geom_hline (yintercept = 15, linetype = "dashed", color = "grey",size = 1)+
  geom_hline (yintercept = -15, linetype = "dashed", color = "grey", size =1)+
  geom_segment(aes(x = 0.5, xend = 0.5, y = -15, yend = 15), linetype = "solid", color = "#6A3D9A", size =1) +
  geom_segment(aes(x = 23.5, xend = 23.5, y = -15, yend = 15), linetype = "dashed", color = "#6A3D9A", size =1) +
  geom_segment(aes(x = 27.5, xend = 27.5, y = -15, yend = 15), linetype = "dashed", color = "#6A3D9A", size =1) +
  geom_segment(aes(x = 38.5, xend = 38.5, y = -15, yend = 15), linetype = "solid", color = "#6A3D9A", size =1) +
  geom_segment(aes(x = 40.5, xend = 40.5, y = -15, yend = 15), linetype = "solid", color = "#FFA500", size =1) +
  geom_segment(aes(x = 45.5, xend = 45.5, y = -15, yend = 15), linetype = "solid", color = "#FFA500", size =1) +
  xlab("")+
  ylab(mrs_of_name[1]) +
  ylim(-60,60)+
  theme_bw(base_size = base_size) +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(size = 1, color = "black"),
    axis.text = element_text(size = 18, face = "bold", family = "sans", color = "black"),
    axis.title = element_text(size = 20, face = "bold"),
    axis.line.x = element_blank(), # remove x-axis line
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(), # remove x-axis tick marks
    legend.position = "none" # remove color legend
  )
p1_cr


p2_pcr <- ggplot() +
  geom_point(data=MRS_data_sub1, aes(x = seq(1, 45), y = (MRS_data_sub1[,mrs_of_interest[2]]/MRS_mean_sub1[mrs_of_interest[2]] - 1)*100, color = "#D40000")) +
  geom_point(data=MRS_data_sub2, aes(x = seq(1, 45), y = (MRS_data_sub2[,mrs_of_interest[2]]/MRS_mean_sub2[mrs_of_interest[2]] - 1)*100, color = "#377EB8")) +
  geom_point(data=MRS_data_sub3, aes(x = seq(1, 45), y = (MRS_data_sub3[,mrs_of_interest[2]]/MRS_mean_sub3[mrs_of_interest[2]] - 1)*100, color = "#4DAF4A")) +
  geom_hline (yintercept = 15, linetype = "dashed", color = "grey",size = 1)+
  geom_hline (yintercept = -15, linetype = "dashed", color = "grey", size =1)+
  geom_segment(aes(x = 0.5, xend = 0.5, y = -15, yend = 15), linetype = "solid", color = "#6A3D9A", size =1) +
  geom_segment(aes(x = 23.5, xend = 23.5, y = -15, yend = 15), linetype = "dashed", color = "#6A3D9A", size =1) +
  geom_segment(aes(x = 27.5, xend = 27.5, y = -15, yend = 15), linetype = "dashed", color = "#6A3D9A", size =1) +
  geom_segment(aes(x = 38.5, xend = 38.5, y = -15, yend = 15), linetype = "solid", color = "#6A3D9A", size =1) +
  geom_segment(aes(x = 40.5, xend = 40.5, y = -15, yend = 15), linetype = "solid", color = "#FFA500", size =1) +
  geom_segment(aes(x = 45.5, xend = 45.5, y = -15, yend = 15), linetype = "solid", color = "#FFA500", size =1) +
  xlab("")+
  ylab(mrs_of_name[2]) +
  ylim(-60,60)+
  theme_bw(base_size = base_size) +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(size = 1, color = "black"),
    axis.text = element_text(size = 18, face = "bold", family = "sans", color = "black"),
    axis.title = element_text(size = 20, face = "bold"),
    axis.line.x = element_blank(), # remove x-axis line
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(), # remove x-axis tick marks
    legend.position = "none" # remove color legend
  )
p2_pcr

p3_tcr <- ggplot() +
  geom_point(data=MRS_data_sub1, aes(x = seq(1, 45), y = (MRS_data_sub1[,mrs_of_interest[3]]/MRS_mean_sub1[mrs_of_interest[3]] - 1)*100, color = "#D40000")) +
  geom_point(data=MRS_data_sub2, aes(x = seq(1, 45), y = (MRS_data_sub2[,mrs_of_interest[3]]/MRS_mean_sub2[mrs_of_interest[3]] - 1)*100, color = "#377EB8")) +
  geom_point(data=MRS_data_sub3, aes(x = seq(1, 45), y = (MRS_data_sub3[,mrs_of_interest[3]]/MRS_mean_sub3[mrs_of_interest[3]] - 1)*100, color = "#4DAF4A")) +
  geom_hline (yintercept = 15, linetype = "dashed", color = "grey",size = 1)+
  geom_hline (yintercept = -15, linetype = "dashed", color = "grey", size =1)+
  geom_segment(aes(x = 0.5, xend = 0.5, y = -15, yend = 15), linetype = "solid", color = "#6A3D9A", size =1) +
  geom_segment(aes(x = 23.5, xend = 23.5, y = -15, yend = 15), linetype = "dashed", color = "#6A3D9A", size =1) +
  geom_segment(aes(x = 27.5, xend = 27.5, y = -15, yend = 15), linetype = "dashed", color = "#6A3D9A", size =1) +
  geom_segment(aes(x = 38.5, xend = 38.5, y = -15, yend = 15), linetype = "solid", color = "#6A3D9A", size =1) +
  geom_segment(aes(x = 40.5, xend = 40.5, y = -15, yend = 15), linetype = "solid", color = "#FFA500", size =1) +
  geom_segment(aes(x = 45.5, xend = 45.5, y = -15, yend = 15), linetype = "solid", color = "#FFA500", size =1) +
  xlab("")+
  ylab(mrs_of_name[3]) +
  ylim(-60,60)+
  theme_bw(base_size = base_size) +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(size = 1, color = "black"),
    axis.text = element_text(size = 18, face = "bold", family = "sans", color = "black"),
    axis.title = element_text(size = 20, face = "bold"),
    axis.line.x = element_blank(), # remove x-axis line
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(), # remove x-axis tick marks
    legend.position = "none" # remove color legend
  )
p3_tcr

p4_glu <- ggplot() +
  geom_point(data=MRS_data_sub1, aes(x = seq(1, 45), y = (MRS_data_sub1[,mrs_of_interest[4]]/MRS_mean_sub1[mrs_of_interest[4]] - 1)*100, color = "#D40000")) +
  geom_point(data=MRS_data_sub2, aes(x = seq(1, 45), y = (MRS_data_sub2[,mrs_of_interest[4]]/MRS_mean_sub2[mrs_of_interest[4]] - 1)*100, color = "#377EB8")) +
  geom_point(data=MRS_data_sub3, aes(x = seq(1, 45), y = (MRS_data_sub3[,mrs_of_interest[4]]/MRS_mean_sub3[mrs_of_interest[4]] - 1)*100, color = "#4DAF4A")) +
  geom_hline (yintercept = 15, linetype = "dashed", color = "grey",size = 1)+
  geom_hline (yintercept = -15, linetype = "dashed", color = "grey", size =1)+
  geom_segment(aes(x = 0.5, xend = 0.5, y = -15, yend = 15), linetype = "solid", color = "#6A3D9A", size =1) +
  geom_segment(aes(x = 23.5, xend = 23.5, y = -15, yend = 15), linetype = "dashed", color = "#6A3D9A", size =1) +
  geom_segment(aes(x = 27.5, xend = 27.5, y = -15, yend = 15), linetype = "dashed", color = "#6A3D9A", size =1) +
  geom_segment(aes(x = 38.5, xend = 38.5, y = -15, yend = 15), linetype = "solid", color = "#6A3D9A", size =1) +
  geom_segment(aes(x = 40.5, xend = 40.5, y = -15, yend = 15), linetype = "solid", color = "#FFA500", size =1) +
  geom_segment(aes(x = 45.5, xend = 45.5, y = -15, yend = 15), linetype = "solid", color = "#FFA500", size =1) +
  # geom_text(aes(x = 12, y = -55, label = "Jan.>May"), face = "bold", color = "#6A3D9A", size = 5) +
  # geom_text(aes(x = 33, y = -55, label = "Nov.&Dec."), face = "bold", color = "#6A3D9A", size = 5) +
  # geom_text(aes(x = 43.5, y = -55, label = "Feb."), face = "bold", color = "#FFA500", size = 5) +
  xlab("")+
  ylab(mrs_of_name[4]) +
  ylim(-60,60)+
  theme_bw(base_size = base_size) +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(size = 1, color = "black"),
    axis.text = element_text(size = 18, face = "bold", family = "sans", color = "black"),
    axis.title = element_text(size = 20, face = "bold"),
    axis.line.x = element_blank(), # remove x-axis line
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(), # remove x-axis tick marks
    legend.position = "none" # remove color legend
  )
p4_glu


p5_glx <- ggplot() +
  geom_point(data=MRS_data_sub1, aes(x = seq(1, 45), y = (MRS_data_sub1[,mrs_of_interest[5]]/MRS_mean_sub1[mrs_of_interest[5]] - 1)*100, color = "#D40000")) +
  geom_point(data=MRS_data_sub2, aes(x = seq(1, 45), y = (MRS_data_sub2[,mrs_of_interest[5]]/MRS_mean_sub2[mrs_of_interest[5]] - 1)*100, color = "#377EB8")) +
  geom_point(data=MRS_data_sub3, aes(x = seq(1, 45), y = (MRS_data_sub3[,mrs_of_interest[5]]/MRS_mean_sub3[mrs_of_interest[5]] - 1)*100, color = "#4DAF4A")) +
  geom_hline (yintercept = 15, linetype = "dashed", color = "grey",size = 1)+
  geom_hline (yintercept = -15, linetype = "dashed", color = "grey", size =1)+
  geom_segment(aes(x = 0.5, xend = 0.5, y = -15, yend = 15), linetype = "solid", color = "#6A3D9A", size =1) +
  geom_segment(aes(x = 23.5, xend = 23.5, y = -15, yend = 15), linetype = "dashed", color = "#6A3D9A", size =1) +
  geom_segment(aes(x = 27.5, xend = 27.5, y = -15, yend = 15), linetype = "dashed", color = "#6A3D9A", size =1) +
  geom_segment(aes(x = 38.5, xend = 38.5, y = -15, yend = 15), linetype = "solid", color = "#6A3D9A", size =1) +
  geom_segment(aes(x = 40.5, xend = 40.5, y = -15, yend = 15), linetype = "solid", color = "#FFA500", size =1) +
  geom_segment(aes(x = 45.5, xend = 45.5, y = -15, yend = 15), linetype = "solid", color = "#FFA500", size =1) +
  # geom_text(aes(x = 12, y = -55, label = "Jan.>May"), face = "bold", color = "#6A3D9A", size = 5) +
  # geom_text(aes(x = 33, y = -55, label = "Nov.&Dec."), face = "bold", color = "#6A3D9A", size = 5) +
  # geom_text(aes(x = 43.5, y = -55, label = "Feb."), face = "bold", color = "#FFA500", size = 5) +
  xlab("")+
  ylab(mrs_of_name[5]) +
  ylim(-60,60)+
  theme_bw(base_size = base_size) +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(size = 1, color = "black"),
    axis.text = element_text(size = 18, face = "bold", family = "sans", color = "black"),
    axis.title = element_text(size = 20, face = "bold"),
    axis.line.x = element_blank(), # remove x-axis line
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(), # remove x-axis tick marks
    legend.position = "none" # remove color legend
  )
p5_glx

p6_naa <- ggplot() +
  geom_point(data=MRS_data_sub1, aes(x = seq(1, 45), y = (MRS_data_sub1[,mrs_of_interest[6]]/MRS_mean_sub1[mrs_of_interest[6]] - 1)*100, color = "#D40000")) +
  geom_point(data=MRS_data_sub2, aes(x = seq(1, 45), y = (MRS_data_sub2[,mrs_of_interest[6]]/MRS_mean_sub2[mrs_of_interest[6]] - 1)*100, color = "#377EB8")) +
  geom_point(data=MRS_data_sub3, aes(x = seq(1, 45), y = (MRS_data_sub3[,mrs_of_interest[6]]/MRS_mean_sub3[mrs_of_interest[6]] - 1)*100, color = "#4DAF4A")) +
  geom_hline (yintercept = 15, linetype = "dashed", color = "grey",size = 1)+
  geom_hline (yintercept = -15, linetype = "dashed", color = "grey", size =1)+
  geom_segment(aes(x = 0.5, xend = 0.5, y = -15, yend = 15), linetype = "solid", color = "#6A3D9A", size =1) +
  geom_segment(aes(x = 23.5, xend = 23.5, y = -15, yend = 15), linetype = "dashed", color = "#6A3D9A", size =1) +
  geom_segment(aes(x = 27.5, xend = 27.5, y = -15, yend = 15), linetype = "dashed", color = "#6A3D9A", size =1) +
  geom_segment(aes(x = 38.5, xend = 38.5, y = -15, yend = 15), linetype = "solid", color = "#6A3D9A", size =1) +
  geom_segment(aes(x = 40.5, xend = 40.5, y = -15, yend = 15), linetype = "solid", color = "#FFA500", size =1) +
  geom_segment(aes(x = 45.5, xend = 45.5, y = -15, yend = 15), linetype = "solid", color = "#FFA500", size =1) +
  xlab("")+
  ylab(mrs_of_name[6]) +
  ylim(-60,60)+
  theme_bw(base_size = base_size) +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(size = 1, color = "black"),
    axis.text = element_text(size = 18, face = "bold", family = "sans", color = "black"),
    axis.title = element_text(size = 20, face = "bold"),
    axis.line.x = element_blank(), # remove x-axis line
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(), # remove x-axis tick marks
    legend.position = "none" # remove color legend
  )
p6_naa

p7_tnaa <- ggplot() +
  geom_point(data=MRS_data_sub1, aes(x = seq(1, 45), y = (MRS_data_sub1[,mrs_of_interest[7]]/MRS_mean_sub1[mrs_of_interest[7]] - 1)*100, color = "#D40000")) +
  geom_point(data=MRS_data_sub2, aes(x = seq(1, 45), y = (MRS_data_sub2[,mrs_of_interest[7]]/MRS_mean_sub2[mrs_of_interest[7]] - 1)*100, color = "#377EB8")) +
  geom_point(data=MRS_data_sub3, aes(x = seq(1, 45), y = (MRS_data_sub3[,mrs_of_interest[7]]/MRS_mean_sub3[mrs_of_interest[7]] - 1)*100, color = "#4DAF4A")) +
  geom_hline (yintercept = 15, linetype = "dashed", color = "grey",size = 1)+
  geom_hline (yintercept = -15, linetype = "dashed", color = "grey", size =1)+
  geom_segment(aes(x = 0.5, xend = 0.5, y = -15, yend = 15), linetype = "solid", color = "#6A3D9A", size =1) +
  geom_segment(aes(x = 23.5, xend = 23.5, y = -15, yend = 15), linetype = "dashed", color = "#6A3D9A", size =1) +
  geom_segment(aes(x = 27.5, xend = 27.5, y = -15, yend = 15), linetype = "dashed", color = "#6A3D9A", size =1) +
  geom_segment(aes(x = 38.5, xend = 38.5, y = -15, yend = 15), linetype = "solid", color = "#6A3D9A", size =1) +
  geom_segment(aes(x = 40.5, xend = 40.5, y = -15, yend = 15), linetype = "solid", color = "#FFA500", size =1) +
  geom_segment(aes(x = 45.5, xend = 45.5, y = -15, yend = 15), linetype = "solid", color = "#FFA500", size =1) +
  xlab("")+
  ylab(mrs_of_name[7]) +
  ylim(-60,60)+
  theme_bw(base_size = base_size) +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(size = 1, color = "black"),
    axis.text = element_text(size = 18, face = "bold", family = "sans", color = "black"),
    axis.title = element_text(size = 20, face = "bold"),
    axis.line.x = element_blank(), # remove x-axis line
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(), # remove x-axis tick marks
    legend.position = "none" # remove color legend
  )
p7_tnaa


p8_tch <- ggplot() +
  geom_point(data=MRS_data_sub1, aes(x = seq(1, 45), y = (MRS_data_sub1[,mrs_of_interest[8]]/MRS_mean_sub1[mrs_of_interest[8]] - 1)*100, color = "#D40000")) +
  geom_point(data=MRS_data_sub2, aes(x = seq(1, 45), y = (MRS_data_sub2[,mrs_of_interest[8]]/MRS_mean_sub2[mrs_of_interest[8]] - 1)*100, color = "#377EB8")) +
  geom_point(data=MRS_data_sub3, aes(x = seq(1, 45), y = (MRS_data_sub3[,mrs_of_interest[8]]/MRS_mean_sub3[mrs_of_interest[8]] - 1)*100, color = "#4DAF4A")) +
  geom_hline (yintercept = 15, linetype = "dashed", color = "grey",size = 1)+
  geom_hline (yintercept = -15, linetype = "dashed", color = "grey", size =1)+
  geom_segment(aes(x = 0.5, xend = 0.5, y = -15, yend = 15), linetype = "solid", color = "#6A3D9A", size =1) +
  geom_segment(aes(x = 23.5, xend = 23.5, y = -15, yend = 15), linetype = "dashed", color = "#6A3D9A", size =1) +
  geom_segment(aes(x = 27.5, xend = 27.5, y = -15, yend = 15), linetype = "dashed", color = "#6A3D9A", size =1) +
  geom_segment(aes(x = 38.5, xend = 38.5, y = -15, yend = 15), linetype = "solid", color = "#6A3D9A", size =1) +
  geom_segment(aes(x = 40.5, xend = 40.5, y = -15, yend = 15), linetype = "solid", color = "#FFA500", size =1) +
  geom_segment(aes(x = 45.5, xend = 45.5, y = -15, yend = 15), linetype = "solid", color = "#FFA500", size =1) +
  xlab("")+
  ylab(mrs_of_name[8]) +
  ylim(-60,60)+
  theme_bw(base_size = base_size) +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(size = 1, color = "black"),
    axis.text = element_text(size = 18, face = "bold", family = "sans", color = "black"),
    axis.title = element_text(size = 20, face = "bold"),
    axis.line.x = element_blank(), # remove x-axis line
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(), # remove x-axis tick marks
    legend.position = "none" # remove color legend
  )
p8_tch

p9_mi <- ggplot() +
  geom_point(data=MRS_data_sub1, aes(x = seq(1, 45), y = (MRS_data_sub1[,mrs_of_interest[9]]/MRS_mean_sub1[mrs_of_interest[9]] - 1)*100, color = "#D40000")) +
  geom_point(data=MRS_data_sub2, aes(x = seq(1, 45), y = (MRS_data_sub2[,mrs_of_interest[9]]/MRS_mean_sub2[mrs_of_interest[9]] - 1)*100, color = "#377EB8")) +
  geom_point(data=MRS_data_sub3, aes(x = seq(1, 45), y = (MRS_data_sub3[,mrs_of_interest[9]]/MRS_mean_sub3[mrs_of_interest[9]] - 1)*100, color = "#4DAF4A")) + 
  geom_hline (yintercept = 15, linetype = "dashed", color = "grey",size = 1)+
  geom_hline (yintercept = -15, linetype = "dashed", color = "grey", size =1)+
  geom_segment(aes(x = 0.5, xend = 0.5, y = -15, yend = 15), linetype = "solid", color = "#6A3D9A", size =1) +
  geom_segment(aes(x = 23.5, xend = 23.5, y = -15, yend = 15), linetype = "dashed", color = "#6A3D9A", size =1) +
  geom_segment(aes(x = 27.5, xend = 27.5, y = -15, yend = 15), linetype = "dashed", color = "#6A3D9A", size =1) +
  geom_segment(aes(x = 38.5, xend = 38.5, y = -15, yend = 15), linetype = "solid", color = "#6A3D9A", size =1) +
  geom_segment(aes(x = 40.5, xend = 40.5, y = -15, yend = 15), linetype = "solid", color = "#FFA500", size =1) +
  geom_segment(aes(x = 45.5, xend = 45.5, y = -15, yend = 15), linetype = "solid", color = "#FFA500", size =1) +
  # geom_text(aes(x = 12, y = -55, label = "Jan.>May"), face = "bold", color = "#6A3D9A", size = 5) +
  # geom_text(aes(x = 33, y = -55, label = "Nov.&Dec."), face = "bold", color = "#6A3D9A", size = 5) +
  # geom_text(aes(x = 43.5, y = -55, label = "Feb."), face = "bold", color = "#FFA500", size = 5) +
  xlab("")+
  ylab(mrs_of_name[9]) +
  ylim(-60,60)+
  theme_bw(base_size = base_size) +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(size = 1, color = "black"),
    axis.text = element_text(size = 18, face = "bold", family = "sans", color = "black"),
    axis.title = element_text(size = 20, face = "bold"),
    axis.line.x = element_blank(), # remove x-axis line
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(), # remove x-axis tick marks
    legend.position = "none" # remove color legend
  )
p9_mi


## change the limits for each subplots
p1_cr <- p1_cr + scale_y_continuous(limits = c(-60, 60), breaks = c(-50,0, 50)) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_line(),
        axis.title.y = element_text(vjust = 0))
p2_pcr <- p2_pcr + scale_y_continuous(limits = c(-60, 60), breaks = c(-50,0, 50))
p3_tcr <- p3_tcr + scale_y_continuous(limits = c(-25, 25), breaks = c(-15,0, 15))
p4_glu <- p4_glu + scale_y_continuous(limits = c(-25, 25), breaks = c(-15,0, 15))+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_line(),
        axis.title.y = element_text(vjust = -6))
p5_glx <- p5_glx + scale_y_continuous(limits = c(-25, 25), breaks = c(-15,0, 15))
p6_naa <- p6_naa + scale_y_continuous(limits = c(-25, 25), breaks = c(-15,0, 15))+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_line(),
        axis.title.y = element_text(vjust = -6))
p7_tnaa <- p7_tnaa + scale_y_continuous(limits = c(-25, 25), breaks = c(-15,0, 15))
p8_tch <- p8_tch + scale_y_continuous(limits = c(-25, 25), breaks = c(-15,0, 15))+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_line(),
        axis.title.y = element_text(vjust = 0))
p9_mi <- p9_mi + scale_y_continuous(limits = c(-25, 25), breaks = c(-15,0, 15))+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_line(),
        axis.title.y = element_text(vjust = 0))

library(patchwork)
(p3_tcr + p2_pcr + p1_cr) / (p8_tch + p9_mi + plot_spacer())/ (p7_tnaa + p6_naa + plot_spacer())/ (p5_glx + p4_glu + plot_spacer())   

p3_tcr + p2_pcr + p1_cr +  p7_tnaa + p6_naa + p8_tch +  p5_glx + p4_glu  + p9_mi

