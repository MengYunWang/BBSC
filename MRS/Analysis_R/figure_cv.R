# This code is to plot the distribution of metabolites in Sub1

# WANG. 28-April-2023

# Set the working directory to the path where your files are located
setwd("/Users/wang/Desktop/Research_projects/BBSC/MRS/Output_v2_4/QuantifyResults")


# Remove all objects created before to prevent clusing
rm()

# Load required packages
library(ggplot2)
library(gridExtra)
library(grid)

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

# Add NA values to sub1
n_missing_rows <- nrow(MRS_data_sub2) - nrow(MRS_data_sub1)
temp <- data.frame(matrix(NA, ncol = 9, nrow = n_missing_rows))
names(temp) <- mrs_of_interest
MRS_data_sub1 <- rbind(MRS_data_sub1, temp)

# Add NA values to sub3
n_missing_rows <- nrow(MRS_data_sub2) - nrow(MRS_data_sub3)
temp <- data.frame(matrix(NA, ncol = 9, nrow = n_missing_rows))
names(temp) <- mrs_of_interest
MRS_data_sub3 <- rbind(MRS_data_sub3, temp)

# calculate mean and limits
MRS_mean_sub1 <- colMeans(MRS_data_sub1)
MRS_mean_sub2 <- colMeans(MRS_data_sub2)
MRS_mean_sub3 <- colMeans(MRS_data_sub3)


# Set base font size
base_size <- 12

p1_cr <- ggplot() +
  geom_point(data=MRS_data_sub1, aes(x = seq(1, 39), y = MRS_data_sub1[,mrs_of_interest[1]], color = "#E41A1C")) +
  geom_point(data=MRS_data_sub2, aes(x = seq(1, 39), y = MRS_data_sub2[,mrs_of_interest[1]], color = "#377EB8")) +
  geom_point(data=MRS_data_sub3, aes(x = seq(1, 39), y = MRS_data_sub3[,mrs_of_interest[1]], color = "#4DAF4A")) +
  xlab("")+
  ylab(mrs_of_name[1]) +
  ylim(0,15)+
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
  geom_point(data=MRS_data_sub1, aes(x = seq(1, 39), y = MRS_data_sub1[,mrs_of_interest[2]], color = "#E41A1C")) +
  geom_point(data=MRS_data_sub2, aes(x = seq(1, 39), y = MRS_data_sub2[,mrs_of_interest[2]], color = "#377EB8")) +
  geom_point(data=MRS_data_sub3, aes(x = seq(1, 39), y = MRS_data_sub3[,mrs_of_interest[2]], color = "#4DAF4A")) +
  xlab("")+
  ylab(mrs_of_name[2]) +
  ylim(0,15)+
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
  geom_point(data=MRS_data_sub1, aes(x = seq(1, 39), y = MRS_data_sub1[,mrs_of_interest[3]], color = "#E41A1C")) +
  geom_point(data=MRS_data_sub2, aes(x = seq(1, 39), y = MRS_data_sub2[,mrs_of_interest[3]], color = "#377EB8")) +
  geom_point(data=MRS_data_sub3, aes(x = seq(1, 39), y = MRS_data_sub3[,mrs_of_interest[3]], color = "#4DAF4A")) +
  geom_hline(yintercept = MRS_mean_sub1[3], linetype = "solid") +

  xlab("")+
  ylab(mrs_of_name[3]) +
  ylim(0,15)+
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
  geom_point(data=MRS_data_sub1, aes(x = seq(1, 39), y = MRS_data_sub1[,mrs_of_interest[4]], color = "#E41A1C")) +
  geom_point(data=MRS_data_sub2, aes(x = seq(1, 39), y = MRS_data_sub2[,mrs_of_interest[4]], color = "#377EB8")) +
  geom_point(data=MRS_data_sub3, aes(x = seq(1, 39), y = MRS_data_sub3[,mrs_of_interest[4]], color = "#4DAF4A")) +
  geom_hline(yintercept = MRS_mean_sub1[4], linetype = "solid") +

  xlab("")+
  ylab(mrs_of_name[4]) +
  ylim(15,30)+
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
  geom_point(data=MRS_data_sub1, aes(x = seq(1, 39), y = MRS_data_sub1[,mrs_of_interest[5]], color = "#E41A1C")) +
  geom_point(data=MRS_data_sub2, aes(x = seq(1, 39), y = MRS_data_sub2[,mrs_of_interest[5]], color = "#377EB8")) +
  geom_point(data=MRS_data_sub3, aes(x = seq(1, 39), y = MRS_data_sub3[,mrs_of_interest[5]], color = "#4DAF4A")) +
  geom_hline(yintercept = MRS_mean_sub1[5], linetype = "solid") +

  xlab("")+
  ylab(mrs_of_name[5]) +
  ylim(15,30)+
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
  geom_point(data=MRS_data_sub1, aes(x = seq(1, 39), y = MRS_data_sub1[,mrs_of_interest[6]], color = "#E41A1C")) +
  geom_point(data=MRS_data_sub2, aes(x = seq(1, 39), y = MRS_data_sub2[,mrs_of_interest[6]], color = "#377EB8")) +
  geom_point(data=MRS_data_sub3, aes(x = seq(1, 39), y = MRS_data_sub3[,mrs_of_interest[6]], color = "#4DAF4A")) +
  geom_hline(yintercept = MRS_mean_sub1[6], linetype = "solid") +

  xlab("")+
  ylab(mrs_of_name[6]) +
  ylim(0,30)+
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
  geom_point(data=MRS_data_sub1, aes(x = seq(1, 39), y = MRS_data_sub1[,mrs_of_interest[7]], color = "#E41A1C")) +
  geom_point(data=MRS_data_sub2, aes(x = seq(1, 39), y = MRS_data_sub2[,mrs_of_interest[7]], color = "#377EB8")) +
  geom_point(data=MRS_data_sub3, aes(x = seq(1, 39), y = MRS_data_sub3[,mrs_of_interest[7]], color = "#4DAF4A")) +
  geom_hline(yintercept = MRS_mean_sub1[7], linetype = "solid") +

  xlab("")+
  ylab(mrs_of_name[7]) +
  ylim(0,30)+
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
  geom_point(data=MRS_data_sub1, aes(x = seq(1, 39), y = MRS_data_sub1[,mrs_of_interest[8]], color = "#E41A1C")) +
  geom_point(data=MRS_data_sub2, aes(x = seq(1, 39), y = MRS_data_sub2[,mrs_of_interest[8]], color = "#377EB8")) +
  geom_point(data=MRS_data_sub3, aes(x = seq(1, 39), y = MRS_data_sub3[,mrs_of_interest[8]], color = "#4DAF4A")) +
  geom_hline(yintercept = MRS_mean_sub1[8], linetype = "solid") +

  xlab("")+
  ylab(mrs_of_name[8]) +
  ylim(0,30)+
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
  geom_point(data=MRS_data_sub1, aes(x = seq(1, 39), y = MRS_data_sub1[,mrs_of_interest[9]], color = "#E41A1C")) +
  geom_point(data=MRS_data_sub2, aes(x = seq(1, 39), y = MRS_data_sub2[,mrs_of_interest[9]], color = "#377EB8")) +
  geom_point(data=MRS_data_sub3, aes(x = seq(1, 39), y = MRS_data_sub3[,mrs_of_interest[9]], color = "#4DAF4A")) + 
  geom_hline(yintercept = MRS_mean_sub1[9], linetype = "solid") +
  xlab("")+
  ylab(mrs_of_name[9]) +
  ylim(0,30)+
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


# change the limits for each subplots
p1_cr <- p1_cr + scale_y_continuous(limits = c(0, 15), breaks = c(0, 15))
p2_pcr <- p2_pcr + scale_y_continuous(limits = c(0, 15), breaks = c(0, 15))
p3_tcr <- p3_tcr + scale_y_continuous(limits = c(0, 15), breaks = c(0, 15))
p4_glu <- p4_glu + scale_y_continuous(limits = c(15, 30), breaks = c(15, 30))
p5_glx <- p5_glx + scale_y_continuous(limits = c(15, 30), breaks = c(15, 30))
p6_naa <- p6_naa + scale_y_continuous(limits = c(10, 25), breaks = c(10, 25))
p7_tnaa <- p7_tnaa + scale_y_continuous(limits = c(10, 25), breaks = c(10, 25))
p8_tch <- p8_tch + scale_y_continuous(limits = c(0, 15), breaks = c(0, 15))
p9_mi <- p9_mi + scale_y_continuous(limits = c(0, 15), breaks = c(0, 15))

library(patchwork)
(p3_tcr + p2_pcr + p1_cr) / (p8_tch + p9_mi + plot_spacer())/ (p7_tnaa + p6_naa + plot_spacer())/ (p5_glx + p4_glu + plot_spacer())   


