# This code is to plot the Fig.1 
# a. distribution of the IQM
# b. correlations between IQMs
# c. the percentage change of the mean value of the different pipelines

# WANG. 27-Mar-2024


# Remove all objects created before to prevent clash
rm(list = ls())

library(dplyr)
library(ggplot2)
library(patchwork)
library(tidyr)
library(readr)
library(reshape2)
library(RColorBrewer)  # For color palettes

# Set the working directory to the path where your files are located
setwd("/Users/wang/Desktop/Research_projects/BBSC/Functional/fmri/derivatives/freesurfer-v7.2.0/")
# setwd("//Users/joeywang/Desktop/BBSC/Functional/fmri/derivatives/freesurfer-v7.2.0/")


# --load the qc data and replace 1 2 3 with 01 02 03,
# so that it can corresponds to the data collection sequence
qc_data <- read_tsv(file = "group_T1w_v23.tsv") %>%
  mutate(bids_name = sub("ses-([1-9])([^0-9])", "ses-0\\1\\2", bids_name)) %>%
  arrange(bids_name)

# IQMs that we are interested in
qc_indx <-
  c("cjv",
    "efc",
    "fwhm_avg",
    "inu_med",
    "fber",
    "wm2max",
    "rpve_gm",
    "snr_total",
    "cnr")
qc_2_control <- qc_data[, qc_indx] %>%
  mutate(
    subject = case_when(
      row_number() >= 1 & row_number() <= 38 ~ "sub1",
      row_number() >= 39 & row_number() <= 78 ~ "sub2",
      row_number() >= 79 ~ "sub3"
    )
  )
qc_2_control[c(34,36),"fber"] <- 35000 # an FBER=-1.0 indicates that there is 
#no signal outside the head mask, here very high quality image. 

subject_colors <-
  c(
    "sub1" = "#D40000",
    "sub2" = "#377EB8",
    "sub3" = "#4DAF4A"
  )

## -------------------------------Plot the correlaiton of IQMs
# Remove the 'subject' column as it's categorical
correlation_matrix  <- dplyr::select(qc_2_control, -subject) %>%
  cor(.)


# Perform hierarchical clustering
hc <- hclust(dist(correlation_matrix))
hc_order <- hc$order

# Reorder the correlation matrix based on the hierarchical clustering
correlation_matrix_ordered <- correlation_matrix[hc_order, hc_order]

# Melt the correlation matrix for visualization
correlation_matrix_melt <- melt(correlation_matrix_ordered)

# Define a reversed RdBu color palette
colors <- rev(colorRampPalette(brewer.pal(n = 11, name = "RdBu"))(100))


# Plotting with ggplot2
plots_qc_matrix <- ggplot(data = correlation_matrix_melt, aes(x = Var1, y = Var2)) +
  geom_tile(aes(fill = value)) +
  geom_text(aes(label = sprintf("%.2f", value)), color = "black", size = 4.5) + # Add correlation coefficients
  scale_fill_gradientn(colors = colors, limits = c(-1, 1), breaks = c(-0.5, 0, 0.5), labels = c("-0.5", "0", "0.5"), name = "") +
  theme_minimal() +
  labs(x = "", y = "", title = "") +
  theme(
    plot.margin = margin(t = 10, r = 10, b = 10, l = 30, unit = "pt"),
    legend.text = element_text(size = 12, family = "Arial", face = "bold", colour = "black"),
    # plot.title = element_text(hjust = 0.5, size = 18, family = "Arial", face = "bold", colour = "black"),
    axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1,
                               size = 16, family = "Arial", face = "bold", colour = "black"),
    axis.text.y = element_text(size = 16, family = "Arial", face = "bold", colour = "black")
  ) +
  coord_fixed()


# Print the plot
plots_qc_matrix


## -------------------------------Plot the distribution of IQMs

plots_qc_half_violin <-  lapply(qc_indx, function(index) {
  qc_data_long <- qc_2_control %>%
    dplyr::select(all_of(index), subject) %>%
    pivot_longer(
      cols = index,
      names_to = "qc_measure",
      values_to = "value"
    )

 ggplot(qc_data_long) +
    geom_violin(
      aes(
        x = subject,
        y = value,
        fill = NULL,
        color = subject
      ),
      trim = FALSE,
      position = position_dodge(0.8),
      draw_quantiles = NULL
    ) +
    geom_dotplot(
      aes(
        x = subject,
        y = value,
        fill = subject,
        color = subject
      ),
      binaxis = "y",
      stackdir = "down",
      position = position_dodge(0.8),
      dotsize = 0.5
    ) +
    scale_fill_manual(values = subject_colors) +
    scale_color_manual(values = subject_colors) +
    theme_minimal() +
    theme(
      plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt"),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_text(
        angle = 30,
        hjust = 1,
        size = 14,
        family = "Arial",
        face = "bold",
        colour = "black"
      ),
      plot.title = element_text(
        hjust = 0.5,
        size = 18,
        family = "Arial",
        face = "bold",
        colour = "black"
      ),
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      panel.grid.minor.y = element_blank()
    ) +
    labs(
      title = index,
      x = NULL,
      y = NULL
    )
})

names(plots_qc_half_violin) <- qc_indx

plots_qc_half_violin$fber <- plots_qc_half_violin$fber +
  scale_y_continuous(trans = 'log10') # the outliers are too large

plots_qc_sub <- plots_qc_half_violin$cnr + plots_qc_half_violin$snr_total + plots_qc_half_violin$wm2max +
  plots_qc_half_violin$cjv + plots_qc_half_violin$efc + plots_qc_half_violin$inu_med +
  plots_qc_half_violin$fber + plots_qc_half_violin$fwhm_avg + plots_qc_half_violin$rpve_gm

#-------------------------load the plots_per_change
load("plots_per_change.RData") # generated from the script 6.1_Total_Brain_CV


# Now combine them together
Fig1 <- (plots_qc_sub| plots_qc_matrix) / plots_per_change

ggsave("/Users/wang/Library/CloudStorage/OneDrive-UniversityofBergen/Desktop/UiB/Manuscripts/BBSC/4.T1w_Image/Figures/Fig1.tiff", 
       plot = Fig1, width = 13, height = 13, units = 'in', dpi = 300, compression = "lzw")

# ggsave("/Users/joeywang/Library/CloudStorage/OneDrive-UniversityofBergen/Desktop/UiB/Manuscripts/BBSC/4.T1w_Image/Figures/Fig1.tiff",
        # plot = Fig1, width = 13, height = 13, units = 'in', dpi = 300, compression = "lzw")




