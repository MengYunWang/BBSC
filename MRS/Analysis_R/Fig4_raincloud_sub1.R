# This code is to plot the distribution of metabolites in Sub1

# WANG. 28-April-2023

# Set the working directory to the path where your files are located
setwd("/Users/wang/Desktop/Research_projects/BBSC/MRS/Output_v2_4/QuantifyResults")

# install and load the raincloudplots packages

# if (!require(remotes)) {
#   install.packages("remotes")
# }
# remotes::install_github('jorvlan/raincloudplots')

library(raincloudplots)

# create a vector of row indices for morning and afternoon
morning <- c(1:4, 6, 8, 9, 10, 12, 15, 17, 18, 20, 23, 24, 27, 30, 31, 32, 34, 35, 37)
afternoon <- c(5, 7, 11, 13, 14, 16, 19, 21, 22, 25, 26, 28, 29, 33, 36, 38)

# load the dataset
MRS_water_tc <- read.delim("A_TissCorrWaterScaled_Voxel_1_Basis_1.tsv")

# the following metabolites were analysized
mrs_of_interest <- c("Cr", "PCr", "Cr_PCr", "Glu", "Glu_Gln", "NAA", "NAA_NAAG", "PCh_GPC", "mI")

MRS_water_tc_morning <- MRS_water_tc[morning, mrs_of_interest]
MRS_water_tc_afternoon <- MRS_water_tc[afternoon,mrs_of_interest]

# Add NA values to sub1
n_missing_rows <- nrow(MRS_water_tc_morning) - nrow(MRS_water_tc_afternoon)
temp <- data.frame(matrix(NA, ncol = 9, nrow = n_missing_rows))
names(temp) <- mrs_of_interest
MRS_water_tc_afternoon <- rbind(MRS_water_tc_afternoon, temp)

# raincloud plot the group of tCr
df_2x3 <- data_2x2(
  array_1 = MRS_water_tc_morning[,'Cr_PCr'],
  array_2 = MRS_water_tc_afternoon[,'Cr_PCr'], # afternoon
  array_3 = MRS_water_tc_morning[,'PCr'],
  array_4 = MRS_water_tc_afternoon[,'PCr'], # afternoon
  array_5 = MRS_water_tc_morning[,'Cr'],
  array_6 = MRS_water_tc_afternoon[,'Cr'], # afternoon
  labels = (c('Morning','Afternoon')),
  jit_distance = .05,
  jit_seed = 321) 


raincloud_2x3_vertical <- raincloud_2x3_repmes(
  data = df_2x3,
  colors = (c('dodgerblue', 'darkorange', 'dodgerblue',
              'darkorange', 'dodgerblue', 'darkorange')),
  fills = (c('dodgerblue', 'darkorange', 'dodgerblue',
             'darkorange', 'dodgerblue', 'darkorange')),
  size = 1,
  alpha = .6,
  ort = 'v'
  ) +
  
  scale_x_continuous(breaks=c(1,2,3), labels=c("tCr", "PCr", "Cr"), limits=c(0.8, 4)) +
  xlab("") +
  ylab("") +
  ylim(0,20)+
  theme_bw(base_size = 12) +
  theme(
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 1, color = "black"),
        axis.text = element_text(size = 18, face = "bold", family = "sans", color = "black"),
        axis.title = element_text(size = 20, face = "bold")
        )

raincloud_2x3_vertical

# next 
df_Glu <- data_2x2(
  array_1 = MRS_water_tc_morning[,'Glu_Gln'],
  array_2 = MRS_water_tc_afternoon[,'Glu_Gln'], # afternoon
  array_3 = MRS_water_tc_morning[,'Glu'],
  array_4 = MRS_water_tc_afternoon[,'Glu'], # afternoon
  array_5 = MRS_water_tc_morning[,'Glu'],
  array_6 = MRS_water_tc_afternoon[,'Glu'], # afternoon
  labels = (c('Morning','Afternoon')),
  jit_distance = .05,
  jit_seed = 321) 


raincloud_Glu <- raincloud_2x3_repmes(
  data = df_Glu,
  colors = (c('dodgerblue', 'darkorange', 'dodgerblue',
              'darkorange', 'dodgerblue', 'darkorange')),
  fills = (c('dodgerblue', 'darkorange', 'dodgerblue',
             'darkorange', 'dodgerblue', 'darkorange')),
  size = 1,
  alpha = .6,
  ort = 'v'
) +
  
  scale_x_continuous(breaks=c(1,2,3), labels=c("Glx", "Glu", "Glu"), limits=c(0.8, 4)) +
  xlab("") +
  ylab("") +
  ylim(15,31)+
  theme_bw(base_size = 12) +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(size = 1, color = "black"),
    axis.text = element_text(size = 18, face = "bold", family = "sans", color = "black"),
    axis.title = element_text(size = 20, face = "bold")
  )
raincloud_Glu


# next 
df_NAA <- data_2x2(
  array_1 = MRS_water_tc_morning[,'NAA_NAAG'],
  array_2 = MRS_water_tc_afternoon[,'NAA_NAAG'], # afternoon
  array_3 = MRS_water_tc_morning[,'NAA'],
  array_4 = MRS_water_tc_afternoon[,'NAA'], # afternoon
  array_5 = MRS_water_tc_morning[,'NAA'],
  array_6 = MRS_water_tc_afternoon[,'NAA'], # afternoon
  labels = (c('Morning','Afternoon')),
  jit_distance = .05,
  jit_seed = 321) 


raincloud_NAA <- raincloud_2x3_repmes(
  data = df_NAA,
  colors = (c('dodgerblue', 'darkorange', 'dodgerblue',
              'darkorange', 'dodgerblue', 'darkorange')),
  fills = (c('dodgerblue', 'darkorange', 'dodgerblue',
             'darkorange', 'dodgerblue', 'darkorange')),
  size = 1,
  alpha = .6,
  ort = 'v'
) +
  
  scale_x_continuous(breaks=c(1,2,3), labels=c("tNAA", "NAA", ""), limits=c(0.8, 4)) +
  xlab("") +
  ylab("") +
  ylim(13,20)+
  theme_bw(base_size = 12) +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(size = 1, color = "black"),
    axis.text = element_text(size = 18, face = "bold", family = "sans", color = "black"),
    axis.title = element_text(size = 20, face = "bold")
  )
raincloud_NAA


# next 
df_Cho <- data_2x2(
  array_1 = MRS_water_tc_morning[,'PCh_GPC'],
  array_2 = MRS_water_tc_afternoon[,'PCh_GPC'], # afternoon
  array_3 = MRS_water_tc_morning[,'PCh_GPC'],
  array_4 = MRS_water_tc_afternoon[,'PCh_GPC'], # afternoon
  array_5 = MRS_water_tc_morning[,'PCh_GPC'],
  array_6 = MRS_water_tc_afternoon[,'PCh_GPC'], # afternoon
  labels = (c('Morning','Afternoon')),
  jit_distance = .05,
  jit_seed = 321) 


raincloud_Cho <- raincloud_2x3_repmes(
  data = df_Cho,
  colors = (c('dodgerblue', 'darkorange', 'dodgerblue',
              'darkorange', 'dodgerblue', 'darkorange')),
  fills = (c('dodgerblue', 'darkorange', 'dodgerblue',
             'darkorange', 'dodgerblue', 'darkorange')),
  size = 1,
  alpha = .6,
  ort = 'v'
) +
  
  scale_x_continuous(breaks=c(1,2,3), labels=c("tCho", "", ""), limits=c(0.5, 4)) +
  xlab("") +
  ylab("") +
  ylim(2,5)+
  theme_bw(base_size = 12) +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(size = 1, color = "black"),
    axis.text = element_text(size = 18, face = "bold", family = "sans", color = "black"),
    axis.title = element_text(size = 20, face = "bold")
  )
raincloud_Cho


# next 
df_mI <- data_2x2(
  array_1 = MRS_water_tc_morning[,'mI'],
  array_2 = MRS_water_tc_afternoon[,'mI'], # afternoon
  array_3 = MRS_water_tc_morning[,'mI'],
  array_4 = MRS_water_tc_afternoon[,'mI'], # afternoon
  array_5 = MRS_water_tc_morning[,'mI'],
  array_6 = MRS_water_tc_afternoon[,'mI'], # afternoon
  labels = (c('Morning','Afternoon')),
  jit_distance = .05,
  jit_seed = 321) 


raincloud_mI <- raincloud_2x3_repmes(
  data = df_mI,
  colors = (c('dodgerblue', 'darkorange', 'dodgerblue',
              'darkorange', 'dodgerblue', 'darkorange')),
  fills = (c('dodgerblue', 'darkorange', 'dodgerblue',
             'darkorange', 'dodgerblue', 'darkorange')),
  size = 1,
  alpha = .6,
  ort = 'v'
) +
  
  scale_x_continuous(breaks=c(1,2,3), labels=c("mI", "", ""), limits=c(0.5, 4)) +
  xlab("") +
  ylab("") +
  ylim(8,13)+
  theme_bw(base_size = 12) +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(size = 1, color = "black"),
    axis.text = element_text(size = 18, face = "bold", family = "sans", color = "black"),
    axis.title = element_text(size = 20, face = "bold")
  )
raincloud_mI
