# -*- coding: utf-8 -*-

import os
from nilearn import datasets, plotting 
from nilearn.maskers import NiftiLabelsMasker
from nilearn.connectome import ConnectivityMeasure
import numpy as np
import glob

os.chdir('/Users/wang/Desktop/Research_projects/BBSC/Functional/fMRI_Prep')
data_path='/Users/wang/Desktop/Research_projects/BBSC/Functional/fMRI_Prep'

#--------- load atlas
atlas_yeo = datasets.fetch_atlas_yeo_2011()
atlas_filepath = atlas_yeo.thick_7
atlas_labels = atlas_yeo.colors_7

#--------- build masker to extract the raw data without normalized
masker = NiftiLabelsMasker(labels_img=atlas_filepath, memory='nilearn_cache',
                              low_pass=.08, high_pass=.009, verbose=5)

# get the path of the data
files_path = glob.glob("*/*/func/*_task-rest_space-MNI152NLin6Asym_res-2_desc-preproc_bold.nii.gz")
# conf_path = glob.glob("*/*/func/mc/prefiltered_func_data_mcf.par")
# print(len(files_path))
files_path = sorted(files_path)
#conf_path = sorted(conf_path)

sub1_data = [f for f in files_path if 'sub-1' in f]
#sub1_conf = [f for f in conf_path if 'sub-1' in f]
#nr_ses = len(sub1_data)
#print(nr_sub1)

#--------- extract the data
from nilearn.interfaces.fmriprep import load_confounds_strategy

all_ses_data = []
for sesi,data_files in enumerate(sub1_data):
      confounds, sample_mask = load_confounds_strategy(data_files, denoise_strategy="simple",
      motion="basic", global_signal="basic")
      time_series = masker.fit_transform(data_files,
                                   confounds=confounds,
                                   sample_mask=sample_mask)
      all_ses_data.append(time_series)

## denoise the data without regess the global signal
#all_ses_data = []
#for sesi,data_files in enumerate(sub1_data):
#      confounds, sample_mask = load_confounds_strategy(data_files, denoise_strategy="simple",
#      motion="basic")
#      time_series = masker.fit_transform(data_files,
#                                   confounds=confounds,
#                                   sample_mask=sample_mask)
#      all_ses_data.append(time_series)


#-------- concatenate the data and then normalize it and break it into sessions again
all_ses_data_= np.array(all_ses_data) #change the type from list to array so can slice the array
alldata_2D = np.concatenate(all_ses_data) #concatenate the sessions together
from sklearn import preprocessing 
all_ses_data_2d_norm = preprocessing.normalize(all_ses_data_2d) #normalize the data
alldata_3D = np.reshape(all_ses_data_2d_norm,(38,350,7)) # reshape to the 3d again
alldata_avg = np.nanmean(alldata_3D, axis=0) #caculate the average across all sesessions

#-------- compute the correaltion between sessions and sesseion_avg
import pandas as pd
ses_corr = []
for sesi in range(alldata_3D.shape[0]):
      df1 = pd.DataFrame(np.squeeze(alldata_3D[sesi,:,:]))
      df2 = pd.DataFrame(alldata_avg)
      df_corr = df1.corrwith(df2)
      ses_corr.append(df_corr)
ses_corr = np.array(ses_corr)





correlation_measure = ConnectivityMeasure(kind='correlation')
correlation_matrix = correlation_measure.fit_transform([time_series])[0]

# Plot the correlation matrix


# Make a large figure
# Mask the main diagonal for visualization:
np.fill_diagonal(correlation_matrix, 0)
# The labels we have start with the background (0), hence we skip the
# first label
# matrices are ordered for block-like representation
plotting.plot_matrix(correlation_matrix, figure=(10, 8), labels=labels[1:],
                     vmax=0.8, vmin=-0.8, title="Confounds",
                     reorder=True)