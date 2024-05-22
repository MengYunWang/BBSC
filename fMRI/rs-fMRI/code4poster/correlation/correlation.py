# -*- coding: utf-8 -*-

"""
To explore the time fluctuation of the network

clean the data preprocessed from frmriprep
extract the data using yeo7 network
calculate the correlations between each session with average seesion

By M.-Y. Wang 10-Jul-2023
"""

import os
from nilearn import datasets 
from nilearn import plotting as nplot
from nilearn import image as nimg
from nilearn.maskers import NiftiLabelsMasker
from nilearn.connectome import ConnectivityMeasure
from nilearn.interfaces.fmriprep import load_confounds
import nibabel as nib
import numpy as np
import pandas as pd
import bids

fmriprep_dir ='/Users/wang/Desktop/Research_projects/BBSC/Functional/Prep/fMRIprep'
os.chdir(fmriprep_dir)

#--------- load atlas
atlas_yeo = datasets.fetch_atlas_yeo_2011()
atlas_filepath = atlas_yeo.thin_7
atlas_labels = atlas_yeo.colors_7


#--------- build masker to extract the raw data without normalized
masker = NiftiLabelsMasker(labels_img = atlas_filepath, 
                           standardize=True,
                           detrend=True, 
                           low_pass=0.08, 
                           high_pass=0.009, 
                           t_r=2.05,
                           verbose=3)


#--------- extract the preprocessed data from fmriPrep
layout = bids.BIDSLayout(fmriprep_dir,validate=False,
                        config=['bids','derivatives'])

sub_list = layout.get_subjects()
sub_list.sort()

all_data = []
# first loop, pull out all the funcitonal data, mask data, and confounds from subjects
for subi,sub_name in enumerate(sub_list):
    func_files = layout.get(subject=sub_list[subi],
                        datatype='func', task='rest',
                        desc='preproc',
                        space='MNI152NLin6Asym',
                        extension='nii.gz',
                       return_type='file')
    
    confound_files = layout.get(subject=sub_name,
                            datatype='func', task='rest',
                            desc='confounds',
                           extension="tsv",
                           return_type='file')
    
    # second loop over sessions
    all_ses_data = []
    for sesi,ses_file in enumerate(func_files):
        
        func_file = func_files[sesi]
        confound_file = confound_files[sesi]
        
        # read the confound file and select which confound we would like to use
        # to regress out the confounds from our data
        confound_df = pd.read_csv(confound_file, delimiter='\t')
            # Select confounds
        confound_vars = ['trans_x','trans_y','trans_z',
                        'rot_x','rot_y','rot_z',
                        'global_signal',
                        'csf', 'white_matter']
            # Get derivative column names
        derivative_columns = ['{}_derivative1'.format(c) for c
                            in confound_vars]
        final_confounds = confound_vars + derivative_columns
        confound_df = confound_df[final_confounds]
         
        #Drop the first four datapoints
        raw_func_img = nimg.load_img(func_file)
        func_img = raw_func_img.slicer[:,:,:,4:]

        #Drop confound dummy TRs
        drop_confound_df = confound_df.loc[4:]
        confounds_matrix = drop_confound_df.values

        # To extract the data and clean the data
        time_series = masker.fit_transform(func_file,confounds_matrix)
                                   
        all_ses_data.append(time_series)
    all_data = all_data.append(all_ses_data)


#-------- calculate the correlaiton between each session and the average of all sessions

for i in range(len(all_data)):
    sub_data = all_data[i]
    sub_data= np.array(sub_data) #change the type from list to array so can slice the array 
    
    if i==0:
        session2delete = [10,13] # session 11 and 14 were deleted because sleep for sub1
    else:
        session2delete = [9,11,12,13,14,19,21,22,27,29,34,37] # sessions 10,12-15,20,22,23,28,30,35,38 were deleted for sub2
    sub_all_data = np.delete(sub_all_data, session2delete, axis=0)

    #sub_data [session2delet,:,:] = np.nan    
    sub_data_avg = np.nanmean(sub_data, axis=0) #caculate the average across all sesessions
    #sub_data [session2delet,:,:] = 0
 
 # compute the intra network or inter session variablity of the networks
    correlation_measure = ConnectivityMeasure(kind='correlation')
    labels=np.arange(1, sub_data.shape[0]+1).tolist()
    net_corr = []
    for neti in range(sub_data.shape[2]):
        net_time_series = np.transpose(np.squeeze(sub_data[:,:,neti]))
        correlation_matrix = correlation_measure.fit_transform([net_time_series])[0]
        np.fill_diagonal(correlation_matrix, 0)
        nplot.plot_matrix(correlation_matrix, figure=(6, 6), labels=labels,
                     vmax=0.2, vmin=-0.2, title='Intra-Net_'+ str(neti+1),
                     reorder=False)
 
 #-------- compute the correaltion between sessions and sesseion_avg
    ses_corr = []
    for sesi in range(sub_data.shape[0]):
        df1 = pd.DataFrame(np.squeeze(sub_data[sesi,:,:]))
        df2 = pd.DataFrame(sub_data_avg)
        df_corr = df1.corrwith(df2)
        ses_corr.append(df_corr)
    ses_corr = np.array(ses_corr)
    
    import seaborn as sns
    corr_map=pd.DataFrame(ses_corr,columns=['N1', 'N2', 'N3','N4','N5','N6','N7'])
    sns.violinplot(corr_map)
