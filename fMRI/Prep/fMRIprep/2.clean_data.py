# -*- coding: utf-8 -*-

"""
clean the data preprocessed from frmriprep with nilearn

By M.-Y. Wang 23-Nov-2022
"""

import os
from nilearn import image as nimg
from nilearn import plotting as nplot
import matplotlib.pyplot as plt
import numpy as np
import nibabel as nib
import bids
import pandas as pd

fmriprep_dir ='/Users/wang/Desktop/Research_projects/BBSC/Functional/fMRI_Prep'
os.chdir(fmriprep_dir)

#------- extract the preprocessed data from fmriPrep
layout = bids.BIDSLayout(fmriprep_dir,validate=False,
                        config=['bids','derivatives'])

sub_list = layout.get_subjects()
sub_list.sort()

#--- first loop, pull out all the funcitonal data, mask data, and confounds from subjects
for subi,sub_name in enumerate(sub_list):
    func_files = layout.get(subject=sub_list[subi],
                        datatype='func', task='rest',
                        desc='preproc',
                        space='MNI152NLin6Asym',
                        extension='nii.gz',
                       return_type='file')

    mask_files = layout.get(subject=sub_list[subi],
                        datatype='func', task='rest',
                        desc='brain',
                        suffix='mask',
                        space='MNI152NLin6Asym',
                        extension="nii.gz",
                       return_type='file')

    confound_files = layout.get(subject=sub_list[subi],
                            datatype='func', task='rest',
                            desc='confounds',
                           extension="tsv",
                           return_type='file')
    #--- second loop, loop through all sessions from the subject
    for sesi,ses_file in enumerate(func_files):
        func_file = func_files[sesi]
        mask_file = mask_files[sesi]
        confound_file = confound_files[sesi]
        
        # read the confound file and select which confound we would like to use
        # to regress out the confounds from our data
        confound_df = pd.read_csv(confound_file, delimiter='\t')
         # confound_df.head()
         # Select confounds
        confound_vars = ['trans_x','trans_y','trans_z',
                        'rot_x','rot_y','rot_z',
                        'global_signal',
                        'csf', 'white_matter']
         # Get derivative column names
        derivative_columns = ['{}_derivative1'.format(c) for c
                            in confound_vars]
         #print(derivative_columns)
        final_confounds = confound_vars + derivative_columns
         #print(final_confounds)
        confound_df = confound_df[final_confounds]
         #confound_df.head()

         #First we'll load in our data and check the shape
        raw_func_img = nimg.load_img(func_file)
        raw_func_img.shape

        func_img = raw_func_img.slicer[:,:,:,4:]
        func_img.shape

         #Drop confound dummy TRs
        drop_confound_df = confound_df.loc[4:]
         #print(drop_confound_df.shape) #number of rows should match that of the functional image
        drop_confound_df.head()

        confounds_matrix = drop_confound_df.values

         #Confirm matrix size is correct
        confounds_matrix.shape

         #Set some constants
        high_pass= 0.009
        low_pass = 0.08
        t_r = 2.05

         #Clean the data with nimg.clean_img
        clean_img = nimg.clean_img(func_img,confounds=confounds_matrix,detrend=True,standardize=True,
                         low_pass=low_pass,high_pass=high_pass,t_r=t_r, mask_img=mask_file)
                         

        clean_img.to_filename('./sub-%s/ses-%s/func/sub-%s_ses-%s_task-rest_space-MNI152NLin6Asym_res-2_desc-preproc_bold_clean_normalized.nii.gz'\
             % (str(subi+1).zfill(2),str(sesi+1).zfill(2),str(subi+1),str(sesi+1)))


#Let's visualize our result!
#nplot.plot_epi(clean_img.slicer[:,:,:,50])
