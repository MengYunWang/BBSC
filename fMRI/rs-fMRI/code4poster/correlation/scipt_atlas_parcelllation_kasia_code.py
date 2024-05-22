# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.


"""
#resolution of the data 2
#TR=2

import os
import numpy as np
from nilearn import datasets

#read in text files with session

os.chdir('/home/restate/Documents/BBSC/Functional/Pre_processing/sub-2/')
path='/home/restate/Documents/BBSC/Functional/Pre_processing/sub-2/'

f=open("sessions.txt", "r")
lines=f.read().splitlines()

atlas=datasets.fetch_atlas_schaefer_2018(n_rois=400, yeo_networks=7, 
                                           resolution_mm=2, 
                                           base_url=None, resume=True, verbose=1)

from nilearn.maskers import NiftiLabelsMasker


atlas_filename = atlas.maps
labels = atlas.labels

all_sessions_timecourses=[]
all_sessions_static=[]

from nilearn.connectome import ConnectivityMeasure

for i in range(len(lines)):
    print(i)
    file=os.path.join(path, lines[i],'prep.feat/filtered_func_data_clean_standard.nii.gz')
    conf=os.path.join(path, lines[i], 'prep.feat/mc/prefiltered_func_data_mcf.par')
    
    masker = NiftiLabelsMasker(labels_img=atlas_filename, standardize=True)

    masker.fit(file)
    timeseries=masker.transform(file, confounds=conf)
    
    
    correlation_static = ConnectivityMeasure(kind='correlation')
    static_conn=correlation_static.fit_transform([timeseries])
    #static=np.fill_diagonal(static_conn, 0)
    
    all_sessions_timecourses.append(timeseries)
    all_sessions_static.append(static_conn)
    
w=10
o=5

windows=[]
all_correlations=[]

def frange(start, stop, step=1.0):
    ''' "range()" like function which accept float type'''
    i = start
    while i < stop:
        yield i
        i += step
 

for numb in range(0,40):
    for i in frange(0, 300, 10):

        ts=all_sessions_timecourses[numb]
        strip=ts[i:i+w,:]
               
        correlation_measure = ConnectivityMeasure(kind='correlation')
        correlation_matrix = correlation_measure.fit_transform([strip])
        correlation_matrix=np.reshape(correlation_matrix, [400,400])
        np.fill_diagonal(correlation_matrix, 0)

        windows.append(correlation_matrix)
        
    all_correlations.append(windows)




        
