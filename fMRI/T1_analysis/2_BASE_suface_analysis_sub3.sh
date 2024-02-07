#! /bin/bash

export FREESURFER_HOME=/Applications/freesurfer/7.2.0
export SUBJECTS_DIR=/Users/wang/Desktop/Research_projects/BBSC/Functional/T1_analysis/sub-3
source $FREESURFER_HOME/SetUpFreeSurfer.sh

cd SUBJECTS_DIR

recon-all -base sub-3_base -tp sub-3_ses-1_T1w \
                           -tp sub-3_ses-2_T1w \
                           -tp sub-3_ses-3_T1w \
                           -tp sub-3_ses-4_T1w \
                           -tp sub-3_ses-5_T1w \
                           -tp sub-3_ses-6_T1w \
                           -tp sub-3_ses-7_T1w \
                           -tp sub-3_ses-8_T1w \
                           -tp sub-3_ses-9_T1w \
                           -tp sub-3_ses-10_T1w \
                           -tp sub-3_ses-11_T1w \
                           -tp sub-3_ses-12_T1w \
                           -tp sub-3_ses-13_T1w \
                           -tp sub-3_ses-14_T1w \
                           -tp sub-3_ses-15_T1w \
                           -tp sub-3_ses-16_T1w \
                           -tp sub-3_ses-17_T1w \
                           -tp sub-3_ses-18_T1w \
                           -tp sub-3_ses-19_T1w \
                           -tp sub-3_ses-20_T1w \
                           -tp sub-3_ses-21_T1w \
                           -tp sub-3_ses-22_T1w \
                           -tp sub-3_ses-23_T1w \
                           -tp sub-3_ses-24_T1w \
                           -tp sub-3_ses-25_T1w \
                           -all