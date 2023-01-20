#! /bin/bash

export FREESURFER_HOME=/Applications/freesurfer/7.2.0
export SUBJECTS_DIR=/Users/wang/Desktop/BBSC/Functional/anat/sub-1
source $FREESURFER_HOME/SetUpFreeSurfer.sh

cd SUBJECTS_DIR

ls *.nii | parallel --jobs 8 recon-all -s {.} -i {} -all -qcache

