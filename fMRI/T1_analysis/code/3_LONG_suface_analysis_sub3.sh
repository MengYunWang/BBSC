#! /bin/bash

export SUBJECTS_DIR=/Users/wang/Desktop/Research_projects/BBSC/Functional/T1_analysis/sub-3
export FREESURFER_HOME=/Applications/freesurfer/7.2.0
source $FREESURFER_HOME/SetUpFreeSurfer.sh

cd $SUBJECTS_DIR
ls -d */ | grep -i ses- >$SUBJECTS_DIR/sesList.txt #list files with 'ses-' in sub directories to sesList.txt

for ses in $(cat $SUBJECTS_DIR/sesList.txt); do
    recon-all -long $ses sub-3_base -all
done
