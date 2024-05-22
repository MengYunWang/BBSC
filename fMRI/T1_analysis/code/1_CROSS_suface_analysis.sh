#! /bin/bash
filepath=/Users/wang/Desktop/Research_projects/BBSC/Functional/T1_analysis/

cd $filepath

ls | grep -i sub- >$filepath/subjList.txt # command grep is trying to filter the input with 'sub-'. -i ignore the high or low cases

for sub in $(cat $filepath/subjList.txt); do
    export FREESURFER_HOME=/Applications/freesurfer/7.2.0
    export SUBJECTS_DIR=/Users/wang/Desktop/BBSC/Functional/anat/$sub
    source $FREESURFER_HOME/SetUpFreeSurfer.sh

    cd $SUBJECTS_DIR

    ls *.nii | parallel --jobs 8 recon-all -s {.} -i {} -all -qcache
done
