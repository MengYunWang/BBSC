#! /bin/bash

filepath=/Users/wang/Desktop/Research_projects/BBSC/Functional/T1_analysis/

cd $filepath

ls | grep -i sub- >$filepath/subjList.txt # command grep is trying to filter the input with 'sub-'. -i ignore the high or low cases

for sub in $(cat $filepath/subjList.txt); do
    
    export SUBJECTS_DIR=/Users/wang/Desktop/Research_projects/BBSC/Functional/T1_analysis/$sub
    export FREESURFER_HOME=/Applications/freesurfer/7.2.0
    source $FREESURFER_HOME/SetUpFreeSurfer.sh

    cd $SUBJECTS_DIR
    ls -d */| grep -i ses- >$filepath/sesList.txt #list files with 'ses-' in sub directories to sesList.txt

    for ses in $(cat $filepath/sesList.txt); do
        recon-all -long $ses ${sub}_base -all
    done
done
