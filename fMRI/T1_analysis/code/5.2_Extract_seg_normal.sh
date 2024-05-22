#! /bin/bash

filepath=/Users/wang/Desktop/Research_projects/BBSC/Functional/fmri/derivatives/freesurfer-v7.2.0

cd $filepath

ls | grep -i sub- >$filepath/subjList.txt # command grep is trying to filter the input with 'sub-'. -i ignore the high or low cases

for sub in $(cat $filepath/subjList.txt); do

    export SUBJECTS_DIR=/Users/wang/Desktop/Research_projects/BBSC/Functional/fmri/derivatives/freesurfer-v7.2.0/$sub
    export FREESURFER_HOME=/Applications/freesurfer/7.2.0
    source $FREESURFER_HOME/SetUpFreeSurfer.sh

    cd $SUBJECTS_DIR

    SESS_DIR=$(printf '%s\t' $(ls -d "${sub}"_ses-*_T1w/ | grep -v 'long.sub-*_base' | sort -t- -k3n | sed 's/\/$//'))

    asegstats2table --subjects ${SESS_DIR} \
        --tablefile brain_seg_volume.txt

    asegstats2table --subjects ${SESS_DIR} \
        --meas mean \
        --tablefile brain_seg_mean.txt 
    
          
done
