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
    #printf '%s\t' $(ls -d "${sub}"_ses-*_T1w/ | grep -v 'long.sub-*_base' | sort -t- -k3n | sed 's/\/$//') >$filepath/sesList.txt

    aparcstats2table --subjects ${SESS_DIR} \
        --hemi lh \
        --tablefile lh_area_dk.txt

    aparcstats2table --subjects ${SESS_DIR} \
        --hemi rh \
        --tablefile rh_area_dk.txt

    aparcstats2table --subjects ${SESS_DIR} \
        --hemi lh \
        --meas thickness \
        --tablefile lh_thickness_dk.txt

    aparcstats2table --subjects ${SESS_DIR} \
        --hemi rh \
        --meas thickness \
        --tablefile rh_thickness_dk.txt
    
    aparcstats2table --subjects ${SESS_DIR} \
        --hemi lh \
        --meas volume \
        --tablefile lh_volume_dk.txt

    aparcstats2table --subjects ${SESS_DIR} \
        --hemi rh \
        --meas volume \
        --tablefile rh_volume_dk.txt    
    
    
    
    aparcstats2table --subjects ${SESS_DIR} \
        --hemi lh \
        --parc aparc.DKTatlas \
        --tablefile lh_area_dkt.txt

    aparcstats2table --subjects ${SESS_DIR} \
        --hemi rh \
        --parc aparc.DKTatlas \
        --tablefile rh_area_dkt.txt

    aparcstats2table --subjects ${SESS_DIR} \
        --hemi lh \
        --meas thickness \
        --parc aparc.DKTatlas \
        --tablefile lh_thickness_dkt.txt

    aparcstats2table --subjects ${SESS_DIR} \
        --hemi rh \
        --meas thickness \
        --parc aparc.DKTatlas \
        --tablefile rh_thickness_dkt.txt
    
    aparcstats2table --subjects ${SESS_DIR} \
        --hemi lh \
        --meas volume \
        --parc aparc.DKTatlas \
        --tablefile lh_volume_dkt.txt

    aparcstats2table --subjects ${SESS_DIR} \
        --hemi rh \
        --meas volume \
        --parc aparc.DKTatlas \
        --tablefile rh_volume_dkt.txt  




    aparcstats2table --subjects ${SESS_DIR} \
        --hemi lh \
        --parc aparc.a2009s \
        --tablefile lh_area_dest.txt

    aparcstats2table --subjects ${SESS_DIR} \
        --hemi rh \
        --parc aparc.a2009s \
        --tablefile rh_area_dest.txt

    aparcstats2table --subjects ${SESS_DIR} \
        --hemi lh \
        --meas thickness \
        --parc aparc.a2009s \
        --tablefile lh_thickness_dest.txt

    aparcstats2table --subjects ${SESS_DIR} \
        --hemi rh \
        --meas thickness \
        --parc aparc.a2009s \
        --tablefile rh_thickness_dest.txt
    
    aparcstats2table --subjects ${SESS_DIR} \
        --hemi lh \
        --meas volume \
        --parc aparc.a2009s \
        --tablefile lh_volume_dest.txt

    aparcstats2table --subjects ${SESS_DIR} \
        --hemi rh \
        --meas volume \
        --parc aparc.a2009s \
        --tablefile rh_volume_dest.txt             
done
