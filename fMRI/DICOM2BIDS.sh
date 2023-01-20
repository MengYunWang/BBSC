#!/bin/bash
# Created by M.Y. WANG
# 31-03-2022

# Organize the data into BIDS format

filepath=/Users/wang/Desktop/Research_projects/BBSC/Functional/Data
cd $filepath
mkdir ./Reorganized/All

ls ./Raw_data | grep -i sub- > ./Raw_data/subjList.txt # command grep is trying to filter the input with 'sub-'. -i ignore the high or low cases 

for sub in `cat ./Raw_data/subjList.txt`; do 
    mkdir ./Reorganized/All/$sub
    ls ./Raw_data/$sub > ./Raw_data/sesList.txt #list files with 'ses-' in sub directories to sesList.txt 
    ses_num=0
    for ses_name in `cat ./Raw_data/sesList.txt`; do
        
        ((ses_num=$ses_num+1)) 
        # ses_num=$(printf "%02d" $ses_num)
        # # ses_num=ses_num+1 | printf "%02d " 
        sespath=./Reorganized/All/$sub/ses-${ses_num}
        
        mkdir $sespath
        mkdir $sespath/anat
        mkdir $sespath/func


        dcm2niix -b y -z y -f ${sub}_ses-${ses_num}_T1w -o $sespath/anat/ ./Raw_data/$sub/$ses_name/*_Sag_FSPGR_3D/ 

        dcm2niix -b y -z y -f ${sub}_ses-${ses_num}_task-rest_bold -o $sespath/func/ ./Raw_data/$sub/$ses_name/*_fMRI_default_restate/
        cd ./$sespath/func/   
        jq '{"TaskName": "rest"} + .' ${sub}_ses-${ses_num}_task-rest_bold.json > temp.json
        rm ${sub}_ses-${ses_num}_task-rest_bold.json
        mv temp.json ${sub}_ses-${ses_num}_task-rest_bold.json  

        cd $filepath
    done
    
done

