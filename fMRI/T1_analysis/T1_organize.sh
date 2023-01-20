#! /bin/bash

filepath=/Users/wang/Desktop/Research_projects/BBSC/Functional/

cd $filepath
#mkdir -p ./T1_analysis

ls Data/Reorganized/All | grep -i sub- > T1_analysis/subjList.txt # command grep is trying to filter the input with 'sub-'. -i ignore the high or low cases 


for sub in `cat T1_analysis/subjList.txt`; do 
    mkdir -p T1_analysis/$sub
    ls Data/Reorganized/All/$sub | grep -i ses- > T1_analysis/sesList.txt #list files with 'ses-' in sub directories to sesList.txt 
    for ses in `cat T1_analysis/sesList.txt`; do
        gunzip -k  Data/Reorganized/All/$sub/$ses/anat/*.gz # unpack the .gz file with gunzip and keep the .gz file
        mv Data/Reorganized/All/$sub/$ses/anat/*.nii T1_analysis/$sub/ #-x means extraction, -v visible in terminal, -f file name, -c extract to the folder you want
    done
done