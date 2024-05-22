#!/bin/zsh

# This script will creat folder 'generate_profiles_and_ini_params', 
# from which Step1_2 will generate the group.mat which will be used in the next step 

##########################
# Specify output directory

out_dir='/Users/wang/Desktop/Research_projects/BBSC/Functional/Parcellation' 
cd $out_dir

##########################
# Create generate_profiles and ini params directory



ls -d data/Surf_data/*/ | grep -i sub- > subjList.txt # command grep is trying to filter the input with 'sub-'. -i ignore the high or low cases 

sub_num=0
for sub in `cat subjList.txt`; do
        ls -d $sub/*/ | grep -i ses- > sesList.txt	
        ((sub_num=$sub_num+1)) 
        sub_num=$(printf "%02d" $sub_num)
        
        mkdir -p $out_dir/1_generate_profiles_and_ini_params/data_list/fMRI_list
        mkdir -p $out_dir/1_generate_profiles_and_ini_params/data_list/censor_list

        ses_num=0
    for ses in `cat sesList.txt`; do
        ((ses_num=$ses_num+1)) 
        ses_num=$(printf "%02d" $ses_num)
		
        # fMRI data
		lh_fmri="$out_dir/$ses/lh_sub-${sub_num}_ses-${ses_num}_MNI2fsaverage5_with_NaN.nii.gz"
		echo $lh_fmri > $out_dir/1_generate_profiles_and_ini_params/data_list/fMRI_list/lh_sub${sub_num}_sess${ses_num}.txt

		rh_fmri="$out_dir/$ses/rh_sub-${sub_num}_ses-${ses_num}_MNI2fsaverage5_with_NaN.nii.gz"
        echo $rh_fmri > $out_dir/1_generate_profiles_and_ini_params/data_list/fMRI_list/rh_sub${sub_num}_sess${ses_num}.txt
		
		## censor list
		#censor_file="_DVARS50_motion_outliers.txt"
		
		#echo $censor_file >> $out_dir/generate_profiles_and_ini_params/data_list/censor_list/sub${sub}_sess${sess}.txt
		
	done
done
