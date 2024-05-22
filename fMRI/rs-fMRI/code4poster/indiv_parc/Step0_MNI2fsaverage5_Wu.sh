#!/bin/zsh

## This cript is trying to conver the volume data into suface data in fsaverage5
# 1. transfer the data into fsaverage; 2. smooth the data with kernel 6; 
# 3. downsample the data into fsaverge5; 4. add nan to medialwall 

# This function is used to convert the MNI data into surface data using WU method (Wu et al., HBP2018).
# Adapted from the CBIG pipeline 
# created by M.-Y. WANG 20-DEC.-2022

out_dir='/Users/wang/Desktop/Research_projects/BBSC/Functional/Parcellation/data' 
cd $out_dir

proj_mesh=fsaverage
down_mesh=fsaverage5
sm=6
regfile="$FREESURFER_HOME/subjects/fsaverage/mri.2mm/reg.2mm.dat"
MATLAB=/Applications/MATLAB_R2022b.app/bin/matlab

ls -d cleaned_data/*/ | grep -i sub- > subjList.txt # command grep is trying to filter the input with 'sub-'. -i ignore the high or low cases 

sub_num=0
for sub in `cat subjList.txt`; do
        ls -p $sub/ | grep -i ses- > sesList.txt	
        ((sub_num=$sub_num+1)) 
        sub_num_double=$(printf "%02d" $sub_num)
        
        ses_num=0
    for ses in `cat sesList.txt`; do
        ((ses_num=$ses_num+1)) 
        ses_num_double=$(printf "%02d" $ses_num)
        
        file2convert=$sub/$ses
      
        mkdir -p Surf_data/sub-${sub_num_double}/ses-${ses_num_double}

        for hemi in lh rh; do
            
            surf_name1=$out_dir/Surf_data/sub-${sub_num_double}/ses-${ses_num_double}/${hemi}_sub-${sub_num_double}_ses-${ses_num_double}_MNI2fsaverage7
            
         ## project data to proj_mesh
            # default method to projection
            #mri_vol2surf --mov $file2convert --reg $regfile --hemi $hemi --projfrac 0.5 --trgsubject $proj_mesh --o $surf_name1.nii.gz --reshape --interp trilinear
            
            # Wu method HBP2018
            matlab_cmd=("addpath(genpath('/Users/wang/Documents/MATLAB/toolbox/MNI2fsaverage'));" \
            "MNI2fsaverage_wang '$hemi' '$file2convert' '$surf_name1.nii.gz';exit;")
            $MATLAB -nodesktop -nodisplay -nosplash -r "$matlab_cmd" |& tee -a $LF

         ## smooth (func from freesurfer)
            mri_surf2surf --hemi $hemi --s $proj_mesh --sval $surf_name1.nii.gz --cortex --fwhm-trg $sm --tval ${surf_name1}_sm${sm}.nii.gz --reshape
 
         ## mri_surf2surf with --cortex flag will make medial wall vertices to be zeros
            # add removed medial wall values back
            matlab_cmd=("addpath('/Users/wang/CBIG/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/utilities');" \
            "CBIG_preproc_fsaverage_medialwall_fillin '$hemi' $proj_mesh '$surf_name1.nii.gz' '${surf_name1}_sm${sm}.nii.gz' '${surf_name1}_sm${sm}_fillin_medialwall.nii.gz'; exit;")
			   $MATLAB -nodesktop -nodisplay -nosplash -r "$matlab_cmd" |& tee -a $LF

         ## downsample to fsaverage5
            surf_name2=$out_dir/Surf_data/sub-${sub_num_double}/ses-${ses_num_double}/${hemi}_sub-${sub_num_double}_ses-${ses_num_double}_MNI2fsaverage5
            mri_surf2surf --hemi $hemi --srcsubject $proj_mesh \
            --sval ${surf_name1}_sm${sm}_fillin_medialwall.nii.gz --nsmooth-in 1 \
            --trgsubject $down_mesh --tval $surf_name2.nii.gz --reshape


         ## Set medial values to be NaN
            matlab_cmd=("addpath('/Users/wang/CBIG/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/utilities');" \
            "CBIG_preproc_set_medialwall_NaN '$hemi' '$down_mesh' '$surf_name2.nii.gz'\
             '${surf_name2}_with_NaN.nii.gz'; exit;")
			   $MATLAB -nodesktop -nodisplay -nosplash -r "$matlab_cmd" |& tee -a $LF
        done
    done
done