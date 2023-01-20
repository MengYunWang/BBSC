
#!/bin/bash
# created by Meng-Yun 
# 29-July-2022

# run the preprocessing on mutisubjects with FSL

filepath=/Users/wang/Desktop/Research_projects/BBSC/Functional
cd $filepath
# mkdir ./Pre_processing

ls ./Data/Reorganized/All | grep -i sub- > ./Pre_processing/subjList.txt # command grep is trying to filter the input with 'sub-'. -i ignore the high or low cases 

for sub in `cat ./Pre_processing/subjList.txt`; do

    mkdir -p ./Pre_processing/$sub
    ls ./Data/Reorganized/All/$sub > ./Pre_processing/$sub/sesList.txt #list files with 'ses-' in sub directories to sesList.txt 
    # ses_num=0
    
    for ses_name in `cat ./Pre_processing/$sub/sesList.txt`; do
     
     echo "===> Starting processing of $sub $ses_name"
     #((ses_num=$ses_num+1)) 
     mkdir -p ./Pre_processing/$sub/$ses_name
     
     echo "Brain extraction..." 
     bet2 $filepath/Data/Reorganized/All/$sub/$ses_name/anat/${sub}_${ses_name}_T1w \
            $filepath/Pre_processing/$sub/$ses_name/${sub}_${ses_name}_T1w_brain  -f 0.5 -g 0 -c 93 119 155 
        # If the brain mask doesn’t exist, create it
        #if [ ! -f anat/${subj}_T1w_brain.nii.gz ]; then
        #    echo "Skull-stripped brain not found, creating it..."
        #    bet2 $filepath/All/sub-1/ses-1/anat/sub-1_ses-1_T1w  \
        #       $filepath/All/sub-1/ses-1/anat/sub-1_ses-1_T1w_brain  -f 0.8 -g 0 -c 93 119 155
        #fi
     
     # Note that we are using the | character to delimit the patterns
     # instead of the usual / character because there are / characters in the pattern.
     cp ./Pre_processing/rs_prep_origin.fsf ./Pre_processing/rs_prep.fsf
     sed -i '' "s|sub-1|${sub}|g" ./Pre_processing/rs_prep.fsf
     sed -i '' "s|ses-1|${ses_name}|g" ./Pre_processing/rs_prep.fsf

     # Now everything is set up to run feat
     echo "===> Starting to run feat"
     feat ./Pre_processing/rs_prep.fsf
     echo "===> Finished! Now next..."
    done
done