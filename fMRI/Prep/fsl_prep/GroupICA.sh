#!/bin/bash
#created by Meng-Yun
#05-Aug-2022


# run group ICA and dual regression

filepath=/Users/wang/Desktop/Research_projects/BBSC/Functional/Pre_processing
cd $filepath

echo "-------------->>>>run group ICA with all data"
ls -1 ./*/*/prep.feat/filtered_func_data_clean_standard.nii.gz >> inputlist_all.txt
# group ICA based on all data including sub1-3
melodic -i inputlist_all.txt -o groupICA \
        -m $FSLDIR/data/standard/MNI152_T1_2mm_brain_mask.nii.gz\
        --report --Oall -d 30 --Ostats --logPower

echo " "

# group ICA based on sessions for each sub
for sub in `cat ./subjList.txt`; do
 mkdir -p ./$sub/ses_groupICA
 ls -1 ./$sub/*/prep.feat/filtered_func_data_clean_standard.nii.gz > ./$sub/inputlist.txt

  for d in 15 30 50 100; do
     mkdir -p ./$sub/ses_groupICA/ICA_$d 
     
     echo "-------------->>>>start to run group ICA on $sub with dimen-$d" 
     melodic -i ./$sub/inputlist.txt -o ./$sub/ses_groupICA/ICA_$d \
             -m $FSLDIR/data/standard/MNI152_T1_2mm_brain_mask.nii.gz\
             --report --Oall -d $d --Ostats --logPower
     dual_regression ./$sub/ses_groupICA/ICA_$d/melodic_IC 1 -1 0 ./$sub/groupICA_back_$d.dr `cat ./$sub/inputlist.txt` # Extract timeseries from ICA --- use group spatial map project back to individual and produce individual time serie
     slices_summary ./$sub/ses_groupICA/ICA_$d/melodic_IC 4 $FSLDIR/data/standard/MNI152_T1_2mm ./$sub/ses_groupICA/ICA_$d.sum -1 # creat a summary img for ICA
     echo "Finished"
     echo " "
    done
done









#------------------------------------------------ dual_regression v0.6

# ***NOTE*** ORDER OF COMMAND-LINE ARGUMENTS IS DIFFERENT FROM PREVIOUS VERSION

# Usage: dual_regression <group_IC_maps> <des_norm> <design.mat> <design.con> <n_perm> [--thr] <output_directory> <input1> <input2> <input3> .........
# e.g.   dual_regression groupICA.gica/groupmelodic.ica/melodic_IC 1 design.mat design.con 500 0 grot `cat groupICA.gica/.filelist`

# <group_IC_maps_4D>            4D image containing spatial IC maps (melodic_IC) from the whole-group ICA analysis
# <des_norm>                    0 or 1 (1 is recommended). Whether to variance-normalise the timecourses used as the stage-2 regressors
# <design.mat>                  Design matrix for final cross-subject modelling with randomise
# <design.con>                  Design contrasts for final cross-subject modelling with randomise
# <n_perm>                      Number of permutations for randomise; set to 1 for just raw tstat output, set to 0 to not run randomise at all.
# [--thr]                       Perform thresholded dual regression to obtain unbiased timeseries for connectomics analyses (e.g., with FSLnets)
# <output_directory>            This directory will be created to hold all output and logfilesg
# <input1> <input2> ...         List all subjects' preprocessed, standard-space 4D datasets

# <design.mat> <design.con>     can be replaced with just
# -1                            for group-mean (one-group t-test) modelling.
# If you need to add other randomise options then edit the line after "EDIT HERE" in the dual_regression script




#  ------------------------------------------------MELODIC (Version 3.15)
#  Multivariate Exploratory Linear Optimised Decomposition into Independent Components

# Author: Christian F. Beckmann 
#  Copyright(c) 2001-2013 University of Oxford

# Usage: 
#  melodic -i <filename> <options>
#  	 	 to run melodic
#  melodic -i <filename> --ICs=melodic_IC --mix=melodic_mix <options>
#  	 	 to run Mixture Model based inference on estimated ICs
#  melodic --help 

# Compulsory arguments (You MUST set one or more of):
# 	-i,--in	input file names (either single file name or comma-separated list or text file)

# Optional arguments (You may optionally specify one or more of):
# 	-o,--outdir	output directory name
# 	--Oall	        output everything
# 	-m,--mask	file name of mask for thresholding
# 	-d,--dim	dimensionality reduction into #num dimensions (default: automatic estimation)
# 	-a,--approach	approach for multi-session/subject data:
# 	   concat	temporally-concatenated group-ICA using MIGP ( default )
# 	   tica  	tensor-ICA
# 	--report	generate Melodic web report
# 	--CIFTI	        input/output as CIFTI (warning: auto-dimensionality estimation for CIFTI data is currently inaccurate)
# 	--vn,--varnorm	switch off variance normalisation
# 	-v,--verbose	switch on diagnostic messages



# 	--nomask	switch off masking
# 	--update_mask	switch off mask updating
# 	--nobet		switch off BET
# 	--bgthreshold	brain / non-brain threshold (only if --nobet selected)

# 	--dimest	use specific dim. estimation technique: lap, bic, mdl, aic, mean (default: lap)
# 	--sep_vn	switch on separate variance nomalisation for each input dataset (off by default)
# 	--disableMigp	switch off MIGP data reduction when using -a concat (full temporal concatenation will be used)
# 	--migpN	Number of internal Eigenmaps
# 	--migp_shuffle	Randomise MIGP file order (default: TRUE)
# 	--migp_factor	Internal Factor of mem-threshold relative to number of Eigenmaps (default: 2)
# 	-n,--numICs	numer of IC's to extract (for deflation approach)
# 	--nl		nonlinearity: gauss, tanh, pow3 (default), pow4
# 	--covarweight	voxel-wise weights for the covariance matrix (e.g. segmentation information)
# 	--eps	minimum error change
# 	--epsS	minimum error change for rank-1 approximation in TICA
# 	--maxit		maximum number of iterations before restart
# 	--maxrestart	maximum number of restarts

# 	--mmthresh	threshold for Mixture Model based inference
# 	--no_mm		switch off mixture modelling on IC maps
 
# 	--ICs		input filename of the IC components file for mixture modelling
# 	--mix		input filename of mixing matrix for mixture modelling / filtering
# 	--smode		input filename of matrix of session modes for report generation
# 	-f,--filter	list of component numbers to remove
 
# 	--bgimage	specify background image for report (default: mean image)
 
# 	--tr		TR in seconds
# 	--logPower	calculate log of power for frequency spectrum

# 	--Tdes	        design matrix across time-domain
# 	--Tcon	        t-contrast matrix across time-domain
# 	--Sdes	        design matrix across subject-domain
# 	--Scon	        t-contrast matrix across subject-domain
# 	--Ounmix	output unmixing matrix
# 	--Ostats	output thresholded maps and probability maps
# 	--Opca		output PCA results
# 	--Owhite	output whitening/dewhitening matrices
# 	--Oorig		output the original ICs
# 	--Omean		output mean volume

# 	-V,--version	prints version information
# 	--copyright	prints copyright information
# 	-h,--help	prints this help message
# 	--debug	        switch on debug messages
# 	--report_maps	control string for spatial map images (see slicer)
# 	--keep_meanvol	do not subtract mean volume


