#bash

# analysis anatomical data only
docker run -ti --rm \
    -v /Users/joeywang/Desktop/BBSC/Functional/fmri/rawdata:/data:ro \
    -v /Users/joeywang/Desktop/BBSC/Functional/fmri/derivatives/fmriprep-23.2.0:/out \
    -v /Users/joeywang/Desktop/license.txt:/opt/freesurfer/license.txt \
    nipreps/fmriprep:23.2.0 \
    /data /out participant \
    --skip_bids_validation --participant_label 1 --anat-only \
    --longitudinal --resource-monitor --verbose

# preprocessing the BOLD data only
docker run -ti --rm \
    -v /Users/joeywang/Desktop/BBSC/Functional/fmri/rawdata:/data:ro \
    -v /Users/joeywang/Desktop/BBSC/Functional/fmri/derivatives/fmriprep-23.2.0:/out \
    -v /Users/joeywang/Desktop/license.txt:/opt/freesurfer/license.txt \
    nipreps/fmriprep:23.2.0 \
    /data /out participant \
    --skip_bids_validation --participant_label 1 \
    --longitudinal --resource-monitor --verbose 
    
    
 # suggestions to use fmriprep:
 # 1. maybe crash because the RAM is small, so limit the process to 4gb per core (using 4 cores even you have 8 cores if the RAM is only 16G)
 # 2. do not be greedy, means if you want all things from fMRIprep (cifiti format, ica, ...), it probabaly crash unless you have a powerful computer;
 #.   so, optimize what you want
