#bash

# crash crash crash, maybe related with RAM 
docker run -ti --rm \
    -v /Users/wang/Desktop/Research_projects/BBSC/Functional/Data/Reorganized:/data:ro \
    -v /Users/wang/Desktop/Research_projects/BBSC/Functional/fMRI_Prep:/out \
    -v /Users/wang/Desktop/Freesurfer/license.txt:/opt/freesurfer/license.txt \
    nipreps/fmriprep:22.0.1 \
    /data /out \
    participant \
    --skip_bids_validation --participant_label 1 --cifti-output 91k --medial-surface-nan \
    --use-aroma --error-on-aroma-warnings \
    --fs-subjects-dir /Users/wang/Desktop/Research_projects/BBSC/Functional/T1_analysis \
    --output-spaces MNI152NLin6Asym:res-2 MNI152NLin2009cAsym \
    --output-layout bids --low-mem --resource-monitor --verbose

# analysis anatomical data only
docker run -ti --rm \
    -v /Users/wang/Desktop/Research_projects/BBSC/Functional/Data/Reorganized:/data:ro \
    -v /Users/wang/Desktop/Research_projects/BBSC/Functional/fMRI_Prep:/out \
    -v /Users/wang/Desktop/Freesurfer/license.txt:/opt/freesurfer/license.txt \
    nipreps/fmriprep:22.0.1 \
    /data /out participant \
    --skip_bids_validation --participant_label 1 --anat-only \
    --fs-subjects-dir /Users/wang/Desktop/Research_projects/BBSC/Functional/T1_analysis \
    --output-spaces MNI152NLin6Asym:res-2 MNI152NLin2009cAsym \ #can only use one
    --output-layout bids --resource-monitor --verbose

# preprocessing the BOLD data only
docker run -ti --rm \
    -v /Users/wang/Desktop/Research_projects/BBSC/Functional/Data/Reorganized:/data:ro \
    -v /Users/wang/Desktop/Research_projects/BBSC/Functional/fMRI_Prep:/out \
    -v /Users/wang/Desktop/Freesurfer/license.txt:/opt/freesurfer/license.txt \
    nipreps/fmriprep:22.0.1 \
    /data /out participant \
    --skip_bids_validation --participant_label 1 --fs-no-reconall \
    --fs-subjects-dir /Users/wang/Desktop/Research_projects/BBSC/Functional/T1_analysis \
    --output-spaces MNI152NLin6Asym:res-2 MNI152NLin2009cAsym \ #can only use one
    --output-layout bids --resource-monitor --verbose  


# project BOLD to surface template?
docker run -ti --rm \
    -v /Users/wang/Desktop/Research_projects/BBSC/Functional/Data/Reorganized:/data:ro \
    -v /Users/wang/Desktop/Research_projects/BBSC/Functional/fMRI_Prep:/out \
    -v /Users/wang/Desktop/Freesurfer/license.txt:/opt/freesurfer/license.txt \
    nipreps/fmriprep:22.0.1 \
    /data /out participant \
    --skip_bids_validation --participant_label 1 --cifti-output 91k --medial-surface-nan \
    --fs-subjects-dir /Users/wang/Desktop/Research_projects/BBSC/Functional/T1_analysis/sub-1 \
    --output-spaces MNI152NLin6Asym:res-2\
    --output-layout bids --nprocs 6 --resource-monitor --verbose