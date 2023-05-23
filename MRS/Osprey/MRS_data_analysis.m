% MRS data analysis





%%

clear all
clc

addpath(genpath('/Users/wang/Documents/MATLAB/toolbox/osprey'))
addpath('/Users/wang/Documents/MATLAB/toolbox/spm12')

MRSCont = OspreyJob('/Users/wang/Desktop/Research_projects/BBSC/Analysis_code/BBSC_code/MRS_analysis/BBSC_GE_v2_4_0.m');
MRSCont = OspreyLoad(MRSCont);
MRSCont = OspreyProcess(MRSCont);
MRSCont = OspreyFit(MRSCont);

MRSCont = OspreyCoreg(MRSCont);
MRSCont = OspreySeg(MRSCont);

MRSCont = OspreyQuantify(MRSCont);
MRSCont = OspreyOverview(MRSCont);


nii_viewer({'/Users/wang/Desktop/Research_projects/BBSC/MRS/Output_v2_4/VoxelMasks/ses-25_P35_space-scanner_mask_VoxelOverlap.nii.gz'})
out = osp_plotVoxelOverlap(MRSCont,[2 32 29])



rmpath(genpath('/Users/wang/Documents/MATLAB/toolbox/osprey'))
rmpath(genpath('/Users/wang/Documents/MATLAB/toolbox/spm12'))