%%    Project individual parcellation to the MNI space

% 1. upsample the lables in fs5 to fsaverge

% 2. use the function (Wu et al. 2018 HBM) to project the surface to mni (Wu et al. 2018 HBM)

%% 
clear all
clc
label_path = '/Users/wang/Desktop/Research_projects/BBSC/Functional/Parcellation/1_generate_profiles_and_ini_params/group';
cd (label_path)

%% 

sub_name = {'group_sub01';'group_sub02'};

for subi = 1: length(sub_name)
    % load data
    data_name = [sub_name{subi} '.mat'];
    sub_par = load(data_name);

    % load the fsaverage data
    lh_mesh7 = CBIG_ReadNCAvgMesh('lh', 'fsaverage', 'white', 'cortex');
    rh_mesh7 = CBIG_ReadNCAvgMesh('rh', 'fsaverage', 'white', 'cortex');
    lh_mesh5 = CBIG_ReadNCAvgMesh('lh', 'fsaverage5', 'white', 'cortex');
    rh_mesh5 = CBIG_ReadNCAvgMesh('rh', 'fsaverage5', 'white', 'cortex');
    
    % upsampling
    lh_labels7 = MARS_NNInterpolate_kdTree(lh_mesh7.vertices, lh_mesh5, sub_par.lh_labels_new');
    rh_labels7 = MARS_NNInterpolate_kdTree(rh_mesh7.vertices, rh_mesh5, sub_par.rh_labels_new');

    output = CBIG_Projectfsaverage2MNI_Ants(lh_labels7', rh_labels7');
    output.vol = int8(output.vol);
    MRIwrite(output, ['sub' num2str(subi,'%02d') 'indi_network.nii.gz']);

end
