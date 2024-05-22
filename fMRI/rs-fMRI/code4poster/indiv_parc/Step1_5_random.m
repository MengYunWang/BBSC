%% --- Generate functional connectivity profiles and estimate clusters


% Adapted by M.-Y. Wang, 16-Jul-2023

%%  1.initial parameters

clear all
clc

project_dir ='/Users/wang/Desktop/Research_projects/BBSC/Functional/Parcellation/1_generate_profiles_and_ini_params';
cd (project_dir)

addpath(genpath('/Users/wang/CBIG/stable_projects/brain_parcellation/Kong2019_MSHBM'))


%%  2.group avaraged profiles
% To obtain the group averaged profiles, we will average the profiles across 
% all sessions under profiles folder.

out_dir = project_dir;

if(~exist(fullfile(out_dir,'profiles','random')))
    mkdir(fullfile(out_dir,'profiles','random'));
end


% make directories
if(~exist(fullfile(out_dir,'group','random')))
    mkdir(fullfile(out_dir,'group','random'));
end

% define varialbes
num_clusters = '7';
num_initialization = '100';
seed_mesh = 'fsaverage3';
targ_mesh = 'fsaverage5';

ref_par = load(fullfile(out_dir,'group','1000subjects_clusters007_ref.mat'));
ref_par_lh_label = ref_par.lh_labels;
ref_par_rh_label = ref_par.rh_labels;

for subi = 2
   
    sub = num2str(subi,'%02d');
    sub_par = load (fullfile(out_dir,'group',['group_sub', sub,'.mat']));
    sub_par_lh_labels = sub_par.lh_labels_new;
    sub_par_rh_labels = sub_par.rh_labels_new;

    if subi == 1
        sess_id =  1:38;
        sess_id([11 14])=[]; % delete the session 11 and 14
        [lh_sub1_dice_rand, rh_sub1_dice_rand] = deal(zeros(7, 100, length(sess_id)));
    else
        sess_id = 1:40;
        sess_id([10, 12:15, 20, 22, 23, 28, 30, 35, 38])=[]; % delete the sessions
        [lh_sub2_dice_rand, rh_sub2_dice_rand] = deal(zeros(7, 100, length(sess_id)));
    end



    for sesi = 1:27
        
        num_sess_cell = cell(1, 100);
        for itera_i = 1:100
            sess_idx = randperm (length(sess_id), sesi);

            num_sess_cell{itera_i} = sess_id(sess_idx);
        end

        parfor itera_i = 1:100

            num_sess = num_sess_cell{itera_i};

            avg_profile = read_profile_data (out_dir, sub, num_sess, targ_mesh, seed_mesh);

            [lh_labels, rh_labels] = indiv_parcel(targ_mesh, 'cortex', num_clusters, avg_profile,...
                0, num_initialization, 0, 100, 1);

            [assign_dice,lh_labels_new, rh_labels_new, lh_dice, rh_dice] = compute_dice2ref (lh_labels, rh_labels,...
                ref_par_lh_label, ref_par_rh_label, sub_par_lh_labels, sub_par_rh_labels);

            if subi == 1
                [lh_dice, rh_dice]= change_net5to7 (lh_labels_new, rh_labels_new, lh_dice, rh_dice, sub_par);
                lh_sub1_dice_rand (:,itera_i,sesi)= lh_dice; % network * randoms * sessions number
                rh_sub1_dice_rand (:,itera_i,sesi) = rh_dice;
            else
                lh_sub2_dice_rand (:,itera_i,sesi)= lh_dice; % network * randoms * sessions number
                rh_sub2_dice_rand (:,itera_i,sesi) = rh_dice;
            end

        end

    end  
end

data_name1 = [out_dir,'/group','/random', '/group_sub01_rand.mat'];
data_name2 = [out_dir,'/group','/random', '/group_sub02_rand.mat'];
% 
save (data_name1, 'lh_sub1_dice_rand', 'rh_sub1_dice_rand');
save (data_name2, 'lh_sub2_dice_rand', 'rh_sub2_dice_rand');

rmpath(genpath('/Users/wang/CBIG/stable_projects/brain_parcellation/Kong2019_MSHBM'))




%% 3. sub functions


% tic
% num_sess_cell = cell(1, 100);
% for itera_i = 1:100
%     sess_idx = randperm (length(sess_id), sesi);
% 
%     num_sess_cell{itera_i} = sess_id(sess_idx);
% 
% end
% 
% parfor itera_i = 1:100
%     
%     num_sess = num_sess_cell{itera_i};
%     disp(['sess: ', num2str(num_sess)]);
%     
% %     avg_profile = read_profile_data (out_dir, sub, num_sess, targ_mesh, seed_mesh);
% % 
% %     [lh_labels, rh_labels] = indiv_parcel(targ_mesh, 'cortex', num_clusters, avg_profile,...
% %         0, num_initialization, 0, 100, 1);
% % 
% %     [assign_dice,lh_labels_new, rh_labels_new, lh_dice, rh_dice] = compute_dice2ref (lh_labels, rh_labels,...
% %         ref_par_lh_label, ref_par_rh_label, sub_par_lh_labels, sub_par_rh_labels);
% % 
% %     if subi == 1
% %         [lh_dice, rh_dice]= change_net5to7 (lh_labels_new, rh_labels_new, lh_dice, rh_dice, sub_par);
% %         lh_sub1_dice_rand (:,itera_i,sesi) = lh_dice; % network * randoms * sessions number
% %         rh_sub1_dice_rand (:,itera_i,sesi) = rh_dice;
% %     else
% %         lh_sub2_dice_rand (:,itera_i,sesi) = lh_dice; % network * randoms * sessions number
% %         rh_sub2_dice_rand (:,itera_i,sesi) = rh_dice;
% %     end
% 
% end
% elapsedTime = toc/60;
% 
% disp(['Elapsed time: ', num2str(elapsedTime), ' mins']);



function avg_profile = read_profile_data (out_dir, sub, num_sess, targ_mesh, seed_mesh)

num_data = 0;
for sesi = 1:length(num_sess)

    sess = num2str(num_sess(sesi),'%02d');

    out_profile_dir = fullfile(out_dir,'profiles',['sub' sub],['sess' sess]);

    lh_profile_file = fullfile(out_profile_dir,['lh.sub' sub '_sess' sess '_' targ_mesh '_roi' seed_mesh '.surf2surf_profile.nii.gz']);
    rh_profile_file = fullfile(out_profile_dir,['rh.sub' sub '_sess' sess '_' targ_mesh '_roi' seed_mesh '.surf2surf_profile.nii.gz']);

    if(exist(lh_profile_file) && exist(rh_profile_file))
        num_data = num_data + 1;

        lh_data = MRIread(lh_profile_file);
        rh_data = MRIread(rh_profile_file);
        if(num_data == 1)
            lh_avg_profile = lh_data;
            rh_avg_profile = rh_data;
        else
            lh_avg_profile.vol = lh_avg_profile.vol + lh_data.vol;
            rh_avg_profile.vol = rh_avg_profile.vol + rh_data.vol;
        end
    else
        fprintf('Skip: %s \n',lh_profile_file);
        fprintf('Skip: %s \n',rh_profile_file);

    end
end

fprintf('Total seesions: %d \n',num_data);

lh_vol = lh_avg_profile.vol./num_data;
rh_vol = rh_avg_profile.vol./num_data;

vol_size = size(lh_vol);
lh_vol = reshape(lh_vol, prod(vol_size(1:3)), vol_size(4));
rh_vol = reshape(rh_vol, prod(vol_size(1:3)), vol_size(4));

avg_profile =[lh_vol;rh_vol];
end


function [lh_labels, rh_labels] = indiv_parcel(mesh_name, mask, num_clusters, series, num_smooth, num_tries, normalize, max_iter, no_silhouette)

% CBIG_VonmisesSeriesClustering_fix_bessel_randnum_bsxfun(mesh_name, mask, num_clusters, output_file, profile1, profile2, num_smooth, num_tries, normalize, max_iter)
%
% Von Mises-Fisher clustering on surface data.
%
% Input arguments:
%     - mesh_name     : e.g. 'fsaverage5'
%     - mask          : e.g. 'cortex' if mesh_name is 'fsaverage*';
%                       otherwise, pass in empty string or 'NONE'
%     - num_clusters  : number of clusters
%     - output_file   : output file name
%     - profile1      : group average profile on left hemisphere for data in
%                       fsaverage* space;
%                       or group average profile on entire
%                       cortex for data in fs_LR_32k space.
%     - profile2      : group average profile on right hemisphere for data
%                       in fsaverage* space;
%                       it is not useful for data in fs_LR_32k space. You
%                       can pass in empty string or 'NONE'.
%     - num_smooth    : how many times that smoothing is performed. It is
%                       only used when mesh_name is 'fsaverage*'
%     - num_tries     : number of difference random initialization
%     - normailize    : 0 or 1, whether z-normalization is performed across vertices
%     - max_iter      : maximum number of iterations for one random initialization
%     - no_silhouette : if 1 is passed in for this flag, silhouette step
%                       will be skipped. Default is 0, meaning silhouette
%                       step will be performed when num_clusters>1 and data
%                       are in fsaverage* space.
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

if(~isempty(strfind(mesh_name, 'fsaverage')))
    lambda = 500;            % function direcClus_fix_bessel_bxfun() treats 0 as not specified, lambda will be set as its default (500)
elseif(strcmp(mesh_name, 'fs_LR_32k'))
    lambda = 650;          % For data in fs_LR_32k space, lambda is set to be 650
else
    error('Unknown mesh name.')
end


if(~exist('max_iter', 'var'))
    max_iter = 100;
else
    if(ischar(max_iter))
        max_iter = str2num(max_iter);
    end
end

if(~exist('no_silhouette', 'var'))
    no_silhouette = 0;     % silouette step will be performed
end


% % We (Hesheng, Mert, Jorge and I) decided that random initialization can be same across groups of subjects.
% rand('twister',5489)

if(ischar(normalize))
    normalize = str2num(normalize);
end

if(ischar(num_clusters))
    num_clusters = str2num(num_clusters);
end

if(ischar(num_tries))
    num_tries = str2num(num_tries);
end

if(ischar(num_smooth))
    num_smooth = str2num(num_smooth);
end

% read mask
if(~isempty(strfind(mesh_name, 'fsaverage')))
    lh_avg_mesh = CBIG_ReadNCAvgMesh('lh', mesh_name, 'inflated', mask);
    l1 = find(lh_avg_mesh.MARS_label == 2); lh_num_verts = size(lh_avg_mesh.vertices, 2);
    rh_avg_mesh = CBIG_ReadNCAvgMesh('rh', mesh_name, 'inflated', mask);
    l2 = find(rh_avg_mesh.MARS_label == 2); rh_num_verts = size(rh_avg_mesh.vertices, 2);
    l = [l1 l2+length(lh_avg_mesh.MARS_label)];
else
    % mesh.l: 0 - medial wall, 1 - cortex
    [mesh.v, mesh.l, mesh.ct] = read_annotation(fullfile(getenv('CBIG_CODE_DIR'), 'data', 'templates', 'surface', 'fs_LR_32k', 'label', 'medialwall.annot'));
    lh_num_verts = length(mesh.l) / 2;
    rh_num_verts = lh_num_verts;
    cort_label = mesh.ct.table(2, 5);
    l1 = 1:lh_num_verts;          l1 = l1(mesh.l(l1)==cort_label);
    l2 = 1:rh_num_verts;          l2 = l2(mesh.l(l2+lh_num_verts)==cort_label);
    l = [l1 l2+lh_num_verts];
end

% We (Hesheng, Mert, Jorge and I) decided that random initialization can be same across groups of subjects.
% rand('twister',5489) was moved to after MRIread is because MRIread calls/apps/arch/Linux_x86_64/freesurfer/4.5.0/matlab/load_nifti.m, which calls rand('state', sum(100*clock));
rand('twister',5489)

% smooth, only applied for data in fsaverage* space
if(~isempty(strfind(mesh_name, 'fsaverage')))
    series(1:end/2, :)     = transpose(MARS_AverageData(lh_avg_mesh, transpose(series(1:end/2, :)), 0, num_smooth));
    series(end/2+1:end, :) = transpose(MARS_AverageData(rh_avg_mesh, transpose(series(end/2+1:end, :)), 0, num_smooth));
end

% extract mask voxels series
series = series(l, :);

% remove voxels that are completely not correlated with any rois.
non_zero_corr_index = (sum(series, 2) ~= 0);
series = series(non_zero_corr_index, :);

% znormalize (series assumed to be voxels x subjects or voxels x profile)
if(normalize)
    mean_series = nanmean(series, 1);
    std_series = nanstd(series, 1, 1);
    series = bsxfun(@minus, series, mean_series);
    series = bsxfun(@times, series, 1./(std_series+eps) );
end

% Perform Kmeans
if(num_clusters > 1)
    % Normalize to zero mean across subjects
    series = bsxfun(@minus, series, mean(series, 2) );
    tic, clustered = direcClus_fix_bessel_bsxfun(series, num_clusters, size(series, 2) - 1, num_tries, lambda, 0, 0, 1e-4, 1, max_iter, 1); toc
    cidx = clustered.clusters;
    lambda = clustered.lambda;
    mtc = clustered.mtc;
    p = clustered.p;
    lowerbound = clustered.likelihood(end);

    if(sum(non_zero_corr_index) < length(non_zero_corr_index)) %there are vertices with no correlation
        new_cidx = zeros(length(non_zero_corr_index), 1);
        new_cidx(non_zero_corr_index) = cidx;
        cidx = new_cidx;
    end
else
    cidx = ones(length(l), 1);
    lambda = 0;
    mtc = 0;
    p = 0;
end

% Write Clustering Results
lh_labels = zeros(lh_num_verts, 1);
lh_labels(l1) = cidx(1:length(l1));

rh_labels = zeros(rh_num_verts, 1);
rh_labels(l2) = cidx(length(l1)+1:end);


if(num_clusters > 1 && ~isempty(strfind(mesh_name, 'fsaverage')) && no_silhouette==0)
    tic; s = silhouette(series, cidx(non_zero_corr_index), 'correlation'); toc
    if(sum(non_zero_corr_index) < length(non_zero_corr_index)) %there are vertices with no correlation
        new_s = ones(length(non_zero_corr_index), 1);
        new_s(non_zero_corr_index) = s;
        s = new_s;
    end
else
    s = zeros(length(l), 1);
end
lh_s = ones(lh_num_verts, 1);
lh_s(l1) = s(1:length(l1));

rh_s = ones(rh_num_verts, 1);
rh_s(l2) = s(length(l1)+1:end);

end


function [assign_dice, lh_labels_new, rh_labels_new, lh_dice, rh_dice]= compute_dice2ref (lh_labels, rh_labels, ref_lh_labels, ref_rh_labels, sub_lh_labels, sub_rh_labels)


% use the following function to match the labels
[~, assign, ~, ~] = CBIG_HungarianClusterMatch([lh_labels; rh_labels], [ref_lh_labels; ref_rh_labels], 1);

lh_labels_new = zeros (length (lh_labels),1);
rh_labels_new = zeros (length (rh_labels),1);
[lh_dice, rh_dice] = deal(zeros(7, 1));

for i = 1:length(assign)


    sub_lh = double(sub_lh_labels ==i);
    sub_rh = double(sub_rh_labels == i);
    
    lh = double(lh_labels == assign(i));
    rh = double(rh_labels == assign(i));

    % compute the dice between sub_par and ref_par
    lh_dice(i) = dice(sub_lh,lh);
    rh_dice(i) = dice(sub_rh,rh);
    
    
    % change the labels as the reference 
    indx_l = (lh_labels == assign(i));
    indx_r = (rh_labels == assign(i));
    lh_labels_new(indx_l) = i;
    rh_labels_new(indx_r) = i;
end

assign_dice = assign;

end


function [lh_dice, rh_dice] = change_net5to7 (lh_labels_new, rh_labels_new, lh_dice, rh_dice, ref_par)

% change the label 5 to label 7
lh_net5_idx = (lh_labels_new == 5);
rh_net5_idx = (rh_labels_new == 5);

lh_labels_new (lh_net5_idx) = 7;
rh_labels_new (rh_net5_idx) = 7;
% change the dice value for net5 to nan
lh_dice(5) = nan;
rh_dice(5) = nan;

% recompute the dice for net7
ref_par_lh = double(ref_par.lh_labels_new ==7); %careful, because we used sub_par as the reference
ref_par_rh = double(ref_par.rh_labels_new == 7); % so should use lables_new

sub_par_lh = double(lh_labels_new == 7);
sub_par_rh = double(rh_labels_new == 7);

% compute the dice between sub_par and ref_par
lh_dice(7) = dice(ref_par_lh,sub_par_lh);
rh_dice(7) = dice(ref_par_rh,sub_par_rh);

end
