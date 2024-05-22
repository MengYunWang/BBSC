%% --- Generate functional connectivity profiles and estimate initialization parameters


% Adapted by M.-Y. Wang, 16-Nov.-2022

%%  1.initial parameters

clear all
clc

project_dir ='/Users/wang/Desktop/Research_projects/BBSC/Functional/Parcellation/1_generate_profiles_and_ini_params';
cd (project_dir)

addpath(genpath('/Users/wang/CBIG/stable_projects/brain_parcellation/Kong2019_MSHBM'))

%%  2.generate profiles

sub_surf = dir(fullfile('../data/Surf_data'));
sub_surf(1:3)=[];

for subi = 1:3 
    
    ses_surf = dir(fullfile(sub_surf(subi).folder,sub_surf(subi).name));
    ses_surf(1:2)=[];

    for sesi = 1:length(ses_surf)

        CBIG_MSHBM_generate_profiles('fsaverage3','fsaverage5',project_dir,num2str(subi,'%02d'), num2str(sesi,'%02d'),'0');
    end
end

%%  3.group avaraged profiles
% To obtain the group averaged profiles, we will average the profiles across 
% all sessions under profiles folder.

out_dir = project_dir;

if(~exist(fullfile(out_dir,'profiles','avg_profile')))
    mkdir(fullfile(out_dir,'profiles','avg_profile'));
end

seed_mesh = 'fsaverage3';
targ_mesh = 'fsaverage5';

% compute the average profile across all subjects and sessions

lh_avg_profile_file = fullfile(out_dir,'profiles','avg_profile',['lh_' targ_mesh '_roi' seed_mesh '_avg_profile.nii.gz']); % change names to may or dec
rh_avg_profile_file = fullfile(out_dir,'profiles','avg_profile',['rh_' targ_mesh '_roi' seed_mesh '_avg_profile.nii.gz']); % change names

num_data = 0;

for subi = 1:2
    
    sub = num2str(subi,'%02d');
    
    if subi == 1
        num_sess = 1:38;
        num_sess([11 14])=[]; % delete the session 11 and 14
    elseif subi == 2
        num_sess = 1:40;
        num_sess([10, 12:15, 20, 22, 23, 28, 30, 35, 38])=[]; % delete the sessions
%     else 
%         num_sess = 1:25;
%         num_sess([1:2, 4:7, 9, 12, 14, 16:18, 21, 25])=[]; % delete the sessions
    end
    
    % compute the average profile for each subject 

    %lh_avg_profile_file = fullfile(out_dir,'profiles','avg_profile',['lh_' targ_mesh '_roi' seed_mesh '_avg_sub' sub '_profile.nii.gz']); % change names to may or dec
    %rh_avg_profile_file = fullfile(out_dir,'profiles','avg_profile',['rh_' targ_mesh '_roi' seed_mesh '_avg_sub' sub '_profile.nii.gz']); % change names

    %num_data = 0;
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
    
%     fprintf('Total seesions: %d \n',num_data);
% 
%     lh_avg_profile.vol = lh_avg_profile.vol./num_data;
%     rh_avg_profile.vol = rh_avg_profile.vol./num_data;
% 
%     MRIwrite(lh_avg_profile,lh_avg_profile_file);
%     MRIwrite(rh_avg_profile,rh_avg_profile_file);
end

    fprintf('Total seesions: %d \n',num_data);

    lh_avg_profile.vol = lh_avg_profile.vol./num_data;
    rh_avg_profile.vol = rh_avg_profile.vol./num_data;

    MRIwrite(lh_avg_profile,lh_avg_profile_file);
    MRIwrite(rh_avg_profile,rh_avg_profile_file);

%% 4. Group or sessions level 
% to compute the parcellation with the yeo2011 clustering algorithm

% define varialbes
num_clusters = '7';
num_initialization = '100';
seed_mesh = 'fsaverage3';
targ_mesh = 'fsaverage5';
out_dir = project_dir;

% make directories
if(~exist(fullfile(out_dir,'group')))
    mkdir(fullfile(out_dir,'group'));
end

% combing 2 subjects together
CBIG_MSHBM_generate_ini_params(seed_mesh, targ_mesh, num_clusters, num_initialization, out_dir)

% % calculate subject seperately
% for subi = 1:2
%     
%     sub = num2str(subi,'%02d');
%     output_file = fullfile(out_dir,'group',['group_sub' sub '.mat']); % change names according summer or winter
% 
%     lh_avg_profile_file = fullfile(out_dir,'profiles','avg_profile',...
%         ['lh_' targ_mesh '_roi' seed_mesh '_avg_sub' sub '_profile.nii.gz']); % changed names
%     rh_avg_profile_file = fullfile(out_dir,'profiles','avg_profile',...
%         ['rh_' targ_mesh '_roi' seed_mesh '_avg_sub' sub '_profile.nii.gz']); % changed names
% 
%     CBIG_VonmisesSeriesClustering_fix_bessel_randnum_bsxfun(targ_mesh, 'cortex', num_clusters, output_file, ...
%         lh_avg_profile_file, rh_avg_profile_file, 0, num_initialization, 0, 100, 1);
% 
%     % Reorganize output variables
%     if(exist(output_file))
%         load(output_file);
%         clustered.mtc = mtc;
%         clustered.lowerbound = lowerbound;
%         clustered.lambda = lambda;
%         save(output_file,'lh_labels','rh_labels','clustered');
%     else
%         error('could not find clustering results group.mat')
%     end
%     
% %     % Visualization
% %     CBIG_DrawSurfaceMaps(lh_labels, rh_labels, targ_mesh, 'inflated');
% end

% alternative methods:

% use the parcellation from the cerebellum parcellation
% clustered=CBIG_IndCBM_generate_MSHBM_params('${lh_profile}', '${rh_profile}', lh_labels, rh_labels);

% l/rh_profile is the avg_profiles of all or summer or winter sessions;%
% l/rh_labels is the labels from the Thomas 2011 parcellation; 

% then use the fuction: CBIG_IndCBM_extract_MSHBM_result('${output_dir}/MSHBM/all');
% https://github.com/ThomasYeoLab/CBIG/tree/master/stable_projects/brain_parcellation/Xue2021_IndCerebellum


