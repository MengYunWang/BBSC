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

if(~exist(fullfile(out_dir,'profiles','sequence')))
    mkdir(fullfile(out_dir,'profiles','sequence'));
end

seed_mesh = 'fsaverage3';
targ_mesh = 'fsaverage5';


for subi = 1:2

    sub = num2str(subi,'%02d');

    if subi == 1
        sess_id =  1:38;
        sess_id([11 14])=[]; % delete the session 11 and 14
        
        sess_idx = 1:3:length(sess_id);

    else
        sess_id = 1:40;
        sess_id([10, 12:15, 20, 22, 23, 28, 30, 35, 38])=[]; % delete the sessions
        sess_idx = 1:3:length(sess_id);
    end

    for i = 1:length(sess_idx)
        
        if sess_idx(i)+5 <= length(sess_id)
        
            num_sess = sess_id(sess_idx(i):sess_idx(i)+5);
        
        else
        
            continue;
        end


        lh_avg_profile_file = fullfile(out_dir,'profiles','sequence',['lh_' targ_mesh '_roi' seed_mesh '_avg_sub' sub '_seq' num2str(i) '_profile.nii.gz']); %
        rh_avg_profile_file = fullfile(out_dir,'profiles','sequence',['rh_' targ_mesh '_roi' seed_mesh '_avg_sub' sub '_seq' num2str(i) '_profile.nii.gz']); %


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

        lh_avg_profile.vol = lh_avg_profile.vol./num_data;
        rh_avg_profile.vol = rh_avg_profile.vol./num_data;

        MRIwrite(lh_avg_profile,lh_avg_profile_file);
        MRIwrite(rh_avg_profile,rh_avg_profile_file);
    end
end


%% 3. Group or sessions level 
% to compute the parcellation with the yeo2011 clustering algorithm

% define varialbes
num_clusters = '7';
num_initialization = '100';
seed_mesh = 'fsaverage3';
targ_mesh = 'fsaverage5';
out_dir = project_dir;

% make directories
if(~exist(fullfile(out_dir,'group','sequence')))
    mkdir(fullfile(out_dir,'group','sequence'));
end


for subi = 1:2
    sub = num2str(subi,'%02d');
    
    if subi == 1
        num_files = sum(arrayfun(@(x) contains(x.name, 'sub01') && contains(x.name, 'lh_'), dir(fullfile(out_dir,'profiles','sequence'))));
    else
        num_files = sum(arrayfun(@(x) contains(x.name, 'sub02') && contains(x.name, 'lh_'), dir(fullfile(out_dir,'profiles','sequence'))));
    end

    for i = 1:num_files

        output_file = fullfile(out_dir,'group','sequence', ['group_sub' sub '_seq' num2str(i) '_.mat']); % 

        lh_avg_profile_file = fullfile(out_dir,'profiles','sequence',...
            ['lh_' targ_mesh '_roi' seed_mesh '_avg_sub' sub '_seq' num2str(i) '_profile.nii.gz']); % changed names
        rh_avg_profile_file = fullfile(out_dir,'profiles','sequence',...
            ['rh_' targ_mesh '_roi' seed_mesh '_avg_sub' sub '_seq' num2str(i) '_profile.nii.gz']); % changed names

        CBIG_VonmisesSeriesClustering_fix_bessel_randnum_bsxfun(targ_mesh, 'cortex', num_clusters, output_file, ...
            lh_avg_profile_file, rh_avg_profile_file, 0, num_initialization, 0, 100, 1);

        % Reorganize output variables
        if(exist(output_file))
            load(output_file);
            clustered.mtc = mtc;
            clustered.lowerbound = lowerbound;
            clustered.lambda = lambda;
            save(output_file,'lh_labels','rh_labels','clustered');
        else
            error('could not find clustering results group.mat')
        end

        
        % Visualization
        CBIG_DrawSurfaceMaps(lh_labels, rh_labels, targ_mesh, 'inflated');
    end
end

% alternative methods:

% use the parcellation from the cerebellum parcellation
% clustered=CBIG_IndCBM_generate_MSHBM_params('${lh_profile}', '${rh_profile}', lh_labels, rh_labels);

% l/rh_profile is the avg_profiles of all or summer or winter sessions;%
% l/rh_labels is the labels from the Thomas 2011 parcellation; 

% then use the fuction: CBIG_IndCBM_extract_MSHBM_result('${output_dir}/MSHBM/all');
% https://github.com/ThomasYeoLab/CBIG/tree/master/stable_projects/brain_parcellation/Xue2021_IndCerebellum


