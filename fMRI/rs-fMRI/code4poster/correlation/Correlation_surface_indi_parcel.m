%% compute the network correlation based on individual parcellation



% Created by M.Y. Wang 31-Jul-2023 

%% define parameters

clear all
clc
datapath = '/Users/wang/Desktop/Research_projects/BBSC/Functional/Parcellation/data/Surf_data/';
parc_path = '/Users/wang/Desktop/Research_projects/BBSC/Functional/Parcellation/1_generate_profiles_and_ini_params/group/';
cd (datapath)

%% extract data and compute correlation on surface data


sub_surf = dir(datapath);
sub_surf(1:3)=[];




data_network = [];

for subi = 1:2

    ses_surf = dir(fullfile(sub_surf(subi).folder,sub_surf(subi).name));
    ses_surf(1:2)=[];

    indi_parc = load ([parc_path 'group_sub' num2str(subi,'%02d') '.mat']);

    for sesi = 1:length(ses_surf)

        surf_data_path = fullfile(ses_surf(sesi).folder,ses_surf(sesi).name);
        % read and extract the data
        [lh_surf_data, rh_surf_data] = read_surf_data(surf_data_path);
        data_network (:,:,sesi,subi) = extract_data_yeo7(lh_surf_data, rh_surf_data, indi_parc); % network * datapoints * sessions * sub

    end
    
    % compute the inter-session correlation
    data2compute = squeeze(data_network(:,:,:,subi));
    
    if subi == 1
        inter_session_corr_sub1 = corr_inter_session (data2compute);
    else
        inter_session_corr_sub2 = corr_inter_session (data2compute);
    end
    
    % compute the session and average session correlation
    session_avg = mean (data2compute, 3);
    
    if subi == 1
        session_avg_corr_sub1 = corr_avg_session (data2compute,session_avg);
    else
        session_avg_corr_sub2 = corr_avg_session (data2compute,session_avg);

    end

end

datasave_path = '/Users/wang/Desktop/Research_projects/BBSC/Functional/Parcellation/correlation';
if (~exist(datasave_path, 'dir'))
    mkdir(datasave_path) 
else 
end

cd (datasave_path)
% save the data
save (fullfile(datasave_path, 'sub1_correlation_indi.mat'), 'inter_session_corr_sub1', 'session_avg_corr_sub1', '-mat');
save (fullfile(datasave_path, 'sub2_correlation_indi.mat'), 'inter_session_corr_sub2', 'session_avg_corr_sub2', '-mat');



%% sub-functions

function [lh_surf_data, rh_surf_data] = read_surf_data (surf_data_path)

lh_fmri = MRIread([surf_data_path '/lh*5_with_NaN.nii.gz']);
rh_fmri = MRIread([surf_data_path '/rh*5_with_NaN.nii.gz']);

vol_size = size(lh_fmri.vol);
lh_surf_data = reshape(lh_fmri.vol, prod(vol_size(1:3)), vol_size(4));
rh_surf_data = reshape(rh_fmri.vol, prod(vol_size(1:3)), vol_size(4));


end

function network_data = extract_data_yeo7 (lh_surf_data, rh_surf_data, parc)

network_data = zeros(7,size (lh_surf_data,2));

for neti = 1:7
    lh_network_indx = (parc.lh_labels_new == neti);
    lh_network_data = lh_surf_data(lh_network_indx,:);

    rh_network_indx = (parc.rh_labels_new == neti);
    rh_network_data = rh_surf_data(rh_network_indx,:);

    data_network = vertcat(lh_network_data, rh_network_data);
    network_data (neti,:)= mean (data_network);
end
end



function inter_session_corr = corr_inter_session (data2compute)

inter_session_corr = zeros(size(data2compute,3), size(data2compute,3), size(data2compute,1));

for neti = 1:7
    data = squeeze (data2compute (neti,:,:));

    inter_session_corr (:,:,neti)= corrcoef (data);

end
end

function   session_avg_corr = corr_avg_session (data2compute, session_avg)

session_avg_corr = zeros (size(data2compute,3),size(data2compute,1));

for neti = 1:7
    data = squeeze (data2compute (neti,:,:));
    data_avg = session_avg (neti,:)';
    data_avg_matrix = repmat (data_avg, 1, size (data,2));
    session_corr = corr (data, data_avg_matrix);
    session_avg_corr (:,neti)= session_corr(:,1);

end
end






