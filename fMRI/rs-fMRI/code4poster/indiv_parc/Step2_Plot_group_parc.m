%% To match the color and plot the parcellation map of subject1 and subject2

% 1. we would like to plot the parcellation with the color same as the yeo7
% 2. calculate the dice coefficient between subject 1/2 and yeo7

% M.Y. WANG 14-Jul-2023


%% 0. change the dir

clear all
clc
cd ('/Users/wang/Desktop/Research_projects/BBSC/Functional/Parcellation/1_generate_profiles_and_ini_params/group')

%% 1. change the color so that they have the same order as the yeo7

data_dir = '/Users/wang/Desktop/Research_projects/BBSC/Functional/Parcellation/1_generate_profiles_and_ini_params/group/';

% sub_name = {'light_night/group_sub01_light';'light_night/group_sub01_night';...
%     'mor_eve/group_sub01_mor';'mor_eve/group_sub01_eve'};

% sub_name = {'light_night/group_sub02_light';'light_night/group_sub02_night';...
%     'mor_eve/group_sub02_mor';'mor_eve/group_sub02_eve'};

% sub_name = {'light_night/group_sub03_light';'light_night/group_sub03_night';...
%     'mor_eve/group_sub03_mor';'mor_eve/group_sub03_eve'};

sub_name = {'group_sub01';'group_sub02';'group_sub03'};

ref_par = load('1000subjects_clusters007_ref.mat');

% plot the yeo7 network
CBIG_DrawSurfaceMaps(ref_par.lh_labels, ref_par.rh_labels, 'fsaverage5', 'inflated', 0, 7, ref_par.colors);

for subi = 1: length(sub_name)
    % load data
    data_name = [data_dir sub_name{subi} '.mat'];
    sub_par = load(data_name);

    lh_labels = sub_par.lh_labels;
    rh_labels = sub_par.rh_labels;
    clustered = sub_par.clustered;
    colors = ref_par.colors;

    % chnage the lables so that to match the color and compute the dice
    [assign_dice,lh_labels_new, rh_labels_new, lh_dice, rh_dice] = compute_dice2ref (sub_par.lh_labels, sub_par.rh_labels, ref_par.lh_labels, ref_par.rh_labels);

    if subi == 1
    [lh_labels_new, rh_labels_new, lh_dice, rh_dice]= change_net5to7 (lh_labels_new, rh_labels_new, lh_dice, rh_dice, ref_par);
    else
    end


    CBIG_DrawSurfaceMaps(lh_labels_new, rh_labels_new, 'fsaverage5', 'inflated', 0, 7, ref_par.colors);

    save (data_name, 'lh_labels', 'rh_labels', 'clustered', 'colors','assign_dice', 'lh_dice', 'rh_dice', 'lh_labels_new','rh_labels_new');

end




% sub functions

function [assign_dice, lh_labels_new, rh_labels_new, lh_dice, rh_dice]= compute_dice2ref (lh_labels, rh_labels, ref_lh_labels, ref_rh_labels)

sub_par.lh_labels = lh_labels;
sub_par.rh_labels = rh_labels;
ref_par.lh_labels = ref_lh_labels;
ref_par.rh_labels = ref_rh_labels;

% use the following function to match the labels
[~, assign, ~, ~] = CBIG_HungarianClusterMatch([sub_par.lh_labels; sub_par.rh_labels], [ref_par.lh_labels; ref_par.rh_labels], 1);

sub_par.lh_labels_new = zeros (length (sub_par.lh_labels),1);
sub_par.rh_labels_new = zeros (length (sub_par.rh_labels),1);

for i = 1:length(assign)


    ref_par.lh = double(ref_par.lh_labels==i);
    ref_par.rh = double(ref_par.rh_labels == i);
    
    sub_par.lh = double(sub_par.lh_labels == assign(i));
    sub_par.rh = double(sub_par.rh_labels == assign(i));

    % compute the dice between sub_par and ref_par
    sub_par.dice_lh(i) = dice(ref_par.lh,sub_par.lh);
    sub_par.dice_rh(i) = dice(ref_par.rh,sub_par.rh);
    
    
    % change the labels as the reference 
    indx_l = find (sub_par.lh_labels == assign(i));
    indx_r = find (sub_par.rh_labels == assign(i));
    sub_par.lh_labels_new(indx_l) = i;
    sub_par.rh_labels_new(indx_r) = i;
end

assign_dice = assign;
lh_dice = sub_par.dice_lh;
rh_dice = sub_par.dice_rh;
lh_labels_new = sub_par.lh_labels_new;
rh_labels_new = sub_par.rh_labels_new;

end

function [lh_labels_new, rh_labels_new, lh_dice, rh_dice] = change_net5to7 (lh_labels_new, rh_labels_new, lh_dice, rh_dice, ref_par)

% change the label 5 to label 7
lh_net5_idx = (lh_labels_new == 5);
rh_net5_idx = (rh_labels_new == 5);

lh_labels_new (lh_net5_idx) = 7;
rh_labels_new (rh_net5_idx) = 7;
% change the dice value for net5 to nan
lh_dice(5) = nan;
rh_dice(5) = nan;

% recompute the dice for net7
ref_par.lh = double(ref_par.lh_labels==7);
ref_par.rh = double(ref_par.rh_labels == 7);

sub_par.lh = double(lh_labels_new == 7);
sub_par.rh = double(rh_labels_new == 7);

% compute the dice between sub_par and ref_par
lh_dice(7) = dice(ref_par.lh,sub_par.lh);
rh_dice(7) = dice(ref_par.rh,sub_par.rh);

end




