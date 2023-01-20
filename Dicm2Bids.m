%% ------Converting the DICM data to BIDS 

% created by M.Y. Wang
% 06-july-2021


%% --------
clear all
clc

cd ('/Users/wang/Desktop/Research_projects/BBSC/Functional/Data')

BBSC_path = '/Users/wang/Desktop/Research_projects/BBSC/Functional/Data';

    
cd ([BBSC_path, '/Raw_data']);
sub_path = dir;
sub_path(1:3,:) = [];

for subi = 1:length (sub_path)
    cd ([BBSC_path, '/Raw_data/', sub_path(subi).name]);

    ses_path = dir;
    ses_path (1:3,:) = [];

    for sesi = 1:length(ses_path)

        %-----anat dicm2nii
        cd ([BBSC_path, '/Raw_data/', sub_path(subi).name, '/',ses_path(sesi).name]);
        str_path = dir ('*3D');
        pf.save_json = getpref('dicm2nii_gui_para', 'save_json', true);
        dicm2nii([BBSC_path, '/Raw_data/',sub_path(subi).name,'/',ses_path(sesi).name,'/',str_path.name],...
            [BBSC_path, '/Reorganized/BIDS/', 'sub-', num2str(subi),'/ses-',num2str(sesi),'/','anat']);

        %------change names
        cd ([BBSC_path, '/Reorganized/BIDS/', 'sub-', num2str(subi),'/ses-',num2str(sesi),'/','anat'])
        temp_file = dir('Sag*');
        movefile (temp_file(1).name, ['sub-', num2str(subi), '_ses-', num2str(sesi),'_T1w','.json'])
        movefile (temp_file(2).name, ['sub-', num2str(subi), '_ses-', num2str(sesi),'_T1w','.nii.gz'])
        delete *.mat

        %----function dicm2nii
        cd ([BBSC_path, '/Raw_data/', sub_path(subi).name, '/',ses_path(sesi).name]);
        fmri_path = dir ('*restate');

        pf.save_json = getpref('dicm2nii_gui_para', 'save_json', true);
        dicm2nii([BBSC_path, '/Raw_data/',sub_path(subi).name,'/',ses_path(sesi).name,'/',fmri_path.name],...
            [BBSC_path, '/Reorganized/BIDS/', 'sub-', num2str(subi),'/ses-',num2str(sesi),'/','func']);

        %------change names
        cd ([BBSC_path, '/Reorganized/BIDS/', 'sub-', num2str(subi),'/ses-',num2str(sesi),'/','func']);
        temp_file = dir('fMRI*');
        movefile (temp_file(1).name, ['sub-', num2str(subi), '_ses-', num2str(sesi),'_task-rest','_bold','.json'])
        movefile (temp_file(2).name, ['sub-', num2str(subi), '_ses-', num2str(sesi),'_task-rest','_bold','.nii.gz'])
        delete *.mat

        %----------add a line in json file
%         jsonText2 = fileread(['sub-', num2str(subi), '_ses-', num2str(sesi),'_task-rest','_bold','.json']);
%          % Convert JSON formatted text to MATLAB data types (3x1 cell array in this example)
%         jsonData = jsondecode(jsonText);
%          % Change HighPrice value in Row 3 from 10000 to 12000
%         jsonData.TaskName = 'rest';
%          % Convert to JSON text
%         jsonText2 = jsonencode(jsonData,'PrettyPrint',true);
%          % Write to a json file
%         fid = fopen('Portfolio2.json', 'w');
%         fprintf(fid, '%s', jsonText2);
%         fclose(fid);
    end

end
