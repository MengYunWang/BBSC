function NODDI_AMICO(fname,fmask,fbvec,fbval)

disp('Start ...');
disp('NODDI and AMICO estimations.');

[fpath, fbody, fext] = fileparts(fname);
addpath(genpath('/cluster/projects/p33/users/ivanma/MATLAB/matlab'));

global AMICO_data_path CONFIG
global KERNELS bMATRIX
global niiSIGNAL niiMASK
AMICO_data_path = './';

% First is AMICO
try
    disp('AMICO is beginning...');
    if ~exist([AMICO_data_path 'AUX_matrices_lmax=12.mat'],'file')
        AMICO_PrecomputeRotationMatrices();
    end
    
        AMICO_SetSubject('NARCO',fname,fbvec,fbval);
        AMICO_LoadData();
        disp('MATLAB: data are loaded');
        AMICO_SetModel('NODDI');
        AMICO_GenerateKernels(true);
        AMICO_ResampleKernels();
        disp('MATLAB: Kernel generated and resampled');
        %start_spams();
        disp('MATLAB: Fit begins from here');
        AMICO_Fit();
    
    disp('AMICO is done.');
catch
    disp('AMICO failed to finish its estimations');
end



disp('It is done.');
end
% END