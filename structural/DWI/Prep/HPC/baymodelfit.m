function baymodelfit(fpath)
% Code is derived from Marco's demo script
% INPUT: path to unziped nifti file and bval/bvec files
% OUTPUT: 11 scalar maps derived from Bayesian model: intra/extra axonal diffusivities,
%         signal fractions for intra/extra axonal water and CSF
%         local/micro DTI metrics derived from simplifications with
%         undirectional metrics (pseudo tensor with three eigenvecs).
%
% version 1.0. Ivan Maximov, 09/12/2021
%

disp('Starting ...');

addpath(genpath('/cluster/projects/p33/users/ivanma/MATLAB'));

disp(['Processing subject: ' fpath]);
fname = [fpath '/data.nii'];
fbvec = [fpath '/data.bvec'];
fbval = [fpath '/data.bval'];

try
    d = bd_loadData_nii(fname,{fbvec,fbval});
catch
    error(['Could not run loadData_nii for ' fpath]);
end

f = nifti([fpath '/brain_mask.nii']);
mask = f.dat(:,:,:);

disp('Data are loaded and sent to fitting function.');

[vols model] = baydiff(d,'noise',20);

disp('Fitting is done. Saving ...');
fmap = {'Dax_intra','Dax_extra','Drad_extra','v_intra','v_csf','v_extra','microADC','microAx','microRd','microFA','turt'};

N = length(fmap);
for i = 1:N
   disp(['Saving: ' fmap{i}]); 
   hdl = nifti;
   hdl = f;
   hdl.dat.fname = [fpath '/data_Bayes_' fmap{i} '.nii'];
   hdl.dat.dtype = 'FLOAT64-LE';
   create(hdl);
   temp = squeeze(vols(:,:,:,i));
   temp(temp>3.1) = 0;
   temp(temp<0) = 0;
   hdl.dat(:,:,:) = temp;
end

disp('Done');
end
% END
