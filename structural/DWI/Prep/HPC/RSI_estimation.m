function RSI_estimation(fname,fbvec,fbval,fmask)
% New version of RSI scripts and output
% Ivan I. Maximov, 22.05.19
% ver 1.0
%

addpath(genpath('/cluster/projects/p33/users/ivanma/MATLAB/matlab'));
     
bvecs = load(fbvec);
bvals = load(fbval);

temp = nifti(fname);
data = double(temp.dat(:,:,:,:));

[fpath, fbody, ~] = fileparts(fname);
if isempty(fpath)
    fpath = '.';
end

temp = nifti(fmask);
mask = logical(temp.dat(:,:,:));

tol = 0.3;
nscales = 5;
norder = 4;

rsiout = rsi_process(data, bvecs, bvals, mask, tol, nscales, norder)


% CI normalised
mname = [fpath '/' fbody '_RSI_CInorm.nii'];
save_metrics(mname,rsiout.CInorm,temp);
% CI
mname = [fpath '/' fbody '_RSI_CI.nii'];
save_metrics(mname,rsiout.CI,temp);

mname = [fpath '/' fbody '_RSI_sbeta.nii'];
save_metrics(mname,rsiout.sbeta,temp);

mname = [fpath '/' fbody '_RSI_sbetaNorm.nii'];
save_metrics(mname,rsiout.sbetaNorm,temp);


mname = [fpath '/' fbody '_RSI_normB0L0.nii'];
save_metrics(mname,rsiout.normspectra.B0.L0,temp);
mname = [fpath '/' fbody '_RSI_normB1L0.nii'];
save_metrics(mname,rsiout.normspectra.B1.L0,temp);
mname = [fpath '/' fbody '_RSI_normB2L0.nii'];
save_metrics(mname,rsiout.normspectra.B2.L0,temp);
mname = [fpath '/' fbody '_RSI_normB3L0.nii'];
save_metrics(mname,rsiout.normspectra.B3.L0,temp);
mname = [fpath '/' fbody '_RSI_normB4L0.nii'];
save_metrics(mname,rsiout.normspectra.B4.L0,temp);
mname = [fpath '/' fbody '_RSI_normB5L0.nii'];
save_metrics(mname,rsiout.normspectra.B5.L0,temp);

% ADCfast
mname = [fpath '/' fbody '_RSI_ADC_fast.nii'];
save_metrics(mname,rsiout.ADCfast,temp);

% ADCslow
mname = [fpath '/' fbody '_RSI_ADC_slow.nii'];
save_metrics(mname,rsiout.ADCslow,temp);

% FAfast
mname = [fpath '/' fbody '_RSI_FA_fast.nii'];
save_metrics(mname,rsiout.FAfast,temp);

% FAslow
mname = [fpath '/' fbody '_RSI_FA_slow.nii'];
save_metrics(mname,rsiout.FAslow,temp);

mname = [fpath '/' fbody '_RSI_rD_FA.nii'];
save_metrics(mname,rsiout.rD.FA,temp);

end
% END

function save_metrics(fname,data,X)
    N = nifti;
    N = X;
    N.dat.fname = fname;
    N.dat.dtype = 'FLOAT64-LE';
    create(N);
    
    N.dat(:,:,:) = data;
end
% END