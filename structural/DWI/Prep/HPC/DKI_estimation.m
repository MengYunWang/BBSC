% Estimation of DKI/DTI parameters 
% Cite original papers: 
% Veraart et al. More Accurate Estimation of Diffusion Tensor Parameters Using Diffusion
% Kurtosis Imaging, MRM 65 (2011): 138-145.
% Fieremans et al., White matter characterization with diffusional kurtosis imaging, 
% Neuroimage 58.1 (2011): 177-188.
%
% Orignal scripts were adopted by
% Ivan I. Maximov, 28/12/2017, Ullevål, Oslo
% Ivan I. Maximov, 24.11.2018, ver 1.1

function DKI_estimation(fbval,fbvec,fname,fmask)

warning off

bval = load(fbval);
d1 = size(bval,1);
if d1 == 1
    bval = bval';
end
bvec = load(fbvec);
d1 = size(bvec,1);
if d1 == 3
    bvec = bvec';
end

grad = horzcat(bvec, bval);

addpath(genpath('/cluster/projects/p33/users/ivanma/MATLAB/matlab'));

[fpath,filename,fext] = fileparts(fname);

if isempty(fpath)
    fpath = '.';
end


f = nifti(fname);
dwi = f.dat(:,:,:,:);
dwi = double(dwi);

f = nifti(fmask);
mask = f.dat(:,:,:);
mask = logical(mask);

[b0, dt] = dki_fit(dwi,grad,mask,[0 1 1]);

[fa, md, rd, ad, fe, mk,  rk, ak, ev1, ev2, ev3] = dki_parameters(dt, mask);

fa(isnan(fa)) = 0;
md(isnan(md)) = 0;
rd(isnan(rd)) = 0;
ad(isnan(ad)) = 0;
mk(isnan(mk)) = 0;
mk(mk < 0) = 0;
mk(mk > 3) = 3;
rk(isnan(rk)) = 0;
rk(rk < 0) = 0;
rk(rk > 3) = 3;
ak(isnan(ak)) = 0;
ak(ak < 0) = 0;
ak(ak > 3) = 3;
ev1(isnan(ev1)) = 0;
ev2(isnan(ev2)) = 0;
ev3(isnan(ev3)) = 0;

[awf, eas, ias] = wmti_parameters(dt, mask);
awf(isnan(awf)) = 0;

AxEAD = eas.de1;
RadEAD = (eas.de2 + eas.de3)./2;
Tort = eas.tort;

AxEAD(isnan(AxEAD)) = 0;
RadEAD(isnan(RadEAD)) = 0;
Tort(isnan(Tort)) = 0;

AxIAD = ias.da1;
RadIAD = (ias.da2 + ias.da3)./2;

AxIAD(isnan(AxIAD)) = 0;
RadIAD(isnan(RadIAD)) = 0;

filesave = [fpath '/' filename];

mname = [fpath '/' filename '_DKI_fa.nii'];
save_metrics(mname,fa,f);

mname = [fpath '/' filename '_DKI_md.nii'];
save_metrics(mname,md,f);

mname = [fpath '/' filename '_DKI_ad.nii'];
save_metrics(mname,ad,f);

mname = [fpath '/' filename '_DKI_rd.nii'];
save_metrics(mname,rd,f);

mname = [fpath '/' filename '_DKI_mk.nii'];
save_metrics(mname,mk,f);

mname = [fpath '/' filename '_DKI_ak.nii'];
save_metrics(mname,ak,f);

mname = [fpath '/' filename '_DKI_rk.nii'];
save_metrics(mname,rk,f);

% WMTI
mname = [fpath '/' filename '_WMTI_awf.nii'];
save_metrics(mname,awf,f);

mname = [fpath '/' filename '_WMTI_axEAD.nii'];
save_metrics(mname,AxEAD,f);

mname = [fpath '/' filename '_WMTI_radEAD.nii'];
save_metrics(mname,RadEAD,f);

mname = [fpath '/' filename '_WMTI_axIAD.nii'];
save_metrics(mname,AxIAD,f);

mname = [fpath '/' filename '_WMTI_radIAD.nii'];
save_metrics(mname,RadIAD,f);

% eigenvalues of DTI
mname = [fpath '/' filename '_DKI_ev1.nii'];
save_metrics(mname,ev1,f);

mname = [fpath '/' filename '_DKI_ev2.nii'];
save_metrics(mname,ev2,f);

mname = [fpath '/' filename '_DKI_ev3.nii'];
save_metrics(mname,ev3,f);

disp('Done');


end
% END

% save metric with a given name using the same info as from input file (x/qform)
function save_metrics(fname,data,X)
    N = nifti;
    N = X;
    N.dat.fname = fname;
    N.dat.dtype = 'FLOAT64-LE';
    create(N);
    
    N.dat(:,:,:) = data;
end












