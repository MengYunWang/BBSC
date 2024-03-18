% Original article:
% J. Veraart, E. Fieremans, DS Novikov
% Diffusion MRI noise mapping using random matrix theory
% MRM 76 (2016) 1582-1593
% Orignal script was modified by
% Ivan I. Maximov 22.12.2017 Ullevål/Oslo
% Ivan I. Maximov 22.11.2018, ver. 1.1

function denoise(fname)

addpath(genpath('/cluster/projects/p33/users/ivanma/MATLAB/matlab'));
f = nifti(fname);
img = f.dat(:,:,:,:);

signal = MPdenoising(img);

nfname = [fname(1:(end-4)) '_n.nii'];
N = nifti;
N = f;
N.dat.fname = nfname;
N.dat.dtype = 'FLOAT64-LE';
create(N);

N.dat(:,:,:,:) = signal;
end
% END