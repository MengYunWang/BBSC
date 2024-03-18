% CITE APPROPRIATE ARTICLE
% E. Kellner, B. Dhital, VG Kiselev, M. Reisert
% Gibbs-ringing artifact removal based on local subvoxel-shifts
% MRM 76 (2016) 1574-1581
% unring - tool for removal of the Gibbs ringing artefact
% Usage: outvol = unring(invol,params)
% Options: invol - input volume 
%          params - 3x1 array with [minW maxW nsh]
%                     nsh discretization of subpixel spaceing (default 20)
%                     minW  left border of window used for TV computation (default 1)
%                     maxW  right border of window used for TV computation (default 3)
% The original script was modified by
% Ivan I. Maximov 22.12.2017, Ullevål/Oslo
% Ivan I. Maximov 22.11.2018, ver. 1.1

function v = unring(v,params)

if nargin == 1
    params = [1 3 20];
end

addpath(genpath('/cluster/projects/p33/users/ivanma/MATLAB/matlab'));

unix('export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/cluster/projects/p33/users/ivanma/MATLAB/matlab/unring');

f = nifti(v);
img = f.dat(:,:,:,:);
data = ringRm(img,params);

nfname = [v(1:(end-5)) 'r.nii'];
N = nifti;
N = f;
N.dat.fname = nfname;
N.dat.dtype = 'FLOAT64-LE';
create(N);

N.dat(:,:,:,:) = data;

end
% END