function Nii = crop(Nii)

fprintf('Cropping...')
N        = numel(Nii{1});
prefix   = 'cr';
for n=1:N
    f = Nii{1}(n).dat.fname;
    
    atlas_crop(f,prefix);
    
    [pth,nam,ext] = fileparts(f);
    nf            = fullfile(pth,[prefix nam ext]);
    
    delete(f);
    Nii{1}(n) = nifti(nf);
end
fprintf('done!\n')
%==========================================================================

%==========================================================================
function [Affine,bb] = atlas_crop(P,prefix)
% Crop image to SPM atlas size
% FORMAT [Affine,bb] = atlas_crop(P,Affine,prefix,rem_neck)
% P        - Path to NIfTI file
% prefix   - File prefix (if empty -> overwrites) ['']
% rem_neck - Remove neck/spine [false]
% bb - Computed bounding box
%
% This function rigidly registers the SPM atlas to an image and then
% removes image data outside of the atlas.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging
if nargin<2, prefix   = 'cr'; end

% Locate TPM.nii in SPM
pthmu = fullfile(spm('dir'),'tpm','TPM.nii,');
Vf    = spm_vol(P);
Vmu   = spm_vol(pthmu);
matf  = Vf(1).mat;
matmu = Vmu(1).mat;
dmf   = Vf(1).dim;
dmmu  = Vmu(1).dim;

mu = spm_load_priors8(Vmu);    

c                = (Vf(1).dim+1)/2;
Vf(1).mat(1:3,4) = -matf(1:3,1:3)*c(:);
[Affine1,ll1]    = spm_maff8(Vf(1),8,(0+1)*16,mu,[],'mni'); % Closer to rigid
Affine1          = Affine1*(Vf(1).mat/matf);

% Run using the origin from the header
Vf(1).mat      = matf;
[Affine2,ll2] = spm_maff8(Vf(1),8,(0+1)*16,mu,[],'mni'); % Closer to rigid

% Pick the result with the best fit
if ll1>ll2, Affine  = Affine1; else Affine  = Affine2; end

Affine = spm_maff8(P,4,32,mu,Affine,'mni');
Affine = spm_maff8(P,4,1,mu,Affine,'mni');

% Voxel locations in input
T  = matf\(Affine\matmu);

% Corners
crnmu = [1    1    1    1
        1    1    dmmu(3) 1
        1    dmmu(2) 1    1
        1    dmmu(2) dmmu(3) 1
        dmmu(1) 1    1    1
        dmmu(1) 1    dmmu(3) 1
        dmmu(1) dmmu(2) 1    1
        dmmu(1) dmmu(2) dmmu(3) 1]';  

crnf = T*crnmu;

% Bounding-box
bb = zeros(2,3);
for i=1:3
    X       = crnf(i,:);
    bb(1,i) = max(X);
    bb(2,i) = min(X);
end

% Do cropping
subvol(Vf,bb,prefix);      
%==========================================================================