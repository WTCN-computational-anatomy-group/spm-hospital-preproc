function Nii = crop(Nii,opt)

keep_neck = opt.keep_neck;

fprintf('Cropping...')
N        = numel(Nii{1});
prefix   = 'cr';
for n=1:N
    f = Nii{1}(n).dat.fname;
    
    atlas_crop(f,prefix,keep_neck);
    
    [pth,nam,ext] = fileparts(f);
    nf            = fullfile(pth,[prefix nam ext]);
    
    delete(f);
    Nii{1}(n) = nifti(nf);
end
fprintf('done!\n')
%==========================================================================

%==========================================================================
function [R,bb] = atlas_crop(P,prefix,keep_neck)
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
if nargin<2, prefix    = 'cr'; end
if nargin<3, keep_neck = true; end

pth_mu = fullfile(spm('dir'),'tpm','TPM.nii,');

% Image params
Vf = spm_vol(P);
Mf = Vf(1).mat;
df = Vf(1).dim;
vf = sqrt(sum(Mf(1:3,1:3).^2));

% Atlas params
Vmu = spm_vol(pth_mu);
Mmu = Vmu(1).mat;
dmu = Vmu(1).dim;
vmu = sqrt(sum(Mmu(1:3,1:3).^2));

% Register atlas to image to get get R (so that Mmu\R*Mf)
mu               = spm_load_priors8(Vmu);    
c                = (Vf(1).dim+1)/2;
Vf(1).mat(1:3,4) = -Mf(1:3,1:3)*c(:);
[Affine1,ll1]    = spm_maff8(Vf(1),8,(0+1)*16,mu,[],'mni'); % Closer to rigid
Affine1          = Affine1*(Vf(1).mat/Mf);

% Run using the origin from the header
Vf(1).mat     = Mf;
[Affine2,ll2] = spm_maff8(Vf(1),8,(0+1)*16,mu,[],'mni'); % Closer to rigid

% Pick the result with the best fit
if ll1>ll2, R  = Affine1; else R = Affine2; end

% Fit final
R = spm_maff8(P,4,32,mu,R,'mni');
R = spm_maff8(P,4,1,mu,R,'mni');
 
% Find bounding box of the projected template
M = Mf\(R\Mmu);
d = dmu;
c = [1    1    1    1
     1    1    d(3) 1
     1    d(2) 1    1
     1    d(2) d(3) 1
     d(1) 1    1    1
     d(1) 1    d(3) 1
     d(1) d(2) 1    1
     d(1) d(2) d(3) 1]';  
c  = M(1:3,1:4)*c;
mn = min(c,[],2)';
mx = max(c,[],2)';
bb = [mn; mx];

if keep_neck
    % Find Dorso-ventral dimension in image space
    c        = [1 1 1    1
                1 1 d(3) 1]';
    c        = M(1:3,1:4)*c;
    [~,dv]   = max(abs(c(:,1) - c(:,2)) .* vf');
    bb(:,dv) = [1; df(dv)];
end

% Do cropping
subvol(Vf,bb,prefix);  
%==========================================================================