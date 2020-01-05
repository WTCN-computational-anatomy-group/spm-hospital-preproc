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

Vf = spm_vol(P);
Vm = spm_vol(pthmu);
Mf = Vf(1).mat;
Mm = Vm(1).mat;
df = Vf(1).dim;
dm = Vm(1).dim;
vf = sqrt(sum(Mf(1:3,1:3).^2));
vm = sqrt(sum(Mm(1:3,1:3).^2));

tpm = spm_load_priors8(Vm);    

c                = (Vf(1).dim+1)/2;
Vf(1).mat(1:3,4) = -Mf(1:3,1:3)*c(:);
[Affine1,ll1]    = spm_maff8(Vf(1),8,(0+1)*16,tpm,[],'mni'); % Closer to rigid
Affine1          = Affine1*(Vf(1).mat/Mf);

% Run using the origin from the header
Vf(1).mat      = Mf;
[Affine2,ll2] = spm_maff8(Vf(1),8,(0+1)*16,tpm,[],'mni'); % Closer to rigid

% Pick the result with the best fit
if ll1>ll2, Affine  = Affine1; else Affine  = Affine2; end

Affine = spm_maff8(P,4,32,tpm,Affine,'mni');
Affine = spm_maff8(P,4,1,tpm,Affine,'mni');

% Template to world (inc affine transformation)
Tm = Affine\Mm;

% Keep neck
fc = [1    1    1    1
      1    1    df(3) 1
      1    df(2) 1    1
      1    df(2) df(3) 1
      df(1) 1    1    1
      df(1) 1    df(3) 1
      df(1) df(2) 1    1
      df(1) df(2) df(3) 1]';  

fc = Mf(1:3,1:4)*fc;

% bounding box
mn  = min(fc,[],2)';
mx  = max(fc,[],2)';
bbf = [mn; mx];

mmn = bbf(1,:);
mmx = bbf(2,:);
    
% Corners
tc = [1    1    1    1
      1    1    dm(3) 1
      1    dm(2) 1    1
      1    dm(2) dm(3) 1
      dm(1) 1    1    1
      dm(1) 1    dm(3) 1
      dm(1) dm(2) 1    1
      dm(1) dm(2) dm(3) 1]';  

tc = Tm(1:3,1:4)*tc;

% bounding box
mn  = min(tc,[],2)';
mx  = max(tc,[],2)';

mn(3) = mmn(3);
% mx(3) = mmx(3);

bbm = [mn; mx];

vmn = bbm(1,:);
vmx = bbm(2,:);
mn(isnan(mn)) = vmn(isnan(mn));
mx(isnan(mx)) = vmx(isnan(mx));

mat    = spm_matrix([mn 0 0 0 vf])*spm_matrix([-1 -1 -1]);
imgdim = ceil(mat \ [mx 1]' - 0.1)';

% Produce new image
VO            = Vf;
[pth,nam,ext] = fileparts(Vf.fname);
VO.fname      = fullfile(pth,[prefix nam ext]);
VO.dim(1:3)   = imgdim(1:3);
VO.mat        = mat;
VO = spm_create_vol(VO);
for i = 1:imgdim(3)
    M = inv(spm_matrix([0 0 -i])*inv(VO.mat)*Vf.mat);
    img = spm_slice_vol(Vf, M, imgdim(1:2), 1); % (linear interp)
    spm_write_plane(VO, img, i);
end
%==========================================================================