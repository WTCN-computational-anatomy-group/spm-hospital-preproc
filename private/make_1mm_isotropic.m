function nf = make_1mm_isotropic(f)

vx  = 1;
deg = 0;
Nii = nifti(f);

samp = sqrt(sum(Nii.mat(1:3,1:3).^2));
samp = samp./vx;

[img,dm,mat] = resample_img(Nii,samp,deg);

% New Nii
f             = Nii.dat.fname;
[nf,nam,ext] = fileparts(f);
nf           = fullfile(nf,['vx' nam ext]);

oNii         = nifti;
oNii.dat     = file_array(nf,dm,Nii.dat.dtype,Nii.dat.offset,Nii.dat.scl_slope,Nii.dat.scl_inter);
oNii.mat     = mat;
oNii.mat0    = mat;
oNii.descrip = 'Resampled';
create(oNii);
oNii.dat(:)  = img(:);
clear img    

delete(f);
%==========================================================================

%==========================================================================
function [img,dm,mat] = resample_img(Nii,samp,deg,bc)
% Resample an image using deg interpolation, with bc boundary conditions.
% If samp < 1, does down-sampling; if samp > 1, does up-sampling.
%_______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

if nargin < 3, deg = 0; end
if nargin < 4, bc  = 0; end

if numel(samp) == 1, samp = samp*ones([1 3]); end
if numel(deg)  == 1, deg  = deg*ones([1 3]);  end
if numel(bc)   == 1, bc   = bc*ones([1 3]);   end

samp(samp == 0) = 1;

% Input image properties
img  = Nii.dat(:,:,:);
mat0 = Nii.mat;
dm0  = size(img);

% Output image properties
vx            = sqrt(sum(mat0(1:3,1:3).^2));
% samp(vx >= 1) = 1;
D             = diag([samp 1]);
mat           = mat0/D;
dm            = floor(D(1:3,1:3)*dm0')';

% Make interpolation grid
[x0,y0,z0] = ndgrid(1:dm(1),1:dm(2),1:dm(3));

T = mat0\mat;    

x1 = T(1,1)*x0 + T(1,2)*y0 + T(1,3)*z0 + T(1,4);
y1 = T(2,1)*x0 + T(2,2)*y0 + T(2,3)*z0 + T(2,4);
z1 = T(3,1)*x0 + T(3,2)*y0 + T(3,3)*z0 + T(3,4);

% Resample
img                 = spm_bsplins(img,x1,y1,z1,[deg bc]);    
img(~isfinite(img)) = 0;
%==========================================================================