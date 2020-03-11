function Nii = resample_images(Nii,opt)

vx  = opt.size;
deg = opt.deg;

fprintf('Changing voxel-sizes...')
N = numel(Nii{1});

for n=1:N
    samp = sqrt(sum(Nii{1}(n).mat(1:3,1:3).^2));
    samp = samp./vx;
    
    [img,dm,mat] = resample_img(Nii{1}(n),samp,deg);
    
    % New Nii
    f             = Nii{1}(n).dat.fname;
    [pth,nam,ext] = fileparts(f);
    pth           = fullfile(pth,['vx' nam ext]);

    oNii         = nifti;
    oNii.dat     = file_array(pth,dm,Nii{1}(n).dat.dtype,Nii{1}(n).dat.offset,Nii{1}(n).dat.scl_slope,Nii{1}(n).dat.scl_inter);
    oNii.mat     = mat;
    oNii.mat0    = mat;
    oNii.descrip = 'Resampled';
    create(oNii);
    oNii.dat(:)  = img(:);
    clear img
    
    Nii{1}(n) = oNii;
    
    delete(f);
end

fprintf('done!\n')
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
mn   = min(img(:));
mx   = max(img(:));
    
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
img                 = spm_bsplinc(img, [deg bc]);
img                 = spm_bsplins(img,x1,y1,z1, [deg bc]);    
img(~isfinite(img)) = 0;
img                 = min(mx, max(mn, img));
%==========================================================================