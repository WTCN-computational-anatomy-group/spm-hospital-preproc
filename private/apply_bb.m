function Nii = apply_bb(Nii)

fprintf('Applying bounding-box...')
N        = numel(Nii{1});
prefix   = 'bb';
for n=1:N
    f = Nii{1}(n).dat.fname;
    
    nf = do_apply_bb(f,prefix);
    
    delete(f);
    Nii{1}(n) = nifti(nf);
end
fprintf('done!\n')
%==========================================================================

%==========================================================================
function nf = do_apply_bb(f,prefix)
Nii0 = nifti(f);
mat0 = Nii0.mat;
dim0 = Nii0.dat.dim;
vx0  = sqrt(sum(mat0(1:3,1:3).^2));

V_tpm = spm_vol(fullfile(spm('dir'),'tpm','TPM.nii,'));
M_tpm = V_tpm(1).mat;
d_tpm = V_tpm(1).dim;
vx_tpm  = sqrt(sum(M_tpm(1:3,1:3).^2));

% Output image properties
D   = diag([vx_tpm./vx0 1]);
mat = M_tpm/D;
dm  = floor(D(1:3,1:3)*d_tpm')';

% Make interpolation grid
[x0,y0,z0] = ndgrid(1:dm(1),1:dm(2),1:dm(3));
T  = mat0\mat;    
x1 = T(1,1)*x0 + T(1,2)*y0 + T(1,3)*z0 + T(1,4);
y1 = T(2,1)*x0 + T(2,2)*y0 + T(2,3)*z0 + T(2,4);
z1 = T(3,1)*x0 + T(3,2)*y0 + T(3,3)*z0 + T(3,4);

% Get image
img = Nii0.dat();
mn  = min(img(:));
mx  = max(img(:));
    
% Resample
img                 = spm_bsplinc(img, [1 1 1  0 0 0]);
img                 = spm_bsplins(img,x1,y1,z1, [1 1 1  0 0 0]);    
img(~isfinite(img)) = 0;
img                 = min(mx, max(mn, img));

% New Nii
[pth,nam,ext] = fileparts(f);
nf            = fullfile(pth,[prefix nam ext]);
oNii          = nifti;
oNii.dat      = file_array(nf,dm,Nii0.dat.dtype,Nii0.dat.offset,Nii0.dat.scl_slope,Nii0.dat.scl_inter);
oNii.mat      = mat;
oNii.mat0     = mat;
oNii.descrip  = 'bb';
create(oNii);
oNii.dat(:)   = img(:);
%==========================================================================