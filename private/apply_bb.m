function Nii = apply_bb(Nii, vx, dim_out)

fprintf('Applying bounding-box...')
N        = numel(Nii{1});
prefix   = 'bb';
for n=1:N
    f  = Nii{1}(n).dat.fname;    
    nf = do_apply_bb(f,prefix,vx,dim_out);
    
    delete(f);
    Nii{1}(n) = nifti(nf);
end
fprintf('done!\n')
%==========================================================================

%==========================================================================
function nf = do_apply_bb(f,prefix,vx_in,dim_bb)
if numel(dim_bb) == 1
    dim_bb = dim_bb*ones(1,3);
end
% get input
Nii_in = nifti(f);
mat_in = Nii_in.mat;
dim_in = Nii_in.dat.dim;
if isempty(vx_in), vx_in  = sqrt(sum(mat_in(1:3,1:3).^2));
else,              vx_in = vx_in(1)*ones(1,3);
end
% get tpm
V_tpm   = spm_vol(fullfile(spm('dir'),'tpm','TPM.nii,'));
mat_tpm = V_tpm(1).mat;
dim_tpm = V_tpm(1).dim;
vx_tpm  = sqrt(sum(mat_tpm(1:3,1:3).^2));
% get output
mat_D   = diag([vx_tpm./vx_in 1]);
mat_out = mat_tpm/mat_D;
dim_out = floor(mat_D(1:3,1:3)*dim_tpm')';
% adjust image dimensions
if all(isfinite(dim_bb))
    mat_bb  = spm_matrix(-round(dim_bb - dim_out)/2);
    mat_out = mat_out*mat_bb;
    dim_out = dim_bb;
end
% make interpolation grid
[x0,y0,z0] = ndgrid(1:dim_out(1),1:dim_out(2),1:dim_out(3));
T  = mat_in\mat_out;    
x1 = T(1,1)*x0 + T(1,2)*y0 + T(1,3)*z0 + T(1,4);
y1 = T(2,1)*x0 + T(2,2)*y0 + T(2,3)*z0 + T(2,4);
z1 = T(3,1)*x0 + T(3,2)*y0 + T(3,3)*z0 + T(3,4);
% image data
img = Nii_in.dat();
mn  = min(img(:));
mx  = max(img(:));
% resample
img                 = spm_bsplinc(img, [1 1 1  0 0 0]);
img                 = spm_bsplins(img,x1,y1,z1, [1 1 1  0 0 0]);    
img(~isfinite(img)) = 0;
img                 = min(mx, max(mn, img));
% make output nifti
[pth,nam,ext]   = fileparts(f);
nf              = fullfile(pth,[prefix nam ext]);
Nii_out         = nifti;
Nii_out.dat     = file_array(nf,dim_out,Nii_in.dat.dtype,Nii_in.dat.offset,Nii_in.dat.scl_slope,Nii_in.dat.scl_inter);
Nii_out.mat     = mat_out;
Nii_out.mat0    = mat_out;
Nii_out.descrip = 'bb';
create(Nii_out);
Nii_out.dat(:) = img(:);
%==========================================================================