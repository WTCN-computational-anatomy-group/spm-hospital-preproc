function P = write_2d(Nii,Seg_pths,dir_out2d,opt)

deg     = opt.deg;
sliceix = opt.sliceix;

nam_slice = {'sag','cor','ax'};

N    = numel(Nii{1});
P    = cell(1,3);
P{1} = cell(1,N);
P{2} = cell(1,N);
P{3} = cell(1,N);
    
fprintf('Writing 2D...')
for ax=1:3 % loop over anatomical axis   
    
    dir_out = fullfile(dir_out2d,['2D-' nam_slice{ax}]);    
    if ~(exist(dir_out,'dir') == 7)  
        % Create 2D output directory
        mkdir(dir_out);  
    end

    for n=1:N
        f           = Nii{1}(n).dat.fname;  
        [~,nam,ext] = fileparts(f);
        nf          = fullfile(dir_out,['2d-' nam_slice{ax} '-' nam ext]);  
        copyfile(f,nf);
        
        nf = make_1mm_isotropic(nf);
        nf = nm_reorient(nf,[],'ro',0);
        nf = atlas_crop(nf);
        nf = do_write2d(nf,deg,ax,sliceix);              
        
        f  = fullfile(dir_out,[nam ext]);
        movefile(nf,f);
        nf = f;
        
        P{1}{n} = nf;
    end

    if numel(Nii) > 1
        % Labels too
        for n=1:N
            if isempty(Nii{2}(n).dat), continue; end

            f           = Nii{2}(n).dat.fname;  
            [~,nam,ext] = fileparts(f);
            nf          = fullfile(dir_out,['2d-' nam_slice{ax} '-' nam ext]);  
            copyfile(f,nf);
                        
            nf = make_1mm_isotropic(nf);
            nf = nm_reorient(nf,[],'ro',0);
            nf = atlas_crop(nf);
            nf = do_write2d(nf,deg,ax,sliceix);     
            
            f  = fullfile(dir_out,[nam ext]);
            movefile(nf,f);
            nf = f;
        
            P{2}{n} = nf;
        end    
    end

    if ~isempty(Seg_pths)
        % Write 2D versions of 'c' segmentations

        Nii_seg = nifti;

        K = numel(Seg_pths{1});
        for k=1:K
            Nii_seg(k) = nifti(Seg_pths{1}{k});
        end

        for k=1:K
            f           = Nii_seg(k).dat.fname;  
            [~,nam,ext] = fileparts(f);
            nf          = fullfile(dir_out,['2d-' nam_slice{ax} '-' nam ext]);  
            copyfile(f,nf);
            
            nf = make_1mm_isotropic(nf);
            nf = nm_reorient(nf,[],'ro',0);
            nf = atlas_crop(nf);
            nf = do_write2d(nf,deg,ax,sliceix);      

            f  = fullfile(dir_out,[nam ext]);
            movefile(nf,f);
            nf = f;
            
            P{3}{k} = nf;
        end
    end
end

fprintf('done!\n')
%==========================================================================

%==========================================================================
function nfname = do_write2d(fname,deg,axis_2d,sliceix)
if nargin < 2, deg      = 0;  end
if nargin < 3, axis_2d  = 3;  end
if nargin < 4, sliceix  = []; end

% Create bounding box
V  = spm_vol(fname);
dm = V.dim;
if axis_2d     == 1
    if isempty(sliceix)        
        sliceix = round(dm(1)/2);
    end
    bb = [sliceix sliceix;-inf inf;-inf inf];   
elseif axis_2d == 2
    if isempty(sliceix) 
        sliceix = round(dm(2)/2);
    end
    bb = [-inf inf;sliceix sliceix;-inf inf];
elseif axis_2d == 3
    if isempty(sliceix) 
        sliceix = round(dm(3) - 95);
    end
    bb = [-inf inf;-inf inf;sliceix sliceix];
end                

% Crop according to bounding-box
subvol(V,bb','sv',deg);      

[pth,nam,ext] = fileparts(fname);
nfname        = fullfile(pth,['sv' nam ext]);  
delete(fname);
    
if axis_2d == 1 || axis_2d == 2
    % Make sure 1D plane is in z dimension
    Nii  = nifti(nfname);
    mat  = Nii.mat;
    
    % Permute image data and apply permutation matrix to orientation matrix
    if axis_2d == 1
        img = permute(Nii.dat(:,:,:),[2 3 1]);            
        P   = [0 1 0 0; 0 0 1 0; 1 0 0 0; 0 0 0 1];
    else
        img = permute(Nii.dat(:,:,:),[1 3 2]);        
        P   = [1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1];
    end   
    mat     = P*mat*P';
    dm      = [size(img) 1];
    
    % Overwrite image data
    VO             = spm_vol(nfname);
    VO.dim(1:3)    = dm(1:3);        
    VO.mat         = mat;
    VO             = spm_create_vol(VO);        
    Nii            = nifti(VO.fname);    
    Nii.dat(:,:,:) = img; 
end
%==========================================================================

%==========================================================================
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

%==========================================================================
function nP = atlas_crop(P,prefix)
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

% tpm = spm_load_priors8(Vm);    
% 
% c                = (Vf(1).dim+1)/2;
% Vf(1).mat(1:3,4) = -Mf(1:3,1:3)*c(:);
% [Affine1,ll1]    = spm_maff8(Vf(1),8,(0+1)*16,tpm,[],'mni'); % Closer to rigid
% Affine1          = Affine1*(Vf(1).mat/Mf);
% 
% % Run using the origin from the header
% Vf(1).mat      = Mf;
% [Affine2,ll2] = spm_maff8(Vf(1),8,(0+1)*16,tpm,[],'mni'); % Closer to rigid
% 
% % Pick the result with the best fit
% if ll1>ll2, Affine  = Affine1; else Affine  = Affine2; end
% 
% Affine = spm_maff8(P,4,32,tpm,Affine,'mni');
% Affine = spm_maff8(P,4,1,tpm,Affine,'mni');

% Template to world (inc affine transformation)
Tm = Mm;

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
nP            = fullfile(pth,[prefix nam ext]);
VO.fname      = nP;
VO.dim(1:3)   = imgdim(1:3);
VO.mat        = mat;
VO = spm_create_vol(VO);
for i = 1:imgdim(3)
    M = inv(spm_matrix([0 0 -i])*inv(VO.mat)*Vf.mat);
    img = spm_slice_vol(Vf, M, imgdim(1:2), 1); % (linear interp)
    spm_write_plane(VO, img, i);
end
delete(Vf.fname);
%==========================================================================