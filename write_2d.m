function [Nii,P] = write_2d(Nii,dir_out2d,deg,axis_2d,sliceix)
if nargin < 3, deg      = 0;    end
if nargin < 4, axis_2d  = 3;    end
if nargin < 5, sliceix    = [];   end

fprintf('Writing 2D...')
N = numel(Nii{1});
for n=1:N
    f           = Nii{1}(n).dat.fname;  
    [~,nam,ext] = fileparts(f);
    nf          = fullfile(dir_out2d,['2d' nam ext]);  
    copyfile(f,nf);
    nf          = nm_reorient(nf);
    
    nf        = do_write2d(nf,deg,axis_2d,sliceix);        
    Nii{1}(n) = nifti(nf);
end

if numel(Nii) > 1
    % Labels too
    for n=1:N
        if isempty(Nii{2}(n).dat), continue; end
        
        f           = Nii{2}(n).dat.fname;  
        [~,nam,ext] = fileparts(f);
        nf          = fullfile(dir_out2d,['2d' nam ext]);  
        copyfile(f,nf);
        nf          = nm_reorient(nf);

        nf        = do_write2d(nf,deg,axis_2d,sliceix);        
        Nii{2}(n) = nifti(nf);
    end    
end

P    = cell(1,2);
P{1} = cell(1,N);
P{2} = cell(1,N);
for i=1:2
    for n=1:N
        if (i == 2 && numel(Nii) == 1) || isempty(Nii{i}(n).dat), continue; end
        
        P{i}{n} = Nii{i}(n).dat.fname;        
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
        sliceix = floor(dm(1)/2) + 1;
    end
    bb = [sliceix sliceix;-inf inf;-inf inf];   
elseif axis_2d == 2
    if isempty(sliceix) 
        sliceix = floor(dm(2)/2) + 1;
    end
    bb = [-inf inf;sliceix sliceix;-inf inf];
elseif axis_2d == 3
    if isempty(sliceix) 
        sliceix = floor(dm(3)/2) + 1;
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