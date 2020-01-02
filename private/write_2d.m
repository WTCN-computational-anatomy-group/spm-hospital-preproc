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
        
        nf = nm_reorient(nf,[],'ro',0);
        nf = make_1mm_isotropic(nf);
        nf = do_write2d(nf,deg,ax,sliceix);                
        
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
            
            nf = nm_reorient(nf,[],'ro',0);
            nf = make_1mm_isotropic(nf);
            nf = do_write2d(nf,deg,ax,sliceix);     
            
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
            
            nf = nm_reorient(nf,[],'ro',0);
            nf = make_1mm_isotropic(nf);
            nf = do_write2d(nf,deg,ax,sliceix);        

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
        sliceix = round(dm(3)/2);
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