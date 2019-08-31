function img = brute_mni_align(img,tpm,flips,init0)
% MNI alignement for badly broken nifti files.
%
% FORMAT img = brute_mni_align(img,[tpm],[flips],[init0])
% img   - Image to fix (filename or spm_vol object)
% tmp   - Reference tissue probability map (filename or spm_vol object) [spm]
% flips - Try to flip dimensions [true]
% init0 - Initialise voxel-to-world with identity [false]
%
% Rigid align input image to MNI and overwrite voxel-to-world matrix in
% file. All possible permutations and inversion of dimensions are tried. 
% The one with highest mutual information is selected. Also compares 
% center-of-mass initialisation and trusting the inital header.
%
% /!\ Warning: this function may not conserve the left-right orientation.

    fname = '';
    if ischar(img)
        fname = img;
        img   = spm_vol(fname);
    end
    if nargin < 2 || isempty(tpm)
        tpm = spm_load_priors8(fullfile(spm('dir'),'tpm','TPM.nii'));
    end
    if ischar(tpm)
        tpm   = spm_vol(tpm);
    end
    if nargin < 3
        flips = true;
    end
    if flips
        fliplist = [0 1];
    else
        fliplist = 0;
    end
    if nargin < 4
        init0 = false;
    end
    
    % Hard-coded parameters
    samp0   = 8;            % Distance between sampled points in mm (initial)
    samp    = 4;            % Distance between sampled points in mm (final)
    fwhm    = 1;            % FWHM for fudge factor (final)
    fwhm0   = (0+1)*16;     % FWHM for fudge factor (initial)
    fwhm1   = (fwhm+1)*16;  % FWHM for fudge factor (coarse)
    affreg  = 'rigid';      % Regularisation
    
    % Build all possible combinations of flips and permutations
    prm = {};
    for flipx=fliplist
    for flipy=fliplist
    for flipz=fliplist
    for com=[0 1]
    for perms={[1 2 3], [1 3 2], [2 3 1], [2 1 3], [3 1 2], [3 2 1]}
        perms = perms{1};
        prm{end+1} = [perms com flipx flipy flipz];
    end
    end
    end
    end
    end
    
    % Try all combinations in parallel
    fprintf('Align: search orientations\n')
    if init0
        M0 = eye(4);
    else
        M0  = img.mat;
    end
    A  = cell(1,numel(prm));
    M  = cell(1,numel(prm));
    ll = zeros(1,numel(prm));
    parfor i=1:numel(prm)
        perms = prm{i}(1:3);
        com   = prm{i}(4);
        flipx = prm{i}(5);
        flipy = prm{i}(6);
        flipz = prm{i}(7);
        
        fprintf('. Permute [%d %d %d].', perms(1), perms(2), perms(3));
        if flips
            fprintf(' Flip [%d %d %d].', flipx, flipy, flipz);
        end
        fprintf(' Centre [%d].', com);
        fprintf('\n');
        
        img1 = img;
        M1   = M0;
        if com
            c = (img1.dim+1)/2;
            M1(1:3,4) = -M1(1:3,1:3)*c(:);
        end
        if flipx, M1(1,:) = -M1(1,:); end
        if flipy, M1(2,:) = -M1(2,:); end
        if flipz, M1(3,:) = -M1(3,:); end
        M1 = M1([perms 4],:);
        
        img1.mat  = M1;
        M{i} = M1;
        [A{i},ll(i)] = spm_maff8(img1, samp0, fwhm0, tpm, [], affreg);
        
    end
    
    % Pick the result with the best fit
    [~,imax] = max(ll);
    A = A{imax};
    M = M{imax};
    img.mat = M;
    prm = prm{imax};
    fprintf('Best parameters: Permute [%d %d %d].', prm(1), prm(2), prm(3));
    if flips
        fprintf(' Flip [%d %d %d].', prm(5), prm(6), prm(7));
    end
    fprintf(' Centre [%d].', prm(4));
    fprintf('\n');

    % Coarse affine registration (initialised with previous guess)
    fprintf('Align: coarse\n')
    A = spm_maff8(img, samp, fwhm1, tpm, A, affreg);
    % Final affine registration (initialised with previous guess)
    fprintf('Align: fine\n')
    A = spm_maff8(img, samp, fwhm, tpm, A, affreg);
    
    % Overwrite header
    img.mat = A*M;
    spm_create_vol(img);
    
    % Check reg
    spm_check_registration(img.fname, [tpm.V(1).fname, ',1']);
    
    % Return updated structure (or filename)
    if ~isempty(fname)
        img = fname;
    end
end