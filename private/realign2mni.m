function [Nii,M] = realign2mni(Nii,M,opt)
fprintf('Realigning to MNI...')

do_rigid = opt.rigid;
ix_realign = opt.ix_realign;
n_classes = opt.n_classes;

N = numel(Nii{1});
if nargin < 2
    M    = cell(1,N);
    M(:) = {eye(4)};
end
R = cell(1,N);
if ix_realign == 0
    for n=1:N
        f         = Nii{1}(n).dat.fname;
        R{n}      = rigid_align(f,do_rigid,n_classes);
        M{n}      = M{n}*R{n};
        Nii{1}(n) = nifti(f);
    end
else
    n         = ix_realign;
    f         = Nii{1}(n).dat.fname;
    R{n}      = rigid_align(f,do_rigid,n_classes);
    M{n}      = M{n}*R{n};
    Nii{1}(n) = nifti(f);
    for n1=1:N
        if n1 == n, continue; end
        f     = Nii{1}(n1).dat.fname;
        R{n1} = R{n};
        M{n1} = M{n1}*R{n1};
        spm_get_space(f, R{n}\Nii{1}(n1).mat);
        Nii{1}(n1) = nifti(f);
    end
end

if numel(Nii) > 1
    % Keep labels in alignment
    for n=1:N
        if n > numel(Nii{2}) || isempty(Nii{2}(n).dat), continue; end
        
        f         = Nii{2}(n).dat.fname;
        mat0      = Nii{2}(n).mat;
        spm_get_space(f,R{n}\mat0); 
        Nii{2}(n) = nifti(f); 
    end    
end
fprintf('done!\n')
%==========================================================================

%==========================================================================
function M = rigid_align(P,do_rigid,n_classes)
% Reposition an image by affine aligning to MNI space and Procrustes adjustment
% FORMAT rigid_align(P)
% P - name of NIfTI image
% M - Affine matrix
%
% OBS: Image will have the matrix in its header adjusted.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

if nargin < 2, do_rigid = true; end
if nargin < 3, n_classes = 6; end

% Load tissue probability data
tpm = fullfile(spm('dir'),'tpm','TPM.nii,');
tpm = [repmat(tpm,[n_classes 1]) num2str((1:n_classes)')];
tpm = spm_load_priors8(tpm);

% Do the affine registration
V = spm_vol(P);

M               = V(1).mat;
c               = (V(1).dim+1)/2;
V(1).mat(1:3,4) = -M(1:3,1:3)*c(:);
[Affine1,ll1]   = spm_maff8(V(1),8,(0+1)*16,tpm,[],'mni'); % Closer to rigid
Affine1         = Affine1*(V(1).mat/M);

% Run using the origin from the header
V(1).mat      = M;
[Affine2,ll2] = spm_maff8(V(1),8,(0+1)*16,tpm,[],'mni'); % Closer to rigid

% Pick the result with the best fit
if ll1>ll2, Affine  = Affine1; else Affine  = Affine2; end

Affine = spm_maff8(P,2,32,tpm,Affine,'mni'); % Heavily regularised
Affine = spm_maff8(P,2,1 ,tpm,Affine,'mni'); % Lightly regularised

% Load header
Nii = nifti(P);

% Generate mm coordinates of where deformations map from
x      = affind(rgrid(size(tpm.dat{1})),tpm.M);

% Generate mm coordinates of where deformation maps to
y1     = affind(x,inv(Affine));

% Weight the transform via GM+WM
weight = single(exp(tpm.dat{1})+exp(tpm.dat{2}));

% Weighted Procrustes analysis
[Affine,R]  = spm_get_closest_affine(x,y1,weight);
    
if do_rigid
    M = R;
else
    M = Affine;
end

% Invert
% R      = inv(R);

% Write the new matrix to the header
Nii.mat = M\Nii.mat;
create(Nii);
%==========================================================================

%==========================================================================
function x = rgrid(d)
x = zeros([d(1:3) 3],'single');
[x1,x2] = ndgrid(single(1:d(1)),single(1:d(2)));
for i=1:d(3)
    x(:,:,i,1) = x1;
    x(:,:,i,2) = x2;
    x(:,:,i,3) = single(i);
end
%==========================================================================

%==========================================================================
function y1 = affind(y0,M)
y1 = zeros(size(y0),'single');
for d=1:3
    y1(:,:,:,d) = y0(:,:,:,1)*M(d,1) + y0(:,:,:,2)*M(d,2) + y0(:,:,:,3)*M(d,3) + M(d,4);
end
%==========================================================================