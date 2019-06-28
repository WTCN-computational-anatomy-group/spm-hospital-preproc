function [P,M] = RunPreproc(Nii,opt)
% Some basic preprocessing of hospital neuroimaging data.
%
% Nii - 1xC nifti struct, or empty to use file selector, of patient images
% opt - Preprocessing options
% P   - Path of preprocessed image(s)
% M   - Orientation matrices to go back to native space orientation as:
%           Mc = spm_get_space(P{c}); 
%           spm_get_space(f,M{c}*Mc); 
%_______________________________________________________________________
%  Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

if nargin == 0
   Nii = nifti(spm_select(Inf,'nifti','Please select images of one patient'));   
end

% Set default options
if nargin < 2, opt = struct; end
% Do
if ~isfield(opt,'do'),          opt.do          = struct; end
if ~isfield(opt.do,'res_orig'), opt.do.res_orig = false; end
if ~isfield(opt.do,'real_mni'), opt.do.real_mni = true; end
if ~isfield(opt.do,'crop'),     opt.do.crop     = true; end
if ~isfield(opt.do,'coreg'),    opt.do.coreg    = true; end
if ~isfield(opt.do,'denoise'),  opt.do.denoise  = true; end
if ~isfield(opt.do,'reslice'),  opt.do.reslice  = true; end
if ~isfield(opt.do,'vx'),       opt.do.vx       = true; end
% Output directory
if ~isfield(opt,'dir_out'),     opt.dir_out     = 'output'; end
% Reslice options
if ~isfield(opt,'reslice'),     opt.reslice     = struct; end
if ~isfield(opt.reslice,'ref'), opt.reslice.ref = 1; end
% Voxel size options
if ~isfield(opt,'vx'),          opt.vx          = struct; end
if ~isfield(opt.vx,'size'),     opt.vx.size     = 1; end
if ~isfield(opt.vx,'def'),      opt.vx.deg      = 4; end
% Crop options
if ~isfield(opt,'crop'),        opt.crop        = struct; end
if ~isfield(opt.crop,'neck'),   opt.crop.neck   = false; end
% Path to denoising toolbox (https://github.com/WCHN/mtv-preproc)
if ~isfield(opt,'pth_mtv'),     opt.pth_mtv     = '/home/mbrud/dev/mbrud/code/matlab/MTV-preproc'; end

% Add denoising toolbox to path
addpath(opt.pth_mtv)

% Because it is possible to include labels (will add more abt this soon)
if ~iscell(Nii)
    Nii = {Nii};
end
C = numel(Nii{1});

if ~(exist(opt.dir_out,'dir') == 7)  
    % Create output directory
    mkdir(opt.dir_out);  
end
% Make sure output directory is encoded by its full path
s           = what(opt.dir_out);
opt.dir_out = s.path;

% Copy so to not overwrite originals
Nii = make_copies(Nii,opt.dir_out);

% Initialise orientation matrix (to go back to native space)
M    = cell(1,C);
M(:) = {eye(4)};

if opt.do.res_orig
    % Reset origin (important for CT)
    [Nii,M] = reset_origin(Nii);
end

if opt.do.real_mni
    % Realing to MNI space
    [Nii,M] = realign2mni(Nii,M);
end

if opt.do.coreg
    % Coreg
    Nii = coreg(Nii);
end

if opt.do.denoise
    % Denoise
    Nii = denoise(Nii);

    % Coreg (one more time after denoising)
    Nii = coreg(Nii);
end

if opt.do.crop
    % Remove uneccesary data
    Nii = crop(Nii,opt.crop.neck);
end

if opt.do.vx
    % Set same voxel size
    Nii = resample_images(Nii,opt.vx.size,opt.vx.deg);
end

if opt.do.reslice
    % Make images same dimensions
    Nii = reslice(Nii,opt.reslice.ref);
end

% Give path of preprocessed image(s) as output and save orientation
% matrices
P = cell(1,C);
for c=1:C
    P{c}      = Nii{1}(c).dat.fname;
    [pth,nam] = fileparts(P{c});
    nP        = fullfile(pth,['mat' nam '.mat']);
    Mc        = M{c};
    save(nP,'Mc')
end
%==========================================================================