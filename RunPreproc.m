function out = RunPreproc(Nii,opt)
% Some basic preprocessing of hospital neuroimaging data.
%
% INPUT
% Nii - 1xC nifti struct (or empty to use file selector) of patient images
% opt - Preprocessing options
%
% OUTPUT
% out.pth.im    - Cell array of paths to preprocessed image(s)
% out.pth.im2d  - Cell array of paths to 2D versions of preprocessed
%                 image(s) (if opt.do.write2d = true)
% out.pth.lab   - Cell array of path to label image (if labels given)
% out.pth.lab2d - Cell array of path to 2D version of label image (if 
%                 labels given and if opt.do.write2d = tru)
% out.mat       - Orientation matrices to go back to native space 
%                 orientation as:
%                   Mc = spm_get_space(P{c}); 
%                   spm_get_space(f,M{c}*Mc); 
%_______________________________________________________________________
%  Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

if nargin == 0
   Nii = nifti(spm_select(Inf,'nifti','Please select patient image(s)'));   
end

% Set options
if nargin < 2, opt = struct; end
opt = get_default_opt(opt);

if opt.do.denoise
    % Add denoising toolbox to path
    addpath(opt.pth_mtv)
end

% Because it is possible to include labels, in the second index of Nii 
% (i.e. Nii{2})
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

% Copy (so to not overwrite originals)
Nii = make_copies(Nii,opt.dir_out);

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
    Nii = coreg(Nii,opt.coreg.ref);
end

if opt.do.denoise
    % Denoise
    Nii = denoise(Nii);

    % Coreg (one more time after denoising)
    Nii = coreg(Nii,opt.coreg.ref);
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

if opt.do.write2d
    % Write 2D versions
    [~,P2d] = write_2d(Nii,opt.dir_out2d,opt.write2d.deg,opt.write2d.axis_2d,opt.write2d.sliceix);
end

% Allocate output
out           = struct;
out.pth.im    = cell(1,C);
out.pth.im2d  = cell(1,C);
out.pth.lab   = cell(1,C);
out.pth.lab2d = cell(1,C);
out.mat       = cell(1,C);
for i=1:2
    for c=1:C
        if (i == 2 && numel(Nii) == 1) || isempty(Nii{i}(c).dat), continue; end
        
        if i == 1
            out.pth.im{c}  = Nii{i}(c).dat.fname;
        else
            out.pth.lab{c} = Nii{i}(c).dat.fname;
        end
                            
        if exist('P2d','var')
            if i == 1
                out.pth.im2d{c}  = P2d{i}{c};
            else
                out.pth.lab2d{c} = P2d{i}{c};
            end
        end
                   
        if opt.do.writemat
            [pth,nam] = fileparts(P{i}{c});   
            nP        = fullfile(pth,['mat' nam '.mat']);
            Mc        = M{c};
            save(nP,'Mc')
        end

        if i == 1
            out.mat{c} = M{c};
        end  
    end
end
%==========================================================================