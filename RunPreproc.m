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

if opt.do.denoise || opt.do.superres
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

if numel(Nii) > 1
    % Collapse labels
    Nii = collapse_labels(Nii,opt.labels.part);    
end

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

if opt.do.crop
    % Remove uneccesary data
    Nii = crop(Nii,opt.crop);
end

if opt.do.coreg
    % Coreg
    Nii = coreg(Nii,opt.coreg);
end

if opt.do.denoise && ~opt.do.superres
    % Denoise
    Nii = denoise(Nii);

    % Coreg (one more time after denoising)
    Nii = coreg(Nii,opt.coreg);
end

% The below steps are for creating images of equal size, either by MTV
% super-resolution, or by just simply reslicing
if opt.do.superres
    % Super-resolve
    Nii = superres(Nii,opt.superres);
    
    % Coreg (one more time after super-resolving)
    Nii = coreg(Nii,opt.coreg);
else
    if opt.do.vx
        % Set same voxel size
        Nii = resample_images(Nii,opt.vx);
    end

    if opt.do.reslice
        % Make images same dimensions
        [Nii,M] = reslice_images(Nii,M,opt.reslice);
    end
end

if numel(Nii) > 1
    % Reslice labels
    Nii = reslice_labels(Nii,opt.reslice);
end

Seg_pths = {};
if opt.do.segment
    % Run SPM12 segmentation
    Seg_pths = segment_preproc8(Nii,opt.segment);
end

P2d = {};
if opt.do.write2d
    % Write 2D versions
    [~,P2d] = write_2d(Nii,Seg_pths,opt.dir_out2d,opt.write2d);
end

% Allocate output
C             = numel(Nii{1});
out           = struct;
out.pth.im    = cell(1,C);
out.pth.im2d  = cell(1,C);
out.pth.lab   = cell(1,C);
out.pth.lab2d = cell(1,C);
out.pth.seg   = {}; 
out.mat       = cell(1,C);
for i=1:2
    for c=1:C
        if (i == 2 && numel(Nii) == 1) || isempty(Nii{i}(c).dat), continue; end
        
        if i == 1
            out.pth.im{c}  = Nii{i}(c).dat.fname;
        else
            out.pth.lab{c} = Nii{i}(c).dat.fname;
        end
                            
        if ~isempty(Seg_pths)
            out.pth.seg = Seg_pths;
        end
                
        if ~isempty(P2d)
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