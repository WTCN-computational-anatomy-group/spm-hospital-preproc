% Demo runs of the RunPreproc function
%
% OBS: 
%   1. Add your own paths on line 24 (single-channel example) or 27 (multi-
%      channel example) or 33 (CT)
%   2. Add path to the denoising code on line 47, downloadable from:
%      '/home/mbrud/dev/mbrud/code/matlab/MTV-preproc'
%_______________________________________________________________________
%  Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

clear;

%----------------------
% Read image data
%----------------------

singlechannel = true; % Single- or multi-channel data?
mri           = true; % MRI or CT?

if mri
    % Data is MRI
    if singlechannel
        % One channel example
        P = '0005-00001-000001-01.nii';
    else
        % More than one channel example
        P = char({'1282601181791321211150584332_T1.nii',...
                  '1282601181791321211150584332_T2.nii',...
                  '1282601181791321211150584332_Flair.nii'});
    end   
else
    % Data is CT
    P = '3_s40271750-0004-00003-000001.nii';
end

Nii = nifti(P);

%----------------------
% Set preprocessing options
%----------------------

opt.do.real_mni = true;
opt.do.coreg    = true;
opt.do.denoise  = true;
opt.do.crop     = true;
opt.do.vx       = false;
opt.pth_mtv     = '/home/mbrud/dev/mbrud/code/matlab/MTV-preproc';
if ~singlechannel && mri
    opt.do.reslice = true;
    opt.do.vx      = true;
end
if ~mri
    opt.do.res_orig = true;
end

%----------------------
% Do preprocessing
%----------------------

[P,M] = RunPreproc(Nii,opt);

%----------------------
% Go back to native space orientation
%----------------------

C = numel(P);
for c=1:C
    Mc = spm_get_space(P{c}); 
    spm_get_space(P{c},M{c}*Mc); 
end