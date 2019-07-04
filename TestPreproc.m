% Demo runs of the RunPreproc function. Four different test-cases.
%
% OBS: To run, replace the hardcoded paths with your own data
%_______________________________________________________________________
%  Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

clear;

%----------------------
% Pick a test-case
% 1. MRI single-channel
% 2. MRI multi-channel
% 3. MRI multi-channel w. labels
% 4. CT
TESTCASE = 1;
%----------------------

if     TESTCASE == 1
    % MRI single-channel
    Nii = nifti(fullfile('example-data','mri-sc','0005-00001-000001-01.nii'));
elseif TESTCASE == 2        
    % MRI multi-channel
    Nii = nifti(char({fullfile('example-data','mri-mc','1282601181791321211150584332_T1.nii'),...
                      fullfile('example-data','mri-mc','1282601181791321211150584332_T2.nii'),...
                      fullfile('example-data','mri-mc','1282601181791321211150584332_Flair.nii')}));
elseif TESTCASE == 3
    % MRI multi-channel w. labels
    Nii       = cell(1,2);
    Nii{1}    = nifti;
    Nii{1}(1) = nifti(fullfile('example-data','labels','1-T1.nii'));
    Nii{1}(2) = nifti(fullfile('example-data','labels','1-FLAIR.nii'));    
    Nii{1}(3) = nifti(fullfile('example-data','labels','1-IR.nii'));    
    Nii{2}    = nifti; % Map label image to intensity image it was annotated on
    Nii{2}(1) = nifti;
    Nii{2}(2) = nifti(fullfile('example-data','labels','1-segm.nii'));
    Nii{2}(3) = nifti;
elseif TESTCASE == 4
    % CT
    Nii = nifti(fullfile('example-data','ct','3_s40271750-0004-00003-000001.nii'));
else
    error('Undefined test-case!')
end

%----------------------
% Set preprocessing options
%----------------------

opt.do.vx       = false;
opt.do.res_orig = false;
opt.do.real_mni = true;
opt.do.coreg    = true;
opt.do.denoise  = true;
opt.do.crop     = true;
opt.do.write2d  = true;
opt.pth_mtv     = fullfile('/home','mbrud','dev','mbrud','code','matlab','MTV-preproc');
if TESTCASE == 2 || TESTCASE == 3
    opt.do.reslice = true;
    opt.do.vx      = true;
end
if TESTCASE == 4
    opt.do.res_orig = true;
end

%----------------------
% Do preprocessing
%----------------------

[P,~,M] = RunPreproc(Nii,opt);

%----------------------
% Go back to native space orientation
%----------------------

C = numel(P{1});
for i=1:2
    for c=1:C
        if (i == 2 && numel(P{i}) == 1) || isempty(P{i}{c}), continue; end
        
        Mc = spm_get_space(P{i}{c}); 
        spm_get_space(P{i}{c},M{c}*Mc); 
    end
end