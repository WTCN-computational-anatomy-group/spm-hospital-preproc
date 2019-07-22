% Demo runs of the RunPreproc function. A couple of different test-cases
% are provided.
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
% 5. Hospital MRI superres
% 6. MRI single-channel with 2D version
% 7. MRI super-res with cell array input
TESTCASE = 7;
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
elseif TESTCASE == 5      
    % Hospital MRI superres
    Nii = nifti(char({fullfile('example-data','mri-ts','0002-00001-000021-01.nii'),...
                      fullfile('example-data','mri-ts','0003-00001-000001-01.nii'),...
                      fullfile('example-data','mri-ts','0004-00001-000001-01.nii')}));
elseif TESTCASE == 6
    % MRI single-channel with 2D version
    Nii = nifti(fullfile('example-data','mri-sc','0005-00001-000001-01.nii'));                  
elseif TESTCASE == 7
    % MRI super-res with cell array input
    ix  = [1 2 1 2 3];
    Nii = nifti(spm_select('FPList','/home/mbrud/dev/mbrud/code/matlab/Patient-Preprocessing/example-data/BrainWebManyN/','^.*\.nii$'));
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
opt.pth_mtv     = fullfile('/home','mbrud','dev','mbrud','code','matlab','MTV-preproc');
if TESTCASE == 2 || TESTCASE == 3
    opt.do.reslice = true;
    opt.do.vx      = true;
end
if TESTCASE == 4
    opt.do.res_orig = true;
end
if TESTCASE == 5 || TESTCASE == 7
    opt.do.denoise  = false;    
    opt.do.superres = true;
    opt.crop.neck   = true;
    if TESTCASE == 7
        opt.superres.ix = ix;
    end
%     opt.do.reslice = true;
%     opt.do.vx      = true;        
end
if TESTCASE == 6
    opt.crop.neck       = true;
    opt.write2d.axis_2d = 3; % 1. Sagittal 2. Coronal 3. Axial
    opt.do.write2d      = true;
end
% opt.do.write2d = true;
% opt.do.denoise = false;

%----------------------
% Do preprocessing
%----------------------

out = RunPreproc(Nii,opt);

if 0
    %----------------------
    % Go back to native space orientation
    %----------------------

    C = numel(out.pth.im);
    for c=1:C    
        Mc = spm_get_space(out.pth.im{c}); 
        spm_get_space(out.pth.im{c},out.mat{c}*Mc); 

        if ~isempty(out.pth.lab{c})
            Mc = spm_get_space(out.pth.lab{c}); 
            spm_get_space(out.pth.lab{c},out.mat{c}*Mc); 
        end
    end
end