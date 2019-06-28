% Demo runs of the RunPreproc function
%_______________________________________________________________________
%  Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

clear;

%----------------------
% Read image data
%----------------------

% One channel example
P = '/home/mbrud/Data/Mo/best_session/001/003_T1_FL2D_SAG_AA/0005-00001-000001-01.nii';

% More than one channel example
% P = char({'/home/mbrud/dev/mbrud/code/matlab/preprocessing-code/1282601181791321211150584332_T1.nii',...
%           '/home/mbrud/dev/mbrud/code/matlab/preprocessing-code/1282601181791321211150584332_T2.nii',...
%           '/home/mbrud/dev/mbrud/code/matlab/preprocessing-code/1282601181791321211150584332_Flair.nii'});
        
Nii = nifti(P);

%----------------------
% Set preprocessing options
%----------------------

opt.do.real_mni = true;
opt.do.coreg    = true;
opt.do.res_orig = false;
opt.do.denoise  = true;
opt.do.crop     = true;
opt.pth_mtv     = '/home/mbrud/dev/mbrud/code/matlab/MTV-preproc';

% If you are using more than once channel, you should uncomment the below
% lines
% opt.do.reslice = true;
% opt.do.vx      = true;

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