% Demo runs of the RunPreproc function. A couple of different test-cases
% are provided.
%
% OBS: To run, replace the hardcoded paths with your own data
%_______________________________________________________________________
%  Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

clear;

%----------------------
% TODO
% Skull-strip
% Warp intensity image
% Write 2D segmentations
%----------------------

%----------------------
% Pick a test-case:
% 1. MRI single-channel with 2D version
% 2. MRI multi-channel with labels
% 3. CT w. labels
% 4. Hospital MRI superres
% 5. MRI super-res with cell array input
% 6. MRI super-res with more than one subject
%----------------------

TESTCASE = 1;
DirInOut = ['TESTCASE-' num2str(TESTCASE)];

% Options to the preprocessing code
opt         = struct;
opt.pth_mtv = fullfile('/home','mbrud','dev','mbrud','code','matlab','MTV-preproc');
opt.dir_out = fullfile('output',DirInOut);   
if exist(opt.dir_out,'dir') == 7    
    rmdir(opt.dir_out,'s');
end

% Define input(s) and test-case specific options
dat = struct;
if     TESTCASE == 1
    % MRI single-channel
    dat.Nii = nifti(fullfile('example-data',DirInOut,'c0003s0008t01.nii'));
    
    % Options    
    opt.do.real_mni     = true;        
    opt.do.crop         = true;   
    opt.crop.neck       = true;
    opt.do.denoise      = true;
    opt.write2d.axis_2d = 3;    % 1. Sagittal 2. Coronal 3. Axial
    opt.do.write2d      = true;
    
    opt.do.segment                  = true;
    opt.segment.write_tc            = false(6,4);
    opt.segment.write_tc(1:3,[3 4]) = true;
    opt.segment.write_df            = [true true];
    opt.segment.write_bf            = [true true];
elseif TESTCASE == 2        
    % MRI multi-channel w. labels
    dat.Nii       = cell(1,2);
    dat.Nii{1}    = nifti;
    dat.Nii{1}(1) = nifti(fullfile('example-data',DirInOut,'1-T1.nii'));
    dat.Nii{1}(2) = nifti(fullfile('example-data',DirInOut,'1-FLAIR.nii'));    
    dat.Nii{1}(3) = nifti(fullfile('example-data',DirInOut,'1-IR.nii'));    
    dat.Nii{2}    = nifti; % Map label image to intensity image it was annotated on
    dat.Nii{2}(1) = nifti;
    dat.Nii{2}(2) = nifti(fullfile('example-data',DirInOut,'1-segm.nii'));
    dat.Nii{2}(3) = nifti;
    
    % Options     
    opt.do.real_mni = true;        
    opt.do.crop     = true; 
    opt.do.reslice  = true;
    opt.reslice.ref = 1;
    opt.do.vx       = true;  
    opt.labels.part = {[1],[2],[3 4],[5],[6],[7],[8]};
elseif TESTCASE == 3
    % CT w. labels
    dat.Nii    = cell(1,2);
    dat.Nii{1} = nifti(fullfile('example-data',DirInOut,'sCROMIS2ICH_01006-0004-00002-000001.nii'));
    dat.Nii{2} = nifti(fullfile('example-data',DirInOut,'sCROMIS2ICH_01006-0004-00002-000001_smask.nii'));
    
    % Options     
    opt.do.res_orig = true;
    opt.do.real_mni = true;        
    opt.do.crop     = true;     
    opt.crop.neck   = true;
    opt.do.vx       = true; 
    opt.vx.size     = 1;                
elseif TESTCASE == 4      
    % Hospital MRI superres    
    dat.Nii    = nifti; 
    dat.Nii(1) = nifti(fullfile('example-data',DirInOut,'0002-00001-000021-01.nii'));
    dat.Nii(2) = nifti(fullfile('example-data',DirInOut,'0003-00001-000001-01.nii'));    
    dat.Nii(3) = nifti(fullfile('example-data',DirInOut,'0004-00001-000001-01.nii'));        
    
    % Options     
    opt.do.real_mni = true;        
    opt.do.crop     = true; 
    opt.crop.neck   = true;
    opt.do.superres = true;        
elseif TESTCASE == 5
    % MRI super-res with cell array input    
    dat.Nii = nifti(spm_select('FPList',fullfile('example-data',DirInOut),'^.*\.nii$'));
        
    % Options     
    opt.do.real_mni = true;        
    opt.do.crop     = true; 
    opt.crop.neck   = true;
    opt.do.superres = true;       
    opt.superres.ix = [1 2 1 2 3];
elseif TESTCASE == 6
    % Hospital MRI superres w. more than one subject
    dat(1).Nii    = nifti; 
    dat(1).Nii(1) = nifti(fullfile('example-data',DirInOut,'1282601181791321211130022749_Flair.nii'));
    dat(1).Nii(2) = nifti(fullfile('example-data',DirInOut,'1282601181791321211130022749_T1.nii'));    
    dat(1).Nii(3) = nifti(fullfile('example-data',DirInOut,'1282601181791321211130022749_T2.nii'));        
    dat(2).Nii    = nifti; 
    dat(2).Nii(1) = nifti(fullfile('example-data',DirInOut,'1282601181791321211130083960_Flair.nii'));
    dat(2).Nii(2) = nifti(fullfile('example-data',DirInOut,'1282601181791321211130083960_T1.nii'));    
    dat(2).Nii(3) = nifti(fullfile('example-data',DirInOut,'1282601181791321211130083960_T2.nii'));            
    
    % Options     
    opt.do.real_mni = true;        
    opt.do.crop     = true; 
    opt.crop.neck   = true;
    opt.do.superres = true;       
else
    error('Undefined test-case!')
end

%----------------------
% Do preprocessing
%----------------------

for s=1:numel(dat) % loop over subjects
    % Run preprocessing
    out = RunPreproc(dat(s).Nii,opt);

    % Go back to native space orientation
    go2native(out);
end