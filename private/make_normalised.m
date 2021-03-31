function pth_norm = make_normalised(Nii,opt)

fprintf('Making normalised...')

N = numel(Nii{1});
    
%----------------------
% Potentially mask using labels
%----------------------

if opt.mask
    % If there are labels, use these to msk input images to segmentation
    buf = cell(1,N);
    for n=1:N
        if isempty(Nii{2}(n).dat) || n ~= opt.vol
            continue; 
        end
        
        mskn                 = Nii{2}(n).dat() > 0;
        im                   = Nii{1}(n).dat();
        buf{n}               = im(mskn);
        im(mskn)             = 0;
        Nii{1}(n).dat(:,:,:) = im;        
    end 
    
%     cnt = false(1,N);
%     buf = cell(1,N);
%     for n=1:N
%         if isempty(Nii{2}(n).dat), continue; end
%         
%         cnt(n)               = true;
%         mskn                 = Nii{2}(n).dat() > 0;
%         im                   = Nii{1}(n).dat();
%         buf{n}               = im(mskn);
%         im(mskn)             = 0;
%         Nii{1}(n).dat(:,:,:) = im;        
%     end 
%     
%     msk = false;
%     for n=1:N
%         if cnt(n)
%             mskn = Nii{2}(n).dat() > 0;
%             msk  = msk | mskn;
%         end        
%     end
%     
%     for n=1:N
%         if any(cnt == true) && ~cnt(n)
%             im                   = Nii{1}(n).dat();
%             buf{n}               = im(msk);
%             im(msk)              = 0;
%             Nii{1}(n).dat(:,:,:) = im;          
%         end        
%     end
end
    
%----------------------
% Make input for normalisataion routine
%----------------------

vol      = {Nii{1}(min(opt.vol,N)).dat.fname};
resample = {};
for n=1:N
    f = Nii{1}(n).dat.fname;
    resample{end + 1} = f;
    
    if isempty(Nii{2}(n).dat), continue; end
    f = Nii{2}(n).dat.fname;
    resample{end + 1} = f;
end
    
%----------------------
% Create normalisation job
%----------------------

job               = struct;
job.subj.vol      = vol;
job.subj.resample = resample';

pth = fileparts(which('spm'));
tpm = fullfile(pth,'tpm','TPM.nii');

job.eoptions.biasreg  = 1e-4;
job.eoptions.biasfwhm = 60;
job.eoptions.tpm      = {tpm};
job.eoptions.affreg   = 'mni';
job.eoptions.reg      = [0 1e-3 0.5 0.05 0.2];
job.eoptions.fwhm     = 0;
job.eoptions.samp     = 3;
job.eoptions.write    = opt.write_df;

job.woptions.bb     = spm_get_defaults('normalise.write.bb');
job.woptions.vox    = spm_get_defaults('normalise.write.vox');
job.woptions.interp = spm_get_defaults('normalise.write.interp');
job.woptions.prefix = spm_get_defaults('normalise.write.prefix');
    
%----------------------
% Run normalisation
%----------------------

pth_norm = my_spm_run_norm(job,Nii,buf,opt);

fprintf('done!\n')
%==========================================================================#

%==========================================================================#
function out = my_spm_run_norm(job,Nii,buf,opt)
% SPM job execution function
% takes a harvested job data structure and call SPM functions to perform
% computations on the data.
% Input:
% job    - harvested job data structure (see matlabbatch help)
% Output:
% out    - computation results, usually a struct variable.
%__________________________________________________________________________
% Copyright (C) 2005-2018 Wellcome Trust Centre for Neuroimaging

% $Id: spm_run_norm.m 7406 2018-08-21 17:29:53Z john $


for i=1:numel(job.subj)

    jobi      = job;
    jobi.subj = job.subj(i);

    %-Normalise: Estimate
    %----------------------------------------------------------------------
    if isfield(jobi,'eoptions')
        normalise(jobi);
    end

    N = numel(Nii{1});
    if opt.mask
        % If there are labels, reapply masked values    
        for n=1:N
            if isempty(buf{n}), continue; end

            mskn                 = Nii{2}(n).dat() > 0;
            im                   = Nii{1}(n).dat();        
            im(mskn)             = buf{n};
            Nii{1}(n).dat(:,:,:) = im;      
        end 
    end

    %-Normalise: Write
    %----------------------------------------------------------------------
    if isfield(jobi,'woptions')
        write_norm(jobi);
    end
end

%-Dependencies
%--------------------------------------------------------------------------
for i=1:numel(job.subj)
    if ~isfield(job.subj(i),'def'),
        out(i).def = {spm_file(char(job.subj(i).vol), 'prefix','y_', 'ext','.nii')};
    end
    
    if isfield(job,'woptions'),
        out(i).files = spm_file(job.subj(i).resample, 'prefix',job.woptions.prefix);
    end
end
%==========================================================================

%==========================================================================
function normalise(job)
% Estimate deformations via Segmentation
preproc8.channel.vols     = '<UNDEFINED>';
preproc8.channel.biasreg  = job.eoptions.biasreg;
preproc8.channel.biasfwhm = job.eoptions.biasfwhm;
preproc8.channel.write    = [0 0];

tpm = job.eoptions.tpm{:};
Nii = nifti(tpm);
for i=1:size(Nii.dat,4),
    preproc8.tissue(i) = struct('tpm',   {{[tpm ',' num2str(i)]}},...
                                'ngaus', Inf,...
                                'native',[0 0],...
                                'warped',[0 0]);
end
preproc8.warp.mrf    = 0;
preproc8.warp.reg    = job.eoptions.reg;
preproc8.warp.affreg = job.eoptions.affreg;
preproc8.warp.fwhm   = job.eoptions.fwhm;
preproc8.warp.samp   = job.eoptions.samp;
preproc8.warp.write  = job.eoptions.write;
preproc8.savemat     = 0;

for i=1:numel(job.subj),
    preproc8.channel.vols = job.subj(i).vol;
    spm_preproc_run(preproc8);
end
%==========================================================================

%==========================================================================
function write_norm(job)
% Write the spatially normalised data

defs.comp{1}.def         = '<UNDEFINED>';
defs.comp{2}.idbbvox.vox = job.woptions.vox;
defs.comp{2}.idbbvox.bb  = job.woptions.bb;
defs.out{1}.pull.fnames  = '';
defs.out{1}.pull.savedir.savesrc = 1;
defs.out{1}.pull.interp  = job.woptions.interp;
defs.out{1}.pull.mask    = 1;
defs.out{1}.pull.fwhm    = [0 0 0];
defs.out{1}.pull.prefix  = job.woptions.prefix;

for i=1:numel(job.subj)
    defs.out{1}.pull.fnames = job.subj(i).resample;
    if ~isfield(job.subj(i),'def')
        defs.comp{1}.def = {spm_file(char(job.subj(i).vol), 'prefix','y_', 'ext','.nii')};
    else
        defs.comp{1}.def = job.subj(i).def;
    end

    Nii = nifti(defs.comp{1}.def);
    vx  = sqrt(sum(Nii.mat(1:3,1:3).^2));
    if det(Nii.mat(1:3,1:3))<0, vx(1) = -vx(1); end

    o   = Nii.mat\[0 0 0 1]';
    o   = o(1:3)';
    dm  = size(Nii.dat);
    bb  = [-vx.*(o-1) ; vx.*(dm(1:3)-o)];

    defs.comp{2}.idbbvox.vox = job.woptions.vox;
    defs.comp{2}.idbbvox.bb  = job.woptions.bb;
    defs.comp{2}.idbbvox.vox(~isfinite(defs.comp{2}.idbbvox.vox)) = vx(~isfinite(defs.comp{2}.idbbvox.vox));
    defs.comp{2}.idbbvox.bb(~isfinite(defs.comp{2}.idbbvox.bb)) = bb(~isfinite(defs.comp{2}.idbbvox.bb));
    spm_deformations(defs);
end
%==========================================================================#
