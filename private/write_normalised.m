function pth_norm = write_normalised(Nii,opt,pth_seg)

fprintf('Writing normalised...')

N = numel(Nii{1});

if isempty(opt.vox)
    opt.vox = spm_get_defaults('normalise.write.vox');
elseif numel(opt.vox) == 1
    opt.vox = opt.vox*ones(1,3);
end

def = pth_seg{5}{1};
if ~exist(def, 'file')
    error('Cannot find forward deformation field: %s\n', def)
end

%----------------------
% Create normalisation job
%----------------------

job               = struct;
job.subj.vol      = [];
job.subj.def      = pth_seg{5}{1};
job.woptions.bb     = spm_get_defaults('normalise.write.bb');
job.woptions.vox    = opt.vox;
job.woptions.interp = spm_get_defaults('normalise.write.interp');
job.woptions.prefix = spm_get_defaults('normalise.write.prefix');

%----------------------
% Make input for normalisataion routine
%----------------------

% Images using default interpolation
images = cell(N,1);
[images{:}] = Nii{1}(:).dat.fname;
job.subj.resample = images;
write_norm(job);

% Labels using nearest neighbour
labels = {};
if numel(Nii) > 1 || ~isempty(Nii{2})
    labels = cell(numel(Nii{2}),1);
    [labels{:}] = Nii{2}(:).dat.fname;
    job.subj.resample = labels;
    job.woptions.interp = 0;
    write_norm(job);
end

% Return filenames
pth_norm.def = {def};
pth_norm.img = spm_file(images, 'prefix',job.woptions.prefix);
if ~isempty(labels)
    pth_norm.lab = spm_file(labels, 'prefix',job.woptions.prefix);
end

%==========================================================================
fprintf('done!\n')

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
