function out = RunPreproc(paths,opt)
% Some basic preprocessing of hospital neuroimaging data.
%
% INPUT
% paths - Can be given as:
%            paths = im.ext, 
%            paths = {im1.ext, ..., imN.ext}
%            paths = {{im1.ext, ..., imN.ext}, 
%                     {lab1.ext, ..., labN.ext}}
%         where valid extensions are .nii and .nii.gz.          
% opt   - Preprocessing options
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

% Set options
if nargin < 2, opt = struct; end
opt = get_default_opt(opt);

% Check MATLAB path
if isempty(fileparts(which('spm')))      
    error('SPM12 not on the MATLAB path!');  % download from https://www.fil.ion.ucl.ac.uk/spm/software/download/
end
if opt.do.denoise || opt.do.superres
    if isempty(fileparts(which('spm_superres')))
        error('spm_superres not on the MATLAB path!');  % download/clone from https://github.com/brudfors/spm_superres
    end
end

if ~(exist(opt.dir_out,'dir') == 7)  
    % Create output directory
    mkdir(opt.dir_out);  
end

% Make sure output directory is encoded by its full path
s           = what(opt.dir_out);
opt.dir_out = s.path;

% Copy (so to not overwrite originals)
[Nii,was_gz,nams] = read_and_copy(paths,opt.dir_out);

C = numel(Nii{1});

if numel(Nii) > 1
    % Collapse labels
    Nii = collapse_labels(Nii,opt.labels.part);    
end

if opt.do.erode
    % Remove a few of the outer voxels
    Nii = erode_im(Nii);
end

M    = cell(1,C);
M(:) = {eye(4)};
if opt.do.res_orig
    % Reset origin (important for CT)
    vx = [];
    if ~opt.do.superres && opt.do.vx
        vx = opt.vx.size;
    end
    [Nii,M] = reset_origin(Nii, opt.reset_orig.only_neg, vx);
end

if opt.do.real_mni
    % Realing to MNI space
    [Nii,M] = realign2mni(Nii,M,opt.realign2mni);
end

if opt.do.nm_reorient
    % Reslice so that image data is in world space
    [Nii,M] = nm_reorient_ims(Nii); % M is set to identity
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
    Nii = superres(Nii,opt.do.coreg,opt.superres);
    
    % Coreg (one more time after super-resolving)
    Nii = coreg(Nii,opt.coreg);
else
    if opt.do.reslice
        % Make images same dimensions
        [Nii,M] = reslice_images(Nii,M,opt.reslice);
    end
    
    if opt.do.vx && ~opt.do.res_orig && ~opt.do.bb_spm
        % Set same voxel size        
        Nii = resample_images(Nii,opt.vx);
    end
end

if opt.do.bb_spm
    % Crop to SPM12 BB
    if opt.do.vx
        vx = opt.vx.size;
    else
        vx = [];
    end
    Nii = apply_bb(Nii, vx);
end

pth_seg = {};
if opt.do.segment
    % Run SPM12 segmentation
    pth_seg = segment_preproc8(Nii,opt.segment);
end

if opt.do.bfcorr
    % Bias field correct (depends on segment_preproc8())
    Nii = bf_correct(Nii,pth_seg);
end

if opt.do.skullstrip
    % Skull-strip (depends on segment_preproc8())
    Nii = skull_strip(Nii,pth_seg);
    for k=1:numel(pth_seg{1}), delete(pth_seg{1}{k}); end
end

if ~isempty(opt.pth_template) && isfile(opt.pth_template)
    % Reslice and affinely register images to a template
    [Nii,M] = reslice2template(Nii,M,opt.pth_template);
end

if numel(Nii) > 1 && (opt.do.nm_reorient || opt.do.crop || opt.do.superres || opt.do.vx || (opt.do.vx && opt.do.res_orig) || opt.do.reslice || (~isempty(opt.pth_template) && isfile(opt.pth_template)))
    % Reslice labels
    Nii = reslice_labels(Nii,opt.reslice);
end

pth_norm = {};
if opt.do.normalise
    % Create normalised versions of Nii
    pth_norm = make_normalised(Nii,opt.normalise);
end

if opt.do.int_norm
    % Normalise image intensities in the range opt.int_norm.rng
    Nii = intensity_normalise(Nii,opt.int_norm);
end

P2d = {};
if opt.do.write2d
    % Write 2D versions
    P2d = write_2d(Nii,pth_seg,opt.dir_out2d,opt.write2d);
end

% Allocate output
C             = numel(Nii{1});
out           = struct;
out.pth.im    = cell(1,C);
out.pth.im2d  = cell(1,C);
out.pth.lab   = cell(1,C);
out.pth.lab2d = cell(1,C);
out.pth.seg   = {}; 
out.pth.norm  = {}; 
out.mat       = cell(1,C);
for i=1:2
    for c=1:C
        if (i == 2 && numel(Nii) == 1) || isempty(Nii{i}(c).dat), continue; end
        
        f = Nii{i}(c).dat.fname;
        
        if i == 1
            out.mat{c} = M{c};
        end  
        
        if i == 1
            out.pth.im{c} = Nii{i}(c).dat.fname;
            
            [pth,~,ext] = fileparts(out.pth.im{c});
            nfname = fullfile(pth, [nams{1, c} ext]);
            movefile(out.pth.im{c}, nfname);
            out.pth.im{c} = nfname;
            
            if opt.do.go2native
                p = out.pth.im{c};    
                Mc = spm_get_space(p); 
                spm_get_space(p,M{c}*Mc); 
            end
            
            if was_gz == true
                % Compress back to .gz
                p = gzip(out.pth.im{c});
                delete(out.pth.im{c});
                out.pth.im{c} = p{1};
            end
        else
            out.pth.lab{c} = Nii{i}(c).dat.fname;
            
            [pth,~,ext] = fileparts(out.pth.lab{c});
            nfname = fullfile(pth, [nams{2, c} ext]);
            movefile(out.pth.lab{c}, nfname);
            out.pth.lab{c} = nfname;
            
            if opt.do.go2native
                p = out.pth.lab{c};    
                Mc = spm_get_space(p); 
                spm_get_space(p,M{c}*Mc); 
            end
            
            if was_gz == true
                % Compress back to .gz
                p = gzip(out.pth.lab{c});
                delete(out.pth.lab{c});
                out.pth.lab{c} = p{1};                
            end
        end
                            
        if ~isempty(pth_seg)
            out.pth.seg = pth_seg;
                                  
            if opt.do.go2native
                for i1=1:numel(out.pth.seg{1})
                    p = deblank(out.pth.seg{1}{i1});    
                    Mc = spm_get_space(p); 
                    spm_get_space(p,M{c}*Mc); 
                end
            end
            
            if was_gz == true
                % Compress back to .gz
                p = gzip(out.pth.seg);
                delete(out.pth.seg);
                out.pth.seg = p{1};                
            end
        end
                            
        if ~isempty(pth_norm)
            out.pth.norm = pth_norm;
            if was_gz == true
                % Compress back to .gz
                p = gzip(out.pth.norm);
                delete(out.pth.norm);
                out.pth.norm = p{1};                
            end
        end
        
        if ~isempty(P2d)
            if i == 1
                out.pth.im2d{c}  = P2d{i}{c};
                
                [pth,~,ext] = fileparts(out.pth.im2d{c});
                nfname = fullfile(pth, [nams{1, c} ext]);
                movefile(out.pth.im2d{c}, nfname);
                out.pth.im2d{c} = nfname;
                
                if was_gz == true
                    % Compress back to .gz
                    p = gzip(oout.pth.im2d{c});
                    delete(out.pth.im2d{c});
                    out.pth.im2d{c} = p{1};                
                end
            else
                out.pth.lab2d{c} = P2d{i}{c};
                
                [pth,~,ext] = fileparts(out.pth.lab2d{c});
                nfname = fullfile(pth, [nams{2, c} ext]);
                movefile(out.pth.lab2d{c}, nfname);
                out.pth.lab2d{c} = nfname;
                
                if was_gz == true
                    % Compress back to .gz
                    p = gzip(out.pth.lab2d{c});
                    delete(out.pth.lab2d{c});
                    out.pth.lab2d{c} = p{1};                
                end
            end
        end
        
        if opt.do.writemat
            [pth,nam] = fileparts(f);   
            nP        = fullfile(pth,['mat' nam '.mat']);
            Mc        = M{c};
            save(nP,'Mc')
        end
    end
end
%==========================================================================