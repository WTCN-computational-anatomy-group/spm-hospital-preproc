function pths = segment_preproc8(Nii,opt)

fprintf('Segmenting...')

write_tc = opt.write_tc;
write_bf = opt.write_bf;
write_df = opt.write_df;
dir_out  = opt.dir_out;
make_4d  = opt.make_4d;

if ~isempty(dir_out) && ~isfolder(dir_out), mkdir(dir_out); end

K = 6;             % Number of tissue classes
N = numel(Nii{1}); % Number of channels
V = spm_vol;
for n=1:N
    V(n) = spm_vol(Nii{1}(n).dat.fname);
end

% Write options
if isempty(write_tc)
    write_tc            = false(K,4);
    write_tc(1:3,[1 3]) = true;
end
if isempty(write_bf)
    write_bf = false(N,2);
end
if numel(write_bf) == 1
    write_bf = repmat(write_bf,1,2);
end
if size(write_bf,1) == 1
    write_bf = repmat(write_bf,N,1);
end
if isempty(write_df)
    write_df = false(1,2);
end
if numel(write_df) == 1
    write_df = repmat(write_df,1,2);
end

% Segmentation options
obj          = struct;
obj.bb       = NaN(2,3);
obj.vox      = NaN;
obj.affreg   = 'mni';
obj.reg      = [0 0.001 0.5 0.05 0.2]*1e-1;
obj.fwhm     = 1;
obj.samp     = 3;
obj.biasfwhm = 60*ones(1,N);
obj.mrf      = 2;    
obj.cleanup  = 1;

% Less bias field if CT
dat  = single(Nii{1}(1).dat.fname);
isct = min(dat(:)) < 0;
if isct
    obj.biasreg  = 10; 
else
    obj.biasreg  = 1e-3*(1/5)*ones(1,N);    
end
clear dat

% Load atlas
tpmname   = fullfile(spm('dir'),'tpm','TPM.nii');
obj.lkp   = [1 1 2 2 3 3 4 4 5 5 5 6 6];
obj.tpm   = spm_load_priors8(tpmname);
obj.image = V;

% Register atlas
M                       = obj.image(1).mat;
c                       = (obj.image(1).dim+1)/2;
obj.image(1).mat(1:3,4) = -M(1:3,1:3)*c(:);
[Affine1,ll1]           = spm_maff8(obj.image(1),8,(0+1)*16,obj.tpm,[],obj.affreg); % Closer to rigid
Affine1                 = Affine1*(obj.image(1).mat/M);

% Run using the origin from the header
obj.image(1).mat = M;
[Affine2,ll2]    = spm_maff8(obj.image(1),8,(0+1)*16,obj.tpm,[],obj.affreg); % Closer to rigid

% Pick the result with the best fit
if ll1>ll2, obj.Affine  = Affine1; else obj.Affine  = Affine2; end

% Initial affine registration.
obj.Affine = spm_maff8(obj.image(1),4,(obj.fwhm+1)*16,obj.tpm, obj.Affine, obj.affreg); % Closer to rigid
obj.Affine = spm_maff8(obj.image(1),3, obj.fwhm,      obj.tpm, obj.Affine, obj.affreg);    

% Run the actual segmentation
res = spm_preproc8(obj);

% Final iteration, so write out the required data.
spm_preproc_write8(res,write_tc,repmat(write_bf,N,1),write_df,obj.mrf,obj.cleanup,obj.bb,obj.vox,dir_out);

% Get paths to segmentations
if isempty(dir_out)
    dir_out = fileparts(V(1).fname);
end
   
[~,nam] = fileparts(V(1).fname);
pths    = cell(1,8);
prefix  = {['c[1-6]' nam],['rc[1-6]' nam],['wc[1-6]' nam], ...
           ['mwc[1-6]' nam],['y_' nam],['iy_' nam],['m' nam],['BiasField_' nam]};
for i=1:numel(prefix)
    files   = spm_select('List',dir_out,['^' prefix{i} '.*\.nii$']);
    if ~isempty(files)
        pths{i} = cellstr(files);
        for i1=1:numel(pths{i})
            pths{i}{i1} = fullfile(dir_out,pths{i}{i1});        
        end
    end
end

K = numel(pths{1});
if K> 0 && K < 6
    % Make background class
    Nii_seg = nifti;
    for k=1:K
        Nii_seg(k) = nifti(pths{1}{k});
    end
        
    bg = 1;    
    for k=1:K        
        im = single(Nii_seg(k).dat());
        bg = bg - im;
    end
    bg(bg < 0) = 0;
    
    f             = Nii_seg(1).dat.fname;  
    [pth,nam,ext] = fileparts(f);
    
    nf = fullfile(pth,['c' num2str(K + 1) nam(3:end) ext]);  
    create_nii(nf,bg,Nii_seg(1).mat,Nii_seg(1).dat.dtype,['Tissue class ' num2str(K + 1)], ...
               Nii_seg(1).dat.offset,Nii_seg(1).dat.scl_slope,Nii_seg(1).dat.scl_inter);
           
    pths{1}{end + 1} = nf;
end


if make_4d
    % Write 4D nifti with segmentations
    K = numel(pths{1});
    if K > 0
        im = [];
        for k=1:K       
            Nii_seg = nifti(pths{1}{k});
            imk     = single(Nii_seg.dat());
            im      = cat(4,im,imk);
        end
        im(im < 0) = 0;
        im         = im./(sum(im,4) + eps('single'));

        f             = Nii_seg(1).dat.fname;  
        [pth,nam,ext] = fileparts(f);

        nf = fullfile(pth,['c0' nam(3:end) ext]);  
        create_nii(nf,im,Nii_seg(1).mat,Nii_seg(1).dat.dtype,'Tissue classes (4D)', ...
                   Nii_seg(1).dat.offset,Nii_seg(1).dat.scl_slope,Nii_seg(1).dat.scl_inter);           
    end
end

fprintf('done!\n')
%==========================================================================