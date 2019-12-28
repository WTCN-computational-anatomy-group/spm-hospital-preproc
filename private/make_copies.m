function Nii = make_copies(Nii,DirOut)
fprintf('Making copies...')

N = numel(Nii{1});
for n=1:N
    f           = Nii{1}(n).dat.fname;
    [~,nam,ext] = fileparts(f);
        
    nf = fullfile(DirOut,[nam ext]);
    if isfile(nf)
        delete(nf);
    end
    copyfile(f,DirOut);
        
    nf                  = fullfile(DirOut,[nam ext]);
    Nii{1}(n).dat.fname = nf;
    
    Nii{1}(n).dat.permission = 'rw';
end

if numel(Nii) > 1
    % Copy labels too
    for n=1:N
        if isempty(Nii{2}(n).dat), continue; end
        
        f           = Nii{2}(n).dat.fname;
        [~,nam,ext] = fileparts(f);
        nf          = fullfile(DirOut,[nam ext]);

        im = Nii{2}(n).dat();
        
        if max(im(:)) > 255, error('Copy labels: values over 255!'); end
        
        create_nii(nf,im,Nii{2}(n).mat,spm_type('uint8'),'labels');        
        
        Nii{2}(n) = nifti(nf);
%         Nii{2}(n).dat.permission = 'rw';
    end    
end
fprintf('done!\n')
%==========================================================================