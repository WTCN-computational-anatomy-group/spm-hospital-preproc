function Nii = make_copies(Nii,DirOut)
fprintf('Making copies...')

N = numel(Nii{1});
for n=1:N
    f       = Nii{1}(n).dat.fname;
    [~,nam] = fileparts(f);        
    nf      = fullfile(DirOut,[nam '.nii']);
    if isfile(nf), delete(nf); end   
            
    create_nii(nf,Nii{1}(n).dat(),Nii{1}(n).mat,Nii{1}(n).dat.dtype,'copy',...
               Nii{1}(n).dat.offset,Nii{1}(n).dat.scl_slope,Nii{1}(n).dat.scl_inter);
               
    Nii{1}(n) = nifti(nf);    
end

if numel(Nii) > 1
    % Copy labels too
    for n=1:N
        if isempty(Nii{2}(n).dat), continue; end
        
        f       = Nii{2}(n).dat.fname;
        [~,nam] = fileparts(f);
        nf      = fullfile(DirOut,[nam '.nii']);
        if isfile(nf), delete(nf); end

        im = Nii{2}(n).dat();        
        if max(im(:)) > 255, error('Copy labels: values over 255!'); end
        
        create_nii(nf,im,Nii{2}(n).mat,spm_type('uint8'),'labels');        
        
        Nii{2}(n) = nifti(nf);
    end    
end
fprintf('done!\n')
%==========================================================================