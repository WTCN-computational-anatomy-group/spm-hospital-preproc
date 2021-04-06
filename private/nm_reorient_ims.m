function [Nii,M] = nm_reorient_ims(Nii)

fprintf('Reslicing voxel-to-world...')
N    = numel(Nii{1});
M    = cell(1,N);
M(:) = {eye(4)};
for n=1:N
    f         = Nii{1}(n).dat.fname;    
    f         = nm_reorient(f);
    Nii{1}(n) = nifti(f);
end
if numel(Nii) > 1
    % Keep labels in alignment
    for n=1:N
        if n > numel(Nii{2}) || isempty(Nii{2}(n).dat), continue; end
        
        f         = Nii{2}(n).dat.fname;    
        f         = nm_reorient(f,[],'',0);
        Nii{2}(n) = nifti(f);
    end    
end
fprintf('done!\n')
%==========================================================================
