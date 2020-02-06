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
fprintf('done!\n')
%==========================================================================