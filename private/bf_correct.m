function Nii = bf_correct(Nii,pth_seg)

fprintf('Bias-field correcting...')
N        = numel(Nii{1});
for n=1:N
    delete(Nii{1}(n).dat.fname);
    Nii{1}(n) = nifti(pth_seg{7}{n});    
end
fprintf('done!\n')
%==========================================================================