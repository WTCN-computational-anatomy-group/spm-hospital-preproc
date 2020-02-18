function Nii = skull_strip(Nii,pth_seg)

fprintf('Skull-stripping...')

% Make mask
msk = 0;
for k=1:3
    Niik = nifti(pth_seg{1}{k});
    im   = Niik.dat();
    msk  = msk + im;
end
msk = msk > 0.5;
    
% Apply mask
N = numel(Nii{1});
for n=1:N
    im                   = Nii{1}(n).dat();    
    im(~msk)             = 0;    
    Nii{1}(n).dat(:,:,:) = im;
end
fprintf('done!\n')
%==========================================================================