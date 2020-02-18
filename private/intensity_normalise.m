function Nii = intensity_normalise(Nii,opt)

rng = opt.rng;

fprintf('Normalising intensities...')
N = numel(Nii{1});
for n=1:N
    im                   = Nii{1}(n).dat();    
    im                   = rescale(im,rng(1),rng(2));
    Nii{1}(n).dat(:,:,:) = im;
end
fprintf('done!\n')
%==========================================================================