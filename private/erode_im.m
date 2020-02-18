function Nii = erode_im(Nii)
fprintf('Eroding...')
N = numel(Nii{1});
for n=1:N
    im = Nii{1}(n).dat();    
    
    % Make binary mask
    bmsk = im ~= 0 & isfinite(im);
    
    % Erode
    se    = strel('sphere',3);
    ebmsk = imerode(bmsk, se);
    ebmsk = ebmsk > 0;
    
%     im1 = bmsk + ebmsk;    
%     figure(11)
%     imagesc(im1(:,:,round(size(im1,3)/2)))
    
    % Update image data
    im(~ebmsk)           = 0;    
    Nii{1}(n).dat(:,:,:) = im;   
end
fprintf('done!\n')
%==========================================================================