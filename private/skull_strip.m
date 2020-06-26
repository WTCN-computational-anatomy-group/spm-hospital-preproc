function Nii = skull_strip(Nii,pth_seg)

fprintf('Skull-stripping...')

% Make mask
gwc = cell(1,3);
for k=1:3
    Niik = nifti(pth_seg{1}{k});
    gwc{k} = Niik.dat();
end

% Build brain mask
msk = gwc{1} > 0.1 | gwc{2} > 0.1 | gwc{3} > 0.999;
r = 3;
msk = imdilate(msk, strel('sphere',r));
msk = imerode(msk, strel('sphere',r));
% msk = imerode(msk, strel('sphere',1));
msk = imfill(msk, 'holes');
% Get largest connected component
C = bwlabeln(msk);
unqC = unique(C);
vls = zeros([1 numel(unqC)]);
for i=1:numel(unqC)
    vls(i) = sum(sum(sum(C == unqC(i))));
end
[~,ix] = sort(vls);
ix = (ix(end - 1) - 1);
msk = C == ix;

if false
    % For testing
    im = Nii{1}(1).dat(); 
    mx = 0.5*max(im(:));
    dm = size(im);
    figure(11)
    subplot(331)
    imagesc(im(:,:,round(0.5*dm(3))), [0 mx])
    axis off
    subplot(334)
    imagesc(squeeze(im(:,round(0.5*dm(2)),:)), [0 mx])
    axis off
    subplot(337)
    imagesc(squeeze(im(round(0.5*dm(1)),:,:)), [0 mx])
    axis off
    colormap(gray)
    mim = msk.*im;
    subplot(332)
    imagesc(mim(:,:,round(0.5*dm(3))), [0 mx])
    axis off
    subplot(335)
    imagesc(squeeze(mim(:,round(0.5*dm(2)),:)), [0 mx])
    axis off
    subplot(338)
    imagesc(squeeze(mim(round(0.5*dm(1)),:,:)), [0 mx])
    axis off
    colormap(gray)
    subplot(333)
    imagesc(msk(:,:,round(0.5*dm(3))))
    axis off
    subplot(336)
    imagesc(squeeze(msk(:,round(0.5*dm(2)),:)))
    axis off
    subplot(339)
    imagesc(squeeze(msk(round(0.5*dm(1)),:,:)))
    axis off
    colormap(gray)
end

% Apply mask
N = numel(Nii{1});
for n=1:N
    im                   = Nii{1}(n).dat();    
    im                   = im.*msk;
    Nii{1}(n).dat(:,:,:) = im;
end
fprintf('done!\n')
%==========================================================================