function Nii = reslice_labels(Nii,opt)

ref_ix = opt.ref;

N = numel(Nii{1});
if numel(Nii) > 1
    fprintf('Reslicing labels...')
    
    % Reslice labels too
    V = spm_vol(Nii{1}(ref_ix).dat.fname);
    for n=1:N
        if isempty(Nii{2}(n).dat), continue; end
        
        Nii{2}(n) = do_reslice(V,Nii{2}(n));
    end    
    
    fprintf('done!\n')
end
%==========================================================================

%==========================================================================
function Niio = do_reslice(V_ref,Nii)
% Parameters
deg = 1;
dt  = [spm_type('uint8') spm_platform('bigend')];

% Get labels, etc
labels        = Nii.dat(:,:,:);
lkp           = unique(labels)';
K             = numel(lkp); % Number of labels
fname         = Nii.dat.fname;
[pth,nam,ext] = fileparts(fname);

fnames   = {V_ref.fname, Nii.dat.fname}';
Vo       = run_reslice(fnames,0);
Niio     = nifti(Vo(2).fname);
delete(Nii.dat.fname);

% % Iterate over each label and create a resliced label image (nlabels)
% dm      = V_ref.dim;
% labelso = zeros([dm K],'single');
% cnt     = 1;
% for k=lkp
%     labels_k = single(labels == k);
% 
%     % Smooth
%     labels_k = convn(labels_k,reshape([0.25 0.5 0.25],[3,1,1]),'same');
%     labels_k = convn(labels_k,reshape([0.25 0.5 0.25],[1,3,1]),'same');
%     labels_k = convn(labels_k,reshape([0.25 0.5 0.25],[1,1,3]),'same');
%                 
%     fname_k = fullfile(pth,['n' nam ext]);    
%     create_nii(fname_k,labels_k,Nii.mat,dt,'Resliced labels');        
%         
%     % Reslice
%     fnames   = {V_ref.fname, fname_k}';
%     Vo       = run_reslice(fnames,deg);
%     labels_k = single(Vo(2).private.dat(:,:,:));
% 
%     labelso(:,:,:,cnt) = cat(4,labels_k);
%     
%     delete(fname_k);
%     
%     cnt = cnt + 1;
% end
% clear labels labels_k
% delete(fname);
% 
% % Get MLs of resliced labels
% [~,ml] = max(labelso,[],4);
% ml     = ml - 1;
% 
% % Write output
% Niio            = nifti(Vo(2).fname);
% Niio.dat(:,:,:) = uint8(ml);
%==========================================================================