function Nii = collapse_labels(Nii,part)
if isempty(part), return; end

N = numel(Nii{1});
if numel(Nii) > 1
    fprintf('Collapsing labels...')
        
    
    for n=1:N
        if n > numel(Nii{2}) || isempty(Nii{2}(n).dat), continue; end
        
        nii     = Nii{2}(n).dat();
        nlabels = do_collapse(nii,part);
        
        Nii{2}(n).dat(:,:,:) = uint8(nlabels);
    end    
    
    fprintf('done!\n')
end
%==========================================================================

%==========================================================================
function nlabels = do_collapse(labels,part)
K       = numel(part);
labels  = round(labels);
if ischar(part)
    p       = unique(labels);
    nlabels = ismember(labels,p(2:end));
else
    nlabels = zeros(size(labels));
    for k=1:K
        if iscell(part)
            p = part{k}';
        else
            p = part(k)';
        end
        msk          = ismember(labels,p);
        nlabels(msk) = k;
    end
end
%==========================================================================