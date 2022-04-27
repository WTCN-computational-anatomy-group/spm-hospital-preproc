function [Nii,M] = reset_origin(Nii,opt,vx)
if nargin < 3, vx = [];    end

only_neg = opt.only_neg;
deg      = opt.deg;


fprintf('Resetting origin...')
N    = numel(Nii{1});
M    = cell(1,N);
M(:) = {eye(4)};
for n=1:N
    if only_neg && min(Nii{1}(n).dat(:)) >= 0
        % No negative values, do not reset
        continue
    end
    f = Nii{1}(n).dat.fname;    
    f = nm_reorient(f,vx,'ro',deg);
    
    M{n} = do_reset_origin(f);
    
    Nii{1}(n) = nifti(f);
end

if numel(Nii) > 1
    % Keep labels in alignment
    for n=1:N
        if n > numel(Nii{2}) || isempty(Nii{2}(n).dat), continue; end
        
        if only_neg && min(Nii{1}(n).dat(:)) >= 0
            % No negative values, do not reset
            continue
        end
    
        f = Nii{2}(n).dat.fname;    
        f = nm_reorient(f,vx,'ro',0);
        do_reset_origin(f);

        Nii{2}(n) = nifti(f);
    end    
end
fprintf('done!\n')
%==========================================================================


%==========================================================================
function Mout = do_reset_origin(pth,orig)
if nargin < 2, orig = []; end

V   = spm_vol(pth);
M   = V.mat;
dim = V.dim;
vx  = sqrt(sum(M(1:3,1:3).^2));

if det(M(1:3,1:3))<0
    vx(1) = -vx(1); 
end

if isempty(orig)
    orig = (dim(1:3)+1)/2;
end

off  = -vx.*orig;
M1   = [vx(1) 0      0         off(1)
           0      vx(2) 0      off(2)
           0      0      vx(3) off(3)
           0      0      0      1];

V    = spm_vol(pth);
M0   = V.mat;
Mout = M0/M1;

spm_get_space(pth,M1);   
%==========================================================================