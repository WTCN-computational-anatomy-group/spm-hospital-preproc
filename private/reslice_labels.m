function Nii = reslice_labels(Nii,opt)

ref_ix = opt.ref;

N = numel(Nii{1});
if numel(Nii) > 1
    fprintf('Reslicing labels...')
    
    % Reslice labels too
    V = spm_vol(Nii{1}(ref_ix).dat.fname);
    for n=1:N
        if n > numel(Nii{2}) || isempty(Nii{2}(n).dat), continue; end
        
        Nii{2}(n) = do_reslice(V,Nii{2}(n));
    end    
    
    fprintf('done!\n')
end
%==========================================================================

%==========================================================================
function oNii = do_reslice(VRef,NiiLabels)
Mn  = VRef(1).mat;
dmn = VRef(1).dim;
Ml  = NiiLabels(1).mat;
lab = round(single(NiiLabels(1).dat()));
fnl = NiiLabels(1).dat.fname;

M = Ml\Mn;
y = Affine(dmn,M);

ul   = unique(lab);
ul   = ul';
l12  = zeros(dmn(1:3),'single');
p1   = zeros(size(l12),'single');
filt = [0.125 0.75 0.125];
for l=ul
    g0        = single(lab == l);
    g0        = convn(g0,reshape(filt,[3,1,1]),'same');
    g0        = convn(g0,reshape(filt,[1,3,1]),'same');
    g0        = convn(g0,reshape(filt,[1,1,3]),'same');
    tmp       = spm_diffeo('bsplins',g0,y,[1 1 1 0 0 0]);
    msk1      = (tmp>p1);
    p1(msk1)  = tmp(msk1);
    l12(msk1) = l;
end     
% Overwrite image data
VO             = spm_vol(fnl);
VO.dim(1:3)    = dmn;        
VO.mat         = Mn;
VO             = spm_create_vol(VO);        
Nii            = nifti(VO.fname);    
Nii.dat(:,:,:) = round(l12); 
% Return label nifti
oNii = nifti(fnl);
%==========================================================================

%==========================================================================
% Affine()
function psi0 = Affine(d,Mat)
id    = Identity(d);
psi0  = reshape(reshape(id,[prod(d) 3])*Mat(1:3,1:3)' + Mat(1:3,4)',[d 3]);
if d(3) == 1, psi0(:,:,:,3) = 1; end
%==========================================================================

%==========================================================================
% Identity()
function id = Identity(d)
id = zeros([d(:)',3],'single');
[id(:,:,:,1),id(:,:,:,2),id(:,:,:,3)] = ndgrid(single(1:d(1)),single(1:d(2)),single(1:d(3)));
%==========================================================================