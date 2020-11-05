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
lab = single(NiiLabels(1).dat());
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
l12 = uint8(l12);

[pth,nam,ext] = fileparts(fnl);
nfnl          = fullfile(pth,['r' nam ext]);

oNii         = nifti;
oNii.dat     = file_array(nfnl,dmn,NiiLabels(1).dat.dtype,NiiLabels(1).dat.offset,NiiLabels(1).dat.scl_slope,NiiLabels(1).dat.scl_inter);
oNii.mat     = Mn;
oNii.mat0    = Mn;
oNii.descrip = 'Labels';
create(oNii);
oNii.dat(:)  = l12(:);

delete(fnl);
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