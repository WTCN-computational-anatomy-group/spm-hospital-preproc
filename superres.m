function Nii = superres(Nii,Verbose)
if nargin < 2, Verbose = 0;     end

fprintf('Super-resolving...')

do_superres(Nii{1},Verbose);

N = numel(Nii{1});
for n=1:N         
    f             = Nii{1}(n).dat.fname;
    [pth,nam,ext] = fileparts(f);
    nf            = fullfile(pth,['sr' nam ext]);
    
    delete(f);
    Nii{1}(n) = nifti(nf);
end
fprintf('done!\n')
%==========================================================================

%==========================================================================
function do_superres(Nii,Verbose)
CoRegister           = false;    
WorkersParfor        = Inf;
Method               = 'superres';
RegScaleSuperResMRI  = 5;
ReadWrite            = false;
SliceGap             = 0;
DecreasingReg        = true;
IterMax              = 10;
IterImage            = 3;
    
fun_args = {'InputImages',Nii, ...
            'Method',Method, ...
            'Verbose',Verbose, ...
            'RegScaleSuperResMRI',RegScaleSuperResMRI, ...
            'ReadWrite',ReadWrite, ...
            'SliceGap',SliceGap, ...
            'CoRegister',CoRegister, ...
            'DecreasingReg',DecreasingReg, ...
            'IterMax',IterMax, ...
            'WorkersParfor',WorkersParfor, ...
            'IterImage',IterImage};
            
spm_mtv_preproc(fun_args{:});
%==========================================================================