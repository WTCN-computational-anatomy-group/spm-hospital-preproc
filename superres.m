function Nii = superres(Nii,Verbose)
if nargin < 2, Verbose = 0;     end

fprintf('Super-resolving...')

do_superres(Nii{1},Verbose);

N = numel(Nii{1});
for n=1:N         
    f             = Nii{1}(n).dat.fname;
    [pth,nam,ext] = fileparts(f);
    nf            = fullfile(pth,['srds' nam ext]);
    
    delete(f);
    Nii{1}(n) = nifti(nf);
end
fprintf('done!\n')
%==========================================================================

%==========================================================================
function do_superres(Nii,Verbose)
CoRegister          = false;    
WorkersParfor       = Inf;
Method              = 'superres';
RegScaleSuperResMRI = 4;
ReadWrite           = false;
SliceGap            = 0;
DecreasingReg       = true;
IterMax             = 12;
IterImage           = 3;
ZeroMissingValues   = false;

pth             = fileparts(Nii(1).dat.fname);
OutputDirectory = pth;

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
            'OutputDirectory',OutputDirectory, ...
            'ZeroMissingValues',ZeroMissingValues, ...
            'IterImage',IterImage};
            
spm_mtv_preproc(fun_args{:});
%==========================================================================
