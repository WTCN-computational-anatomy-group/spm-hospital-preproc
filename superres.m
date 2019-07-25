function Nii = superres(Nii,ix,Verbose)
if nargin < 2, ix      = []; end
if nargin < 3, Verbose = 1;  end

fprintf('Super-resolving...')

if ~isempty(ix)
    % Multiple observations of one, or multiple, channels. Make sure input
    % is correctly formated
    C   = numel(unique(ix));    
    tmp = cell(1,C);
    for c=1:C
        tmp{c} = Nii{1}(ix == c);
    end
    
    oNii = do_superres(tmp,Verbose);
    
    N = numel(Nii{1});
    for n=1:N         
        f = Nii{1}(n).dat.fname;
        delete(f);
    end
    
    Nii{1} = oNii;
else
    % One observation of each channels
    do_superres(Nii{1},Verbose);
    
    N = numel(Nii{1});
    for n=1:N         
        f             = Nii{1}(n).dat.fname;
        [pth,nam,ext] = fileparts(f);
        nf            = fullfile(pth,['srds' nam ext]);

        delete(f);
        Nii{1}(n) = nifti(nf);
    end
end

fprintf('done!\n')
%==========================================================================

%==========================================================================
function oNii = do_superres(Nii,Verbose)
CoRegister          = false;    
WorkersParfor       = Inf;
Method              = 'superres';
RegScaleSuperResMRI = 5;
ReadWrite           = false;
SliceGap            = 0;
DecreasingReg       = true;
IterMax             = 12;
IterImage           = 3;
ZeroMissingValues   = false;

if iscell(Nii)
    pth = fileparts(Nii{1}(1).dat.fname);
else
    pth = fileparts(Nii(1).dat.fname);
end
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
            
oNii = spm_mtv_preproc(fun_args{:});
%==========================================================================