function Nii = superres(Nii,opt)

ix      = opt.ix;
Verbose = opt.Verbose;

fprintf('Super-resolving...')

if ~isempty(ix)
    % Multiple observations of one, or multiple, channels. Make sure input
    % is correctly formated
    C = numel(unique(ix));    
    P = cell(1,C);
    for c=1:C        
        N    = sum(ix == c);
        P{c} = cell(1,N);
        cnt  = 1;
        for n=find(ix == c)
            P{c}{cnt} = Nii{1}(n);
            cnt       = cnt + 1;
        end
    end
    
    oNii = do_superres(P,Verbose);
    
    N = numel(Nii{1});
    for n=1:N         
        f = Nii{1}(n).dat.fname;
        delete(f);
    end
    
    Nii{1} = oNii;
else
    % One observation of each channels
    C = numel(Nii{1});
    P = cell(1,C);
    for c=1:C
        P{c} = Nii{1}(c).dat.fname;
    end
    
    do_superres(P,Verbose);
    
    N = numel(Nii{1});
    for n=1:N         
        f             = Nii{1}(n).dat.fname;
        [pth,nam,ext] = fileparts(f);
        nf            = fullfile(pth,['y' nam ext]);

        delete(f);
        Nii{1}(n) = nifti(nf);
    end
end

fprintf('done!\n')
%==========================================================================

%==========================================================================
function oNii = do_superres(P,Verbose)
% CoRegister          = false;    
% WorkersParfor       = Inf;
% Method              = 'superres';
% RegScaleSuperResMRI = 4;
% ReadWrite           = false;
% SliceGap            = 0;
% DecreasingReg       = true;
% IterMax             = 12;
% IterImage           = 3;
% ZeroMissingValues   = false;
% 
% if iscell(Nii)
%     pth = fileparts(Nii{1}(1).dat.fname);
% else
%     pth = fileparts(Nii(1).dat.fname);
% end
% OutputDirectory = pth;
% 
% fun_args = {'InputImages',Nii, ...
%             'Method',Method, ...
%             'Verbose',Verbose, ...
%             'RegScaleSuperResMRI',RegScaleSuperResMRI, ...
%             'ReadWrite',ReadWrite, ...
%             'SliceGap',SliceGap, ...
%             'CoRegister',CoRegister, ...
%             'DecreasingReg',DecreasingReg, ...
%             'IterMax',IterMax, ...
%             'WorkersParfor',WorkersParfor, ...
%             'OutputDirectory',OutputDirectory, ...
%             'ZeroMissingValues',ZeroMissingValues, ...
%             'IterImage',IterImage};
            
opt.DoCoReg = false;

oNii = spm_superres(P,opt);
%==========================================================================