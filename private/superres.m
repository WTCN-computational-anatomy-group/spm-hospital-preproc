function Nii = superres(Nii,CoRegDone,opt)

ix      = opt.ix;
Verbose = opt.Verbose;
DoCoreg = true;
if CoRegDone
    DoCoreg = false;
end

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
            P{c}{cnt} = Nii{1}(n).dat.fname;
            cnt       = cnt + 1;
        end
        P{c} = char(P{c}); % Needs to be character array not cell array
    end
    
    oNii = do_superres(P,Verbose,DoCoreg);
    
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
    
    do_superres(P,Verbose,DoCoreg);
    
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
function oNii = do_superres(P,Verbose,DoCoReg)
% Options
opt = struct;
opt.DoCoReg = DoCoReg;
% opt.IterMax = 1;  % Uncomment for testing
opt.Verbose = Verbose;
% Run
oNii = spm_superres(P,opt);
%==========================================================================
