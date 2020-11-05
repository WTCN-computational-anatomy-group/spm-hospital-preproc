function [Nii,was_gz,nams] = read_and_copy(Nii,DirOut,prefix)
fprintf('Making copies...')

if ischar(Nii)
    Nii = {Nii};
end
if ischar(Nii{1})
    Nii = {Nii};
end

was_gz = false;
if ~isa(Nii{1}, 'nifti')
    % Input are paths
    N = numel(Nii{1});    
    for n=1:N
        f = Nii{1}{n};
        [dir,~,ext] = fileparts(f);
        if strcmp(ext, '.gz')
            % Input is zipped -> extract
            f = gunzip(f);
            [~,nam,ext] = fileparts(f{1});
            Nii{1}{n} = fullfile(dir, [nam ext]);
            was_gz = true;
        end
    end
    Nii{1} = nifti(Nii{1});
end

N = numel(Nii{1});
nams = cell(2, N);
for n=1:N      
    Nii_n = Nii{1}(n);        
    f     = Nii_n.dat.fname;
    
    [~,nam]    = fileparts(f);
    nam        = [prefix nam];  % Integrate prefix option
    nams{1, n} = nam;
    
    nf = fullfile(DirOut,[nam '.nii']);
    if exist(nf, 'file') == 2
        delete(nf); 
    end   
            
    create_nii(nf,Nii_n.dat(),Nii_n.mat,Nii_n.dat.dtype,'copy',...
               Nii_n.dat.offset,Nii_n.dat.scl_slope,Nii_n.dat.scl_inter);
                   
    Nii{1}(n) = nifti(nf);
   
    if was_gz
        delete(f)
    end
end

if numel(Nii) > 1
    % Deal with labels
            
    if ~isa(Nii{2}, 'nifti')
        % Input are paths
        N = numel(Nii{2});    
        for n=1:N
            f = Nii{2}{n};
            [dir,~,ext] = fileparts(f);
            if strcmp(ext, '.gz')
                % Input is zipped -> extract
                f = gunzip(f);
                [~,nam,ext] = fileparts(f{1});
                Nii{2}{n} = fullfile(dir, [nam ext]);
            end
        end
        Nii{2} = nifti(Nii{2});
    end

    for n=1:N
        if isempty(Nii{2}(n).dat), continue; end
        
        Nii_n   = Nii{2}(n); 
        f       = Nii_n.dat.fname;
        
        [~,nam]    = fileparts(f);
        nam        = [prefix nam];  % Integrate prefix option
        nams{2, n} = nam;
        
        nf = fullfile(DirOut,[nam '.nii']);
        if exist(nf, 'file') == 2
            delete(nf); 
        end

        im = Nii_n.dat();        
        if max(im(:)) > 255, error('Copy labels: values over 255!'); end
        
        create_nii(nf,im,Nii_n.mat,spm_type('uint8'),'labels');        
        
        Nii{2}(n) = nifti(nf);
        
        if was_gz
            delete(f)
        end
    end    
end
fprintf('done!\n')
%==========================================================================