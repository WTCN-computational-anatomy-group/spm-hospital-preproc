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
            
    write_nii(nf, Nii_n.dat(), Nii_n.mat, 'Patient-Preprocessing Copy', 'int16');
                   
    Nii{1}(n) = nifti(nf);
   
    if was_gz
        delete(f)
    end
end

if numel(Nii) > 1
    % Deal with labels
            
    if ~isa(Nii{2}, 'nifti')
        % Input are paths
        N1   = numel(Nii{1});
        Nii2 = nifti;
        for n=1:N1
            if n > numel(Nii{2}) || isempty(Nii{2}{n}), continue; end
            
            f = Nii{2}{n};
            [dir,~,ext] = fileparts(f);
            if strcmp(ext, '.gz')
                % Input is zipped -> extract
                f = gunzip(f);
                [~,nam,ext] = fileparts(f{1});
                f           = fullfile(dir, [nam ext]);
            end
            Nii2(n) = nifti(f);
        end
        Nii{2} = Nii2;
        clear Nii2
    end

    N = numel(Nii{2});
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
        
        write_nii(nf, im, Nii_n.mat, 'Patient-Preprocessing Copy', 'uint8');            
        
        Nii{2}(n) = nifti(nf);
        
        if was_gz
            delete(f)
        end
    end    
end
fprintf('done!\n')
%==========================================================================