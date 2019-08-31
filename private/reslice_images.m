function [Nii,M] = reslice_images(Nii,M,opt)

ref_ix = opt.ref;
deg    = opt.deg;

N = numel(Nii{1});
if N > 1    
    fprintf('Reslicing images...')

    ixs       = 1:N;
    source_ix = ixs(ixs ~= ref_ix);

    % Use SPM batch job to reslice
    fnames = {Nii{1}(ref_ix).dat.fname, Nii{1}(source_ix).dat.fname}';
    run_reslice(fnames,deg);

    % Update M
    for n=source_ix
        M{n} = M{ref_ix};
    end

    % Update NIIs
    for n=source_ix
        f             = Nii{1}(n).dat.fname;
        [pth,nam,ext] = fileparts(f);
        nf            = fullfile(pth,['r' nam ext]);

        delete(f);
        Nii{1}(n) = nifti(nf);
    end
    
    fprintf('done!\n')
end
%==========================================================================