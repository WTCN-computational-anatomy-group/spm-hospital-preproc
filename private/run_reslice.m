function V = run_reslice(fnames,deg)
matlabbatch{1}.spm.spatial.realign.write.data            = fnames;
matlabbatch{1}.spm.spatial.realign.write.roptions.which  = [1 0];
matlabbatch{1}.spm.spatial.realign.write.roptions.interp = deg;
matlabbatch{1}.spm.spatial.realign.write.roptions.wrap   = [0 0 0];
matlabbatch{1}.spm.spatial.realign.write.roptions.mask   = 0;
matlabbatch{1}.spm.spatial.realign.write.roptions.prefix = 'r';
spm_jobman('run',matlabbatch);

V = spm_vol(fnames{1});
for n=2:numel(fnames)
    f             = fnames{n};
    [pth,nam,ext] = fileparts(f);
    nf            = fullfile(pth,['r' nam ext]);
    V(n)          = spm_vol(nf);
end
%==========================================================================