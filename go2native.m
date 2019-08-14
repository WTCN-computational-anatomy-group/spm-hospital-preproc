function go2native(out)
C = numel(out.pth.im);
for c=1:C    
    Mc = spm_get_space(out.pth.im{c}); 
    spm_get_space(out.pth.im{c},out.mat{c}*Mc); 

    if ~isempty(out.pth.lab{c})
        Mc = spm_get_space(out.pth.lab{c}); 
        spm_get_space(out.pth.lab{c},out.mat{c}*Mc); 
    end        
end

if ~isempty(out.pth.seg)
    for i1=1:numel(out.pth.seg{1})
        fname = deblank(out.pth.seg{1}{i1});
        Mc    = spm_get_space(fname);
        spm_get_space(fname,out.mat{1}*Mc);
    end
end
%==========================================================================