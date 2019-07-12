clear;

DirData      = '/data/mbrud/populations/original/RIRE-T1T2CT';
DirPreproc   = '/data/mbrud/populations/preproc/haem/RIRE-T1T2CT';
DirPreproc2d = '/data/mbrud/populations/2d/haem/RIRE-T1T2CT';

if exist(DirPreproc,'dir'), rmdir(DirPreproc,'s'); end; mkdir(DirPreproc);   
if exist(DirPreproc2d,'dir'), rmdir(DirPreproc2d,'s'); end; mkdir(DirPreproc2d);   

files = spm_select('List',DirData,'^.*\.nii$');
N     = size(files,1);

% Get unique patients
Nii  = {};
s    = 1;
c    = 1;
onam = '';
ix   = 11;
for n=1:N
    FileName = deblank(files(n,:));
    FileName = fullfile(DirData,FileName);
    [~,nam]  = fileparts(FileName);
    nam      = nam(1:ix);
    
    if isempty(onam) || strcmp(nam,onam)
        if c == 1
            Nii{s}    = nifti(FileName);
        else        
            Nii{s}(c) = nifti(FileName);
        end
        
        c = c + 1;
    else        
        c = 2;
        s = s + 1;        
        
        Nii{s} = nifti(FileName);                        
    end
    
    onam = nam;
end

%% Preproc

S0      = numel(Nii);
sliceix = [60 48 40 23 55 40 55];
out     = cell(1,S0);
parfor s=1:S0
   
    opt             = struct;
    opt.do.real_mni = true;
    opt.do.coreg    = true;
    opt.do.denoise  = false;
    opt.do.crop     = true;
    opt.pth_mtv     = fullfile('/home','mbrud','dev','mbrud','code','matlab','MTV-preproc');
    opt.do.reslice  = true;
    opt.do.vx       = true;
    opt.dir_out     = DirPreproc;
    opt.reslice.ref = 2;
    opt.coreg.ref   = 2;    
    opt.dir_out2d   = DirPreproc2d;
    opt.do.writemat = false;
    opt.do.res_orig = true;
    
    opt.do.write2d      = true;
    opt.write2d.sliceix = sliceix(s);
    
    out{s} = RunPreproc(Nii{s},opt);
end

%% Make JSON
Channel = {'CT','T1','T2'};

for s=1:S0
    p = out{s}.pth.im2d;
    
    for c=1:numel(p)
        FileName    = p{c};
        [~,nam,ext] = fileparts(FileName);                
        
        j            = struct;
        j.name       = num2str(s);
        j.population = 'RIRE-T1T2CT';
        j.pth        = [nam '.nii'];                                     
        j.modality   = 'MRI';
        j.channel    = Channel{c};
        j.healthy    = true;         

        j = orderfields(j);

        pth_json = fullfile(DirPreproc2d,[nam '.json']);
        spm_jsonwrite(pth_json,j);
    end
end

%% test
dat = spm_json_manager('init_dat',DirPreproc2d);