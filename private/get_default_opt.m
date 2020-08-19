function opt = get_default_opt(opt)
% Do
if ~isfield(opt,'do'),             opt.do             = struct; end
if ~isfield(opt.do,'res_orig'),    opt.do.res_orig    = false;  end
if ~isfield(opt.do,'real_mni'),    opt.do.real_mni    = false;  end
if ~isfield(opt.do,'nm_reorient'), opt.do.nm_reorient = false;  end
if ~isfield(opt.do,'crop'),        opt.do.crop        = false;  end
if ~isfield(opt.do,'coreg'),       opt.do.coreg       = false;  end
if ~isfield(opt.do,'denoise'),     opt.do.denoise     = false;  end
if ~isfield(opt.do,'superres'),    opt.do.superres    = false;  end
if ~isfield(opt.do,'reslice'),     opt.do.reslice     = false;  end
if ~isfield(opt.do,'vx'),          opt.do.vx          = false;  end
if ~isfield(opt.do,'write2d'),     opt.do.write2d     = false;  end
if ~isfield(opt.do,'writemat'),    opt.do.writemat    = false;  end
if ~isfield(opt.do,'segment'),     opt.do.segment     = false;  end
if ~isfield(opt.do,'normalise'),   opt.do.normalise   = false;  end
if ~isfield(opt.do,'erode'),       opt.do.erode       = false;  end
if ~isfield(opt.do,'skullstrip'),  opt.do.skullstrip  = false;  end
if ~isfield(opt.do,'bfcorr'),      opt.do.bfcorr      = false;  end
if ~isfield(opt.do,'int_norm'),    opt.do.int_norm    = false;  end
if ~isfield(opt.do,'bb_spm'),      opt.do.bb_spm      = false;  end
if ~isfield(opt.do,'go2native'),   opt.do.go2native   = true;   end
% Output directory
if ~isfield(opt,'dir_out'),   opt.dir_out   = 'output';    end
if ~isfield(opt,'dir_out2d'), opt.dir_out2d = opt.dir_out; end
% int_norm options
if ~isfield(opt,'int_norm'),     opt.int_norm    = struct; end
if ~isfield(opt.int_norm,'rng'), opt.int_norm.rng = [0 1]; end
% Labels options
if ~isfield(opt,'labels'),      opt.labels      = struct; end
if ~isfield(opt.labels,'part'), opt.labels.part = [];     end
% reset-origin options
if ~isfield(opt,'reset_orig'),          opt.reset_orig          = struct; end
if ~isfield(opt.reset_orig,'only_neg'), opt.reset_orig.only_neg = false;  end
% realign2mni options
if ~isfield(opt,'realign2mni'),       opt.realign2mni       = struct; end
if ~isfield(opt.realign2mni,'rigid'), opt.realign2mni.rigid = true;   end
% coreg options
if ~isfield(opt,'coreg'),     opt.coreg     = struct; end
if ~isfield(opt.coreg,'ref'), opt.coreg.ref = 1;      end
% crop options
if ~isfield(opt,'crop'),           opt.crop           = struct; end
if ~isfield(opt.crop,'keep_neck'), opt.crop.keep_neck = true;   end
% superres options
if ~isfield(opt,'superres'),         opt.superres         = struct; end
if ~isfield(opt.superres,'ix'),      opt.superres.ix      = [];     end
if ~isfield(opt.superres,'Verbose'), opt.superres.Verbose = 0;      end
% Reslice options
if ~isfield(opt,'reslice'),     opt.reslice     = struct; end
if ~isfield(opt.reslice,'ref'), opt.reslice.ref = 1;      end
if ~isfield(opt.reslice,'deg'), opt.reslice.deg = 1;      end
% Voxel size options
if ~isfield(opt,'vx'),         opt.vx          = struct; end
if ~isfield(opt.vx,'size'),    opt.vx.size     = 1;      end
if ~isfield(opt.vx,'deg'),     opt.vx.deg      = 1;      end
% 2D options
if ~isfield(opt,'write2d'),         opt.write2d         = struct; end
if ~isfield(opt.write2d,'deg'),     opt.write2d.deg     = 0;      end
if ~isfield(opt.write2d,'sliceix'), opt.write2d.sliceix = [];     end
% Segment options
if ~isfield(opt,'segment'),          opt.segment          = struct; end
if ~isfield(opt.segment,'write_tc'), opt.segment.write_tc = [];     end
if ~isfield(opt.segment,'write_bf'), opt.segment.write_bf = [];     end
if ~isfield(opt.segment,'write_df'), opt.segment.write_df = [];     end
if ~isfield(opt.segment,'dir_out'),  opt.segment.dir_out  = [];     end
if ~isfield(opt.segment,'make_4d'),  opt.segment.make_4d  = false;  end
if ~isfield(opt.segment,'mask'),     opt.segment.mask     = false;  end
% Normalise options
if ~isfield(opt,'normalise'),      opt.normalise      = struct; end
if ~isfield(opt.normalise,'mask'), opt.normalise.mask = false;  end
if ~isfield(opt.normalise,'vol'),  opt.normalise.vol  = 1;      end
% Path to template (good for using, e.g., VoxelMorph)
if ~isfield(opt,'pth_template'),   opt.pth_template = []; end

if opt.do.bb_spm
    % Cropping to SPM12 atlas requires aligning input to MNI
    opt.do.real_mni = true;
end

if opt.do.segment && (opt.do.skullstrip || opt.do.bfcorr)
    warning('OBS: Not writing segmentations, because opt.do.skullstrip || opt.do.bfcorr'); 
end

opt.do.segment0 = opt.do.segment;
if opt.do.skullstrip || opt.do.bfcorr    
    opt.do.segment = true;
    if opt.do.bfcorr
        opt.segment.write_bf = [false true];
    end
    if opt.do.skullstrip
        opt.segment.write_tc        = false(6,4);
        opt.segment.write_tc(1:3,1) = true;
    end
end
%==========================================================================