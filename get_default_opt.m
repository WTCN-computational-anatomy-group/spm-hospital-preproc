function opt = get_default_opt(opt)
% Do
if ~isfield(opt,'do'),          opt.do          = struct; end
if ~isfield(opt.do,'res_orig'), opt.do.res_orig = true;   end
if ~isfield(opt.do,'real_mni'), opt.do.real_mni = true;   end
if ~isfield(opt.do,'crop'),     opt.do.crop     = true;   end
if ~isfield(opt.do,'coreg'),    opt.do.coreg    = true;   end
if ~isfield(opt.do,'denoise'),  opt.do.denoise  = true;   end
if ~isfield(opt.do,'reslice'),  opt.do.reslice  = true;   end
if ~isfield(opt.do,'vx'),       opt.do.vx       = true;   end
if ~isfield(opt.do,'write2d'),  opt.do.write2d  = true;   end
if ~isfield(opt.do,'writemat'), opt.do.writemat = false;   end
% Output directory
if ~isfield(opt,'dir_out'),     opt.dir_out   = 'output';    end
if ~isfield(opt,'dir_out2d'),   opt.dir_out2d = opt.dir_out; end
% coreg options
if ~isfield(opt,'coreg'),       opt.coreg     = struct; end
if ~isfield(opt.coreg,'ref'),   opt.coreg.ref = 1;      end
% Reslice options
if ~isfield(opt,'reslice'),     opt.reslice     = struct; end
if ~isfield(opt.reslice,'ref'), opt.reslice.ref = 1;      end
% Voxel size options
if ~isfield(opt,'vx'),          opt.vx      = struct; end
if ~isfield(opt.vx,'size'),     opt.vx.size = 1;      end
if ~isfield(opt.vx,'def'),      opt.vx.deg  = 4;      end
% Crop options
if ~isfield(opt,'crop'),        opt.crop      = struct; end
if ~isfield(opt.crop,'neck'),   opt.crop.neck = false;  end
% 2D options
if ~isfield(opt,'write2d'),         opt.write2d         = struct; end
if ~isfield(opt.write2d,'deg'),     opt.write2d.deg     = 0;      end
if ~isfield(opt.write2d,'axis_2d'), opt.write2d.axis_2d = 3;      end
if ~isfield(opt.write2d,'sliceix'), opt.write2d.sliceix = [];     end
% Path to denoising toolbox (https://github.com/WCHN/mtv-preproc)
if ~isfield(opt,'pth_mtv'),     opt.pth_mtv   = '/home/mbrud/dev/mbrud/code/matlab/MTV-preproc'; end
%==========================================================================