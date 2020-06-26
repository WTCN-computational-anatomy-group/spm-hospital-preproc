# Patient-Preprocessing

This is MATLAB code for various preprocessing operations (registration, reslicing, denoising, segmentation, etc.) of neuroimaging data. It takes as input nifti files and produces copies of this data to which the requested preprocessing steps are applied. It additionally handles image data with pairied label masks (e.g., a T1w MRI and a tumour mask), and makes sure that the resulting preprocessed data is consistent. See below for example use cases.

## Dependencies

The algorithm requires that the following package is on the MATLAB path:
* **SPM12:** Download from https://www.fil.ion.ucl.ac.uk/spm/software/spm12/.

## Use cases

### Example 1

Below is a MATLAB snippet that takes as input MR images of multiple sequences and produces images that have been: co-registered, resliced, and made to have 1 mm isotropic voxel size. This code could be run to, for example, produce inputs to a multi-channel segmentation routine.
```
% Paths to images
pth_img = {'MRI_T1w.nii', 'MRI_T2w.nii', 'MRI_PDw.nii'};  % Paths to image data nifti files as cell array

% Format input
data = nifti(pth_img);

% Set preprocessing options
opt             = struct;    
opt.dir_out     = 'output'; % Output directory
opt.do.coreg    = true;     % Co-register using NMI
opt.do.reslice  = true;     % Reslice to have same image grids
opt.reslice.ref = 1;        % Reslice to image pth_img(1)
opt.do.vx       = true;     % Change voxel size (1 mm isotropic by default)  
    
% Run preprocessing
RunPreproc(data, opt);
```

### Example 2

Below is a MATLAB snippet that takes as input an image and a label mask (both as niftis) and produces images that have been: rigidly realigned (to MNI space), cropped of neck and air data, made to have 2 mm isotropic voxel size, and made to have the same field-of-view as the SPM12 atlas. This code could be run on, for example, multiple subjects' images to produce input to some machine learning model.
```
% Format input
data    = cell(1, 2);
data{1} = nifti('img.nii');    % Give image nifti here
data{2} = nifti('labels.nii'); % Give label nifti here

% Set preprocessing options
opt                = struct;    
opt.dir_out        = 'output'; % Output directory
opt.do.real_mni    = true;     % Rigid alignment to MNI
opt.do.crop        = true;     % Remove air data (makes for smaller image)
opt.crop.keep_neck = false;    % Remove also neck (makes for even smaller image)
opt.do.vx          = true;     % Change voxel size
opt.vx.size        = 2;        % What voxel size do you want?
opt.do.bb_spm      = true;     % Make image have same bounding-box as default SPM12 template
  
% Run preprocessing
RunPreproc(data, opt);
```

## License

This software is released under the [GNU General Public License version 3](LICENSE) (GPL v3). As a result, you may copy, distribute and modify the software as long as you track changes/dates in source files. Any modifications to or software including (via compiler) GPL-licensed code must also be made available under the GPL along with build & install instructions.

[TL;DR: GPL v3](https://tldrlegal.com/license/gnu-general-public-license-v3-(gpl-3))
