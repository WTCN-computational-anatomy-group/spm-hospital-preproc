# Patient-Preprocessing

This is MATLAB code for various neuroimaging preprocessing operations (registration, reslicing, denoising, segmentation, etc.), which was originally intended for processing routine clinical data (hence the name) [1]. It takes as input nifti files (as .nii or .nii.gz) and produces copies of this data to which the requested preprocessing steps are applied. It additionally handles image data with paired label masks (e.g., a T1w MRI and a tumour mask (or multiple classes)), and makes sure that the resulting preprocessed data is consistent. See below for some example use cases, which could be run stand-alone or be inspiration for more complicated preprocessing tasks.

Please refer to the [wiki](home) for more details.

## Dependencies

The algorithm requires that the following package is on the MATLAB path:
* **SPM12:** Download from https://www.fil.ion.ucl.ac.uk/spm/software/spm12.
* **spm_superres:** Download/clone from https://github.com/brudfors/spm_superres (if you want to use the denoising or super-resolution options).

## Example use cases

### 1. Multi-channel MRI segmentation

This MATLAB snippet takes as input MR images of multiple sequences and produces images that have been co-registered and resliced. These images are then segmented using the SPM12 unified segmentation routine and native+template (unmodulated) space GM, WM and CSF segmentations are written to disk.
```
% Paths to multi-channel images
paths = {'MRI_T1w.nii', 'MRI_T2w.nii', 'MRI_PDw.nii'};

% Set preprocessing options
opt             = struct;    
opt.dir_out     = 'output'; % Output directory
opt.do.coreg    = true;     % Co-register using NMI
opt.do.reslice  = true;     % Reslice to have same image grids
opt.reslice.ref = 1;        % Reslice to image pth_img(1)
opt.do.segment  = true;     % Enable unified segmentation
% Write GM, WM and CSF segmentations (in native and unmodulated template space)
opt.segment.write_tc            = false(6,4);  
opt.segment.write_tc(1:3,[1 3]) = true;

% Run preprocessing
RunPreproc(paths, opt);
```

### 2. Image with label mask

This MATLAB snippet takes as input an image and a label mask (both as niftis) and produces images that have been: rigidly realigned (to MNI space), cropped of neck and air data, made to have 2 mm isotropic voxel size, and made to have the same field-of-view as the SPM12 atlas. This code could be run on, for example, multiple subjects' images to produce input to some machine learning model.
```
% Path to image and labels
paths = {{'img.nii'}, {'labels.nii'}};

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
RunPreproc(paths, opt);
```

### 3. MRI denoising (requires spm_superres)

This MATLAB snippet takes as input an MR image and applies a total variation denoising routine to it [2].
```
% Path to noisy MRI
paths = 'MRI.nii';

% Set preprocessing options
opt            = struct;    
opt.dir_out    = 'output'; % Output directory
opt.do.denoise = true;     % Enabla denoising
    
% Run preprocessing
RunPreproc(paths, opt);
```

### 4. Multi-channel MRI super-resolution (requires spm_superres)

This MATLAB snippet takes as input thick-sliced, multi-channel MR images and applies a super-resolution routine to it [2]; producing 1 mm isotropic images on the same grid.
```
% Paths to thick-sliced MRIs
paths = {'MRI_T1w.nii', 'MRI_T2w.nii', 'MRI_PDw.nii'};

% Set preprocessing options
opt             = struct;    
opt.dir_out     = 'output'; % Output directory
opt.do.superres = true;     % Enable super-resolution
    
% Run preprocessing
RunPreproc(paths, opt);
```

### 5. Changing the voxel size of a bunch of images

This MATLAB snippet simply changes the voxel size of multiple image volumes (e.g., MRIs or tissue segmentations), assumed to be located in the same folder. The interpolation order can be changed, here nearest neigbour is used.
```
% Path to images
dir_data = '/path/to/folder/with/images'; % Directory with images
ext      = '.nii';                        % File extension of images
prefix    = 'c';                          % Prefix of image filenames

% Get paths and convert to cell array
paths = spm_select('FPList',dir_data,['^' prefix '.*\' ext '$']);
paths = cellstr(paths);

% Set preprocessing options
opt         = struct;    
opt.dir_out = 'output'; % Output directory
opt.do.vx   = true;     % Enables changing the voxel size
opt.vx.size = 3;        % What output voxel size?
opt.vx.deg  = 0;        % What interpolation order (0: nearest neighbour)?

% Run preprocessing
RunPreproc(paths, opt);
```

## References

1. Brudfors, M. (2020). 
Generative Models for Preprocessing of Hospital Brain Scans.
Doctoral dissertation, UCL (University College London).

2. Brudfors M, Balbastre Y, Nachev P, Ashburner J.
MRI Super-Resolution Using Multi-channel Total Variation.
In Annual Conference on Medical Image Understanding and Analysis
2018 Jul 9 (pp. 217-228). Springer, Cham.    

## License

This software is released under the [GNU General Public License version 3](LICENSE) (GPL v3). As a result, you may copy, distribute and modify the software as long as you track changes/dates in source files. Any modifications to or software including (via compiler) GPL-licensed code must also be made available under the GPL along with build & install instructions.

[TL;DR: GPL v3](https://tldrlegal.com/license/gnu-general-public-license-v3-(gpl-3))
