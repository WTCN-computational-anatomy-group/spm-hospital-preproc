function Nii = read_data(Data,DirOut)
if isa(Data,'nifti')
    % Is nifti object
    Nii = {Data};
elseif ischar(Data)
    % Is a char array
    [~,~,ext] = fileparts(Data(1,:));
    if strcmpi(ext,'.nii')
        Nii = {nifti(Data)};
    elseif strcmpi(ext,'.dcm')
        Nii = dcm2nii(Data,DirOut);
    end
else
    error('Input error!')
end
%==========================================================================

%==========================================================================
function Nii = dcm2nii(Data,DirOut)
hdr = spm_dicom_headers(Data);
out = spm_dicom_convert(hdr,'all','flat',spm_get_defaults('images.format'),DirOut);
Nii = {nifti(out.files{1})};
%==========================================================================