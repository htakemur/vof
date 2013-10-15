function  vof_fiberendpointnifti_dorsalventral(fe, keepFG, t1file, fname, skernal, option)

% Extract the fiber endpoint of VOF fibers and create the nifti file
% describing the endpoint density.
%
% INPUT:
% fe: fe structure (created by LIFE code, which is available at https://github.com/francopestilli/life).
% keepFG: one dimensional matrix describing which fibers in fe sturcture are included in the analysis. 1: included, 0: excluded.
% t1file: full path to nifti file of t1-weighted image data
% fname: file name for the output nifti file
% skernal: size smoothing Gaussian kernel for fiber endpoint density. Default is [3 3 3], which means 3 voxels in the resolution of diffusion MRI data
% option: 1, computing dorsal VOF endpoint density. 2, computing ventral VOF endpoint density.
% 
% Hiromasa Takemura, (c) Stanford Vista Team, 2013

%% Argument Checking
if notDefined('skernal')
    skernal= [3 3 3];
end

%% Read fibers in image coordinate
fefiber = feGet(fe,'fibers img');
fgkeep = fgExtract(fefiber, transpose(keepFG), 'keep');
t1 = niftiRead(t1file);

%% Extract fiber endpoints (option=1, dorsal VOF endpoint; option=2, ventral VOF endpoint)
fbsize = size(fgkeep.fibers);
for ph = 1:fbsize(1)
    f_coords = cell2mat(fgkeep.fibers(ph));
    fcoordsize = size(f_coords);
    switch option
        case 1,
            if f_coords(3,1)>f_coords(3,fcoordsize(2))
                fendpoints(1, ph) = round(f_coords(1,1));
                fendpoints(2, ph) = round(f_coords(2,1));
                fendpoints(3, ph) = round(f_coords(3,1));
            else
                fendpoints(1, ph) = round(f_coords(1,fcoordsize(2)));
                fendpoints(2, ph) = round(f_coords(2,fcoordsize(2)));
                fendpoints(3, ph) = round(f_coords(3,fcoordsize(2)));
            end
            
        case 2,
            if f_coords(3,1)>f_coords(3,fcoordsize(2))
                fendpoints(1, ph) = round(f_coords(1,fcoordsize(2)));
                fendpoints(2, ph) = round(f_coords(2,fcoordsize(2)));
                fendpoints(3, ph) = round(f_coords(3,fcoordsize(2)));
                
            else
                fendpoints(1, ph) = round(f_coords(1,1));
                fendpoints(2, ph) = round(f_coords(2,1));
                fendpoints(3, ph) = round(f_coords(3,1));
            end
    end
end

%% Attaching endpoint density data into nifti file
dwi = niftiRead(feGet(fe,'dwifile'));
nii = dwi;
nii.fname = fname;
nii.data = zeros(dwi.dim(1), dwi.dim(2), dwi.dim(3));
nii.ndim = numel(size(nii.data));
nii.dim  = nii.dim(1:nii.ndim);
nii.pixdim = [dwi.pixdim(1:nii.ndim)];
nii.descrip = 'Fiber density';
nii.data_type = niftiClass2DataType(class(nii.data));
nii.nifti_type = 2;

for pp = 1:fbsize(1)
        nii.data(fendpoints(1, pp), fendpoints(2, pp), fendpoints(3, pp)) =  nii.data(fendpoints(1, pp), fendpoints(2, pp), fendpoints(3, pp)) + 1;
end
%% Spatial Smoothing
nii.data = smooth3(nii.data,'gaussian', skernal);

%% Resample to T1 image size, then save file
[finalResNifti, resample_params] = mrAnatResampleToNifti(nii, t1);
niftiWrite(finalResNifti);
