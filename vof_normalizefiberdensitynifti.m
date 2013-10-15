function vof_normalizefiberdensitynifti(niifile,saveFileName)

% Read nifti file of fiber density data and then normalize.
% The fiber endpoint density value of each voxel is divided by maximum value.
% 
% INPUT: 
% niifile: a full path to nifti file including fiber endpoint density data
% saveFileName: a file name for output nifti file including normalized
% fiber endpoint density
%
% (C) Hiromasa Takemura, 2013, Stanford VISTA team

%% Read nifti

nii = niftiRead(niifile);
niisave = nii;

%% Find out maximum density, and then normalize

maxnum = max(nii.data);
maxmax = max(maxnum);
maxmaxmax = max(maxmax);
niisave.data = nii.data/maxmaxmax;

%% Save file
niftiWrite(niisave, saveFileName);