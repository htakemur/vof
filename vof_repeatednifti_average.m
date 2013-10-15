function vof_repeatednifti_average(fibernii, avgniiname)

% Averaging multiple nifti files.
% 
% INPUT:
% fibernii = full path to nifti files (string)
% avgniiname: output nifti file name

% EXAMPLE:
% fibernii{1} = 'LVOF1_rmsexv_maxDist3_ventralfiberendpoint.nii.gz';
% fibernii{2} = 'LVOF2_rmsexv_maxDist3_ventralfiberendpoint.nii.gz';
% fibernii{3} = 'LVOF3_rmsexv_maxDist3_ventralfiberendpoint.nii.gz';
% avgniiname = 'LVOFavg_rmsexv_maxDist3_ventralfiberendpoint.nii.gz';
% vof_repeatednifti_average(fibernii, avgniiname)
%
% (C) Hiromasa Takemura, 2013, Stanford VISTA Team

%% Read nifti files

numniftifile = size(fibernii);

for k = 1:numniftifile(2)
fiberniiload(k) = niftiRead(fibernii{k});
end

%% Averaging nifti files
averagenii = niftiRead(fibernii{1});
for i = 2:numniftifile(2)
averagenii.data = averagenii.data + fiberniiload(i).data;
end
averagenii.data = averagenii.data/numniftifile(2);
averagenii.fname = avgniiname;

%% Save
niftiWrite(averagenii, avgniiname);

