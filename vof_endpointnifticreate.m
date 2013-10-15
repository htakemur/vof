function vof_endpointnifticreate(feFileToLoad, fibername, roisDir, roiNames, roiOperations, t1file, skernal)

% This code segments the VOF from the connectome (stored in fe structure, LIFE) and then create the nifti file containing the dorsal or ventral fiber endpoint density data.
% This code requires LIFE and MBA codes by Franco Pestilli, which are availabe at:
% https://github.com/francopestilli/life
% https://github.com/francopestilli/mba  
%
% INPUT:
% feFileToLoad: a full path to .mat file containing fe structure (generated by LIFE code). Use strings for loading multiple connectome data.
% fibername : a name for fibers. Used for file names of output nifti file.
% roisDir: a directory containing waypoint ROI mat files
% roiNames: names for waypoint ROIs (string)
% roiOperations: default is {'and', 'and'} (segmenting fibers passing through both waypoint ROIs)
% t1file : a full path to t1 nifti file
% skernal:  size smoothing Gaussian kernel for fiber endpoint density. Default is [3 3 3], which means 3 voxels in the resolution of diffusion MRI data
%
% EXAMPLE:
% feFileToLoad{1} = '/home/vof/data/S1/life/LH_1stconnectome_fe.mat';
% feFileToLoad{2} = '/home/vof/data/S1/life/LH_2ndconnectome_fe.mat';
% feFileToLoad{3} = '/home/vof/data/S1/life/LH_3rdconnectome_fe.mat';
% fibername = 'LVOF';
% roisDir = '/home/vof/data/S1/ROIs/Waypoint/'
% roiNames = {'LH_VOFPlane1.mat','LH_VOFPlane2.mat'};
% roiOperations = {'and','and'};
% t1file = '/home/vof/data/S1/t1/t1.nii.gz';
% skernal = [3 3 3];
% vof_endpointnifticreate(feFileToLoad, fibername, roisDir, roiNames, roiOperations, t1file, skernal)

% Copyright (C) Hiromasa Takemura, 2013, Stanford VISTA Team


% Argument checking
if notDefined('roiOperations')
    roiOperations = {'and','and'};
end
if notDefined('skernal')
    skernal= [3 3 3];
end

sizecon = size(feFileToLoad);

% Repeat the process if the input were multiple connectome files
for connum = 1:sizecon(2)
    
    % Define file names
    savedorsalfile = [fibername, num2str(connum),'_dorsalendpoint.nii.gz'] ;
    saveventralfile = [fibername, num2str(connum),'_ventralendpoint.nii.gz'] ;
    
    % Load ROIs
    disp('loading waypoint ROI files...')
    for iroi = 1:length(roiNames)
        rois{iroi} = fullfile(roisDir,roiNames{iroi});
    end
    
    % Load fe structure
    disp('loading the LiFE structure...')
    if ischar(feFileToLoad{connum})
        fprintf('Loading %s ...\n',feFileToLoad{connum})
        load(feFileToLoad{connum});
    else
        fe  =feFileToLoad{connum};
        clear feFileToLoad;
    end
    
    disp('Segmenting VOF from connectome...')
    % Extract the fiber group from the FE structure
    fg = feGet(fe,'fibers acpc');
    
    [fgsegment, keepFascicles] = feSegmentFascicleFromConnectome(fg, rois, roiOperations, 'prob connectome');
    
    disp('Excluding outliers from VOF...')
    
    % Exclude outliers
    [fgsegment2 keepFascicles2] = vof_gradient_removeoutlier(fe, keepFascicles,2,rois{1},rois{2});
    [fgsegment3, keepFascicle3] = mbaComputeFibersOutliers(fgsegment2,3,3,25);
    
    keepFascicles4 = keepFascicles2;
    keepFascicle2matrix = find(keepFascicles2);
    
    for kj = 1:length(keepFascicle3)
        if keepFascicle3(kj)==0,
            keepFascicles4(keepFascicle2matrix(kj)) = 0;
        else
        end
    end
    
    % Compute endpoint density, and then save it in nifti format
    disp('Computing fiber endpoint density...')
    vof_fiberendpointnifti_dorsalventral(fe, keepFascicles4, t1file, savedorsalfile, skernal, 1);
    vof_fiberendpointnifti_dorsalventral(fe, keepFascicles4, t1file, saveventralfile, skernal, 2);
    
end
