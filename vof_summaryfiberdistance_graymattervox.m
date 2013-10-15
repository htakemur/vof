function vof_summaryfiberdistance_graymattervox(VofFiberToLoad, threshold, targetroiDir, targetroiNames, savematfile)

% Compute the proportion of gray matter voxels within a certain distance
% from VOF termination in each gray matter ROIs
%
% INPUT:
% VofFiberToLoad: a full path to fiber files (.pdb or .mat). Defined as strings to allow to load multiple fiber files.
% threshold: a threshold on the squared distance between gray matter voxel and fiber termination point. 
%            Code compute the proportion of gray matter voxel within these distances from VOF termination. 
%            Default is [2.25 9 20.25 36]. (1.5, 3, 4.5, 6 mm in ACPC space)
% targetroiDir: a directory containing gray matter ROI files
% targetroiNames: file names of gray matter ROIs (.mat format, mrDiffusion)
% savematfile: file name of output mat file
%
% EXAMPLE:
% VofFiberToLoad{1} = '/home/vof/data/S1/fibers/LVOF1.pdb';
% VofFiberToLoad{2} = '/home/vof/data/S1/fibers/LVOF2.pdb';
% VofFiberToLoad{3} = '/home/vof/data/S1/fibers/LVOF3.pdb';
% threshold = [2.25 9 20.25 36];
% targetroiDir = '/home/vof/data/S1/ROIs/map/';
% targetroiNames = {'LV3AB.mat','LV3d.mat','LV2d.mat','LV1.mat','LV2v.mat','LV3v.mat','LhV4.mat','LVO1.mat','LLO1.mat','LLO2.mat'};
% savematfile = 'S1LH_grayROIproportion_VOF.mat';
% vof_summaryfiberdistance_graymattervox(VofFiberToLoad, threshold, targetroiDir, targetroiNames, savematfile)
% 
% Copyright (C) Hiromasa Takemura, 2013, Stanford VISTA team


% Argument checking
if notDefined('threshold')
    threshold = [2.25 9 20.25 36];
end

sizevoffile = size(VofFiberToLoad);

for i =1:sizevoffile(2)
    
    % Load fibers
    vof = fgRead(VofFiberToLoad{i});
    voffiberlength(i) = length(vof.fibers);
    
    % Extract ACPC coordinate of dorsal and ventral VOF endpoint
    for kk = 1:voffiberlength(i)
        fibercoordinate = cell2mat(vof.fibers(kk));
        fiberlength = size(fibercoordinate);
        if  fibercoordinate(3,1) > fibercoordinate(3,fiberlength(2))
            dorsalfiberend(kk,1,i) = fibercoordinate(1,1);
            dorsalfiberend(kk,2,i) = fibercoordinate(2,1);
            dorsalfiberend(kk,3,i) = fibercoordinate(3,1);
            ventralfiberend(kk,1,i) = fibercoordinate(1,fiberlength(2));
            ventralfiberend(kk,2,i) = fibercoordinate(2,fiberlength(2));
            ventralfiberend(kk,3,i) = fibercoordinate(3,fiberlength(2));
            
        else
            ventralfiberend(kk,1,i) = fibercoordinate(1,1);
            ventralfiberend(kk,2,i) = fibercoordinate(2,1);
            ventralfiberend(kk,3,i) = fibercoordinate(3,1);
            dorsalfiberend(kk,1,i) = fibercoordinate(1,fiberlength(2));
            dorsalfiberend(kk,2,i) = fibercoordinate(2,fiberlength(2));
            dorsalfiberend(kk,3,i) = fibercoordinate(3,fiberlength(2));
            
        end
        clear fibercoordinate fiberlength
    end
    clear vof
    
    for kk = 1:length(targetroiNames)
        % Load target ROIs
        for troi = 1:length(targetroiNames)
            trois{troi} = fullfile(targetroiDir ,targetroiNames{troi});
        end
        
        % Extract ROI coordinates in ACPC space and the transpose those
        load(trois{kk});
        roicoords = transpose(roi.coords);
        dorsalendcoords = transpose(dorsalfiberend(:,:,i));
        ventralendcoords = transpose(ventralfiberend(:,:,i));
        roicoordsize = size(roicoords);
        
        % Compute the best square distance between particular ROI voxel and
        % any VOF fiber endpoint
        
        for kp = 1:roicoordsize(2)
            
            [indices_d(kp), bestSqDis_d(kp)] = nearpoints(roicoords(:,kp), dorsalendcoords(:,1:voffiberlength(i)));
            [indices_v(kp), bestSqDis_v(kp)] = nearpoints(roicoords(:,kp), ventralendcoords(:,1:voffiberlength(i)));
            
        end
        
        % Find a number of voxels which is within thresholded distance
        
        for ko = 1:length(threshold)
            dorsal_vofvoxelnum(kk, i, ko) = length(find((bestSqDis_d)<threshold(ko)));
            ventral_vofvoxelnum(kk, i, ko) = length(find((bestSqDis_v)<threshold(ko)));
            dorsal_vofvoxelnum_norm(kk, i, ko) = dorsal_vofvoxelnum(kk, i, ko)/length(bestSqDis_d);
            ventral_vofvoxelnum_norm(kk, i, ko) = ventral_vofvoxelnum(kk, i, ko)/length(bestSqDis_v);
            
        end
        clear roicoords dorsalendcoords ventralendcoords roicoordsize indices_d bestSqDis_d indices_v bestSqDis_v
    end
    
    
end
save(savematfile, 'dorsal_vofvoxelnum','ventral_vofvoxelnum','dorsal_vofvoxelnum_norm','ventral_vofvoxelnum_norm');
