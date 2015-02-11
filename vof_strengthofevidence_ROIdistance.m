function vof_strengthofevidence_ROIdistance(feFileToLoad,fibername,roisDir,roiNames,roiOperations,targetroidir,targetroiNames,threshold,nboots, nmontecarlo)

% Compute the Strength of Connection Evidence (S) for the VOF projection to
% particular ROIs. VOF projection is associated with fibers located within
% particular distance (threshold) from gray matter ROIs. We compared RMSE
% in model WITH and WITHOUT subset of VOF fibers within a certain distance
% from each ROIs.
%
% INPUT:
% feFileToLoad: a full path to .mat file containing fe structure (generated by LIFE code). Use strings for loading multiple connectome data.
% fibername: name for fiber groups. Used for file name of fiber file.
% roisDir: a directory containing waypoint ROI mat files
% roiNames: names for waypoint ROIs (string)
% roiOperations: default is {'and', 'and'} (segmenting fibers passing through both waypoint ROIs)
% targetroidir: a full path to directory containing target gray matter ROI
%               .mat files (e.g. 'LhV4.mat')
% targetroiNames: a string containing the name of target gray matter ROI .mat
%                 files
% threshold:  a threshold on the squared distance between gray matter voxel and VOF fiber termination point.
%            Default is 9. (3 mm distance in ACPC space)
% nboots: number of bootstrap. Default is 100000.
% nmontecarlo: number of the repetition of montecarlo simulation. Default is 5.
%
% EXAMPLE:
% feFileToLoad{1} = '/home/vof/data/S1/life/LH_1stconnectome_fe.mat';
% feFileToLoad{2} = '/home/vof/data/S1/life/LH_2ndconnectome_fe.mat';
% feFileToLoad{3} = '/home/vof/data/S1/life/LH_3rdconnectome_fe.mat';
% fibername = 'LVOF';
% roisDir = '/home/vof/data/S1/ROIs/Waypoint/'
% roiNames = {'LH_VOFPlane1.mat','LH_VOFPlane2.mat'};
% roiOperations = {'and','and'};
% targetroidir ='/home/vof/data/S1/ROIs/map/';
% targetroiNames = {'LV3AB.mat','LV3d.mat','LV2d.mat','LV1.mat','LV2v.mat','LV3v.mat','LhV4.mat','LVO1.mat','LLO1.mat','LLO2.mat'};
% vof_strengthofevidence_ROIdistance(feFileToLoad,fibername,roisDir,roiNames,roiOperations,targetroidir,targetroiNames);
%
% Copyright (C) Hiromasa Takemura, 2013, Stanford VISTA Team

% Argument checking
if notDefined('threshold')
    threshold = 9;
end
if notDefined('roiOperations')
    roiOperations = {'and','and'};
end
if notDefined('nboots')
    nboots = 100000;
end
if notDefined('nmontecarlo')
    nmontecarlo = 5;
end

sizecon = size(feFileToLoad);
numtargetROIs = size(targetroiNames);

% Run the test for multiple gray matter ROIs
for kk = 1:numtargetROIs(2)
    
    % Specifying the location and name of target gray matter ROIs
    targetROIfile = fullfile(targetroidir, targetroiNames{kk});
    [targetROI,remain] = strtok(targetroiNames{kk}, '.');
    
    % Run the test for multiple connectomes
    for connum = 1:sizecon(2)
        
        % Define the name of output .mat file
        fibernamenum = [fibername, num2str(connum)];
        savematfile = [fibername, targetROI,'_threshold',num2str(threshold),'_graymatter_S.mat'];
        
        %% Segmentation of the VOF
        % Load waypoint ROIs
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
        
        % Extract the fiber group from the FE structure
        fg = feGet(fe,'fibers acpc');
        
        % Segmenting fibers from connectomes
        fprintf('Segmenting the Tract from Connectome ...\n')
        [fgsegment, keepFascicles] = feSegmentFascicleFromConnectome(fg, rois, roiOperations, 'prob connectome');
        
        % Removing outlier fibers from VOF
        fprintf('Excluding outliers from Fasciculus ...\n')
        [fgsegment2 keepFascicles2] = ht_gradient_removeoutlier(fe, keepFascicles,2,rois{1},rois{2});
        [fgsegment3, keepFascicle3] = mbaComputeFibersOutliers(fgsegment2,3,3,25);
        
        keepFascicles4 = keepFascicles2;
        keepFascicle2matrix = find(keepFascicles2);
        
        for kj = 1:length(keepFascicle3)
            if keepFascicle3(kj)==0,
                keepFascicles4(keepFascicle2matrix(kj)) = 0;
            else
            end
        end
        
    %% Select fibers within a certain distance from each ROIs
    
    voffiberlength(connum) = length(fgsegment3.fibers);
    
    % Extract ACPC coordinate of dorsal and ventral VOF endpoint
    for kk = 1:voffiberlength(connum)
        fibercoordinate = cell2mat(fgsegment3.fibers(kk));
        fiberlength = size(fibercoordinate);
        if  fibercoordinate(3,1) > fibercoordinate(3,fiberlength(2))
            dorsalfiberend(kk,1,connum) = fibercoordinate(1,1);
            dorsalfiberend(kk,2,connum) = fibercoordinate(2,1);
            dorsalfiberend(kk,3,connum) = fibercoordinate(3,1);
            ventralfiberend(kk,1,connum) = fibercoordinate(1,fiberlength(2));
            ventralfiberend(kk,2,connum) = fibercoordinate(2,fiberlength(2));
            ventralfiberend(kk,3,connum) = fibercoordinate(3,fiberlength(2));
            
        else
            ventralfiberend(kk,1,connum) = fibercoordinate(1,1);
            ventralfiberend(kk,2,connum) = fibercoordinate(2,1);
            ventralfiberend(kk,3,connum) = fibercoordinate(3,1);
            dorsalfiberend(kk,1,connum) = fibercoordinate(1,fiberlength(2));
            dorsalfiberend(kk,2,connum) = fibercoordinate(2,fiberlength(2));
            dorsalfiberend(kk,3,connum) = fibercoordinate(3,fiberlength(2));
            
        end
        clear fibercoordinate fiberlength
    end
    
    % Load target ROI
    load(fullfile(targetROIdir,targetROIfile));
    roicoords = transpose(roi.coords);
    dorsalendcoords = transpose(dorsalfiberend(:,:,connum));
    ventralendcoords = transpose(ventralfiberend(:,:,connum));
    roicoordsize = size(roicoords);
    
    % Compute squared distance between ROI and each fascicle endpoints
    for kp = 1:voffiberlength(connum)
        [indices_d(kp), bestSqDis_d(kp)] = nearpoints(dorsalendcoords(:,kp), roicoords);
        [indices_v(kp), bestSqDis_v(kp)] = nearpoints( ventralendcoords(:,kp),roicoords);
    end
    
    
    % Select fibers within the threshold, then create keepFascicles
    % structure
    dorsalfiber = find(bestSqDis_d <threshold);
    ventralfiber = find(bestSqDis_v <threshold);
    fascicleidxtotal = find(keepFascicles4);
    dorsalsize = size(dorsalfiber);
    ventralsize = size(ventralfiber);
    
    fprintf('Testing the dorsal tract projection ...\n')
    if dorsalsize(2)>0
        keepFascicles5d = zeros(size(keepFascicles4));
        for ik = 1:length(dorsalfiber)
            keepFascicles5d(fascicleidxtotal(dorsalfiber(ik))) = 1;
            
        end
        [se_d{connum}] = feVirtualLesion(fe, logical(keepFascicles5d));
    else 
        se_d{connum} = 0;
    end
    
    fprintf('Testing the ventral tract projection ...\n')
    
    if ventralsize(2)>0
        keepFascicles5v = zeros(size(keepFascicles4));
        for ik = 1:length(ventralfiber)
            keepFascicles5v(fascicleidxtotal(ventralfiber(ik))) = 1;
            
        end
        [se_v{connum}] = feVirtualLesion(fe, logical(keepFascicles5v));
    else 
        se_v{connum} = 0;

    end
    
    clear keepFascicles4 keepFascicles3 keepFascicles2 keepFascicles keepFascicles5v keepFascicles5d fascicleidxtotal ventralfiber dorsalfiber bestSqDis_d bestSqDis_v indices_d indices_v    
    end
    % save files
    save(savematfile,'se_d','se_v');
    clear roicoords roicoordsize
end
