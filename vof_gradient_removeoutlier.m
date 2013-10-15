function [fgclean keepFascicles] = vof_gradient_removeoutlier(fe, keepFG,maxSD,roi1,roi2)

% This function removes the outlier based on the orientation of VOF fibers
% between two waypoint ROIs. First, this code computes the gradient of
% fibers in X, Y and Z coordinates. Second, it computes the mean and standard deviation of each gradient. Then the code removes the outlier fibers (e.g. 2SD away from mean gradient in any X, Y, Z coordinate).
%
% INPUT:
% fe: fe structure (created by LIFE code, which is available at https://github.com/francopestilli/life).
% keepFG: one dimensional matrix describing which fibers in fe sturcture are included in the analysis. 1: included, 0: excluded.
% maxSD: a threshod for fiber orientation. The unit is standard deviation of gradient (X, Y and Z coordinates) from the mean gradient. Fibers with higher deviation from this criterion are removed.
% roi1: first waypoint ROI for VOF segmentation. Mat format in mrDiffusion is required.
% roi2: second waypoint ROI for VOF segmentation.
%
% OUTPUT:
% fgclean: cleaned fiber groups.
% keepFascicles: one dimensional matrix describing which fibers are included as "cleaned" VOF. 1: included, 0: excluded.
%
% Hiromasa Takemura, (c) Stanford Vista Team, 2013

%% Extract fibers
fefiber = feGet(fe,'fibers acpc');
fg= fgExtract(fefiber, transpose(keepFG), 'keep');
keepFGindice = find(keepFG);

%% Read waypoint rois and extract z-coordinate for reference
roifirst = dtiReadRoi(roi1);
roisecond = dtiReadRoi(roi2);
roifirstcoord = roifirst.coords(1,:);
roisecondcoord = roisecond.coords(1,:);

%% Define nearest point between fibers and reference z-coordinate
for i=1:size(fg.fibers)
    
    currfiber = fg.fibers{i};
    currfiberz = currfiber(3,:);
    sizecurrfiberz = size(currfiberz);
    if currfiberz(1) > currfiberz(sizecurrfiberz(2))
        currfiber = fliplr(currfiber);
        currfiberz = fliplr(currfiberz);
    else
    end
    
    % Define nearest point
    temp1 = abs(currfiberz - roifirstcoord(3));
    [idx1(i) idx1(i)] = min(temp1);
    temp2 = abs(currfiberz - roisecondcoord(3));
    [idx2(i) idx2(i)] = min(temp2);
    
    % Compute gradients between two reference point
    gradientcurr =  gradient(currfiber);
    if idx1(i) < idx2(i)
        gradientcurrx = gradientcurr(1,idx1(i):idx2(i));
        gradientcurry = gradientcurr(2,idx1(i):idx2(i));
        gradientcurrz = gradientcurr(3,idx1(i):idx2(i));
    else
        gradientcurrx = gradientcurr(1,idx2(i):idx1(i));
        gradientcurry = gradientcurr(2,idx2(i):idx1(i));
        gradientcurrz = gradientcurr(3,idx2(i):idx1(i));
    end
    meangradientx(i) = mean(gradientcurrx);
    meangradienty(i) = mean(gradientcurry);
    meangradientz(i) = mean(gradientcurrz);
end

%% Calculate mean of the gradient across fibers, calculate standard deviation, remove outlier
median_x = median(meangradientx);
median_y = median(meangradienty);
median_z = median(meangradientz);
std_x = std(meangradientx);
std_y = std(meangradienty);
std_z = std(meangradientz);

dev_x = abs(meangradientx - median_x);
dev_y = abs(meangradienty - median_y);
dev_z = abs(meangradientz - median_z);

index_x = find(dev_x > maxSD*std_x);
index_y = find(dev_y > maxSD*std_y);
index_z = find(dev_z > maxSD*std_z);

keepmatrix = ones(size(fg.fibers));
for jj=1:length(index_x)
    keepmatrix(index_x(jj)) = 0;
end
for pp=1:length(index_y)
    keepmatrix(index_y(pp)) = 0;
end
for rr=1:length(index_z)
    keepmatrix(index_z(rr)) = 0;
end

%% Compute proper keepindices for subsequent analyses
keepFascicles = zeros(size(fefiber.fibers));
for kp = 1:length(keepFGindice)
    if keepmatrix(kp) == 1,
        keepFascicles(keepFGindice(kp)) = 1;
    else
    end
end

keepFascicles = logical(keepFascicles);

%% Write fibers
fgclean= fgExtract(fefiber, transpose(keepFascicles), 'keep');
