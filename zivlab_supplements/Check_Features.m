function [filtersFeatures,probICFlag,circularity,distance,isRegionInEdge] = Check_Features(filterMatrix,linIndsOfEdges,options,edges_flag)

% Noa- edges_flag was added in order to not discard cells in the edges

% filterMatrix - nxm matrix with the IC filter.
% Using region props, identify certain featurs to automaticly decide if the Ic is indeed a cell

properties ='all'; % {'Centroid','Area','BoundingBox'};

level = graythresh(filterMatrix); % Using otso method
bw = im2bw(filterMatrix ,level);
filtersFeatures = regionprops(bw,properties);

numOfROIs = numel(filtersFeatures);
if numOfROIs>1
    areas =cell2mat( {filtersFeatures.Area});
    [~,idxOfLargestArea] = max(areas);
    max_area=areas(idxOfLargestArea);
    other_areas=areas;
    other_areas(idxOfLargestArea)=[];
    if max_area*options.secRoiThrshld<sum(other_areas)
        numofROIsError = true;
    else
        numofROIsError = false;
    end
    filtersFeatures=filtersFeatures(idxOfLargestArea);
else
    numofROIsError = false;
end
lindIdxOfPixels =  filtersFeatures.PixelIdxList ;
if edges_flag
    isRegionInEdge = sum(ismember(lindIdxOfPixels,linIndsOfEdges))>0;
else
   isRegionInEdge = false;  
end

%|| madTraceVal > options.madThreshold ;

%            filtersFeatures.Eccentricity>options.maxEcc ||...\

% For cells at the periphery larger and less circular filters are allowed
center_x=round(size(filterMatrix,2)/2);
center_y=round(size(filterMatrix,1)/2);
max_distance=sqrt((2*center_x)^2+(2*center_y)^2);
locs=filtersFeatures.Centroid;
distance=sqrt((center_x-locs(1))^2+(center_y-locs(2))^2);

allowed_max_area=options.maxCellSize+options.maxCellSize*distance/max_distance;
allowed_min_area=options.minCellSize+options.minCellSize*distance/max_distance;

% circularity computed for a thresholded footprint:
level = 0.5; % Using otso method
thresholded_filterMatrix=filterMatrix>max(max(filterMatrix))*level;
% bw = im2bw(filterMatrix ,level);
thresholded_filtersFeatures = regionprops(thresholded_filterMatrix,properties);

thresholded_numOfROIs = numel(thresholded_filtersFeatures);
if thresholded_numOfROIs>1
    thresholded_areas =cell2mat( {thresholded_filtersFeatures.Area});
    [~,thresholded_idxOfLargestArea] = max(thresholded_areas);
    thresholded_filtersFeatures=thresholded_filtersFeatures(thresholded_idxOfLargestArea);
end
circularity = (4*pi*thresholded_filtersFeatures.Area)./((thresholded_filtersFeatures.Perimeter)^2);





allowed_circularity=options.circularityThreshold-options.circularityThreshold/3*distance/max_distance;

probICFlag = numofROIsError || isRegionInEdge || filtersFeatures.Area > allowed_max_area ...
    || filtersFeatures.Area < allowed_min_area || circularity < allowed_circularity;


end