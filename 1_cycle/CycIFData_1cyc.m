function [featureData] = CycIFData(FOVstack, nuclei, nucleiShrink, cells,...
    cytoplasm, bugs, bugsCellLabel, saveDirectory, name, FOV, channels,...
    experiment, timepoint, row, column,field) 
%% extract data
% background subtraction - current manual value
for c = 1:3
    FOVstack(:,:, c) = imsubtract(FOVstack(:,:,c), 110);
end

% slice that is GFP/mCherry
%bugRatio = double(FOVstack(:,:,2))./double(FOVstack(:,:,3));

% object property data, once per object

numCells = max(nuclei(:));
properties = {'Area', 'BoundingBox', 'Centroid', 'Eccentricity',...
    'MajorAxisLength', 'MinorAxisLength', 'Perimeter', 'Solidity'};
nucleiStats = regionprops(nuclei, properties);
cellStats = regionprops(cells, properties);
cytoplasmStats = regionprops(cytoplasm, properties);
bugsCellStats = regionprops(bugsCellLabel, properties);
%pad end of bugsCellStats in case last few cells have 0 bugs
% need to fix so it works if there is nuclei but no cell
for i = 1:numCells
    if length(bugsCellStats) < length(nucleiStats)
        bugsCellStats(end+1).Area = 0;
    else
        break
    end
end
%same for cellStats and cytoplasmStats in case nuclei is on very edge of image
for i = 1:numCells
    if length(cellStats) < length(nucleiStats)
        cellStats(end+1).Area = 0;
    else
        break
    end
end
for i = 1:numCells
    if length(cytoplasmStats) < length(nucleiStats)
        cytoplasmStats(end+1).Area = 0;
    else
        break
    end
end

% add additional columns to morphology tables
% add form factor = perimeter/area for each compartment
for i = 1:numCells
    nucleiStats(i).FormFactor = nucleiStats(i).Perimeter/nucleiStats(i).Area;
    nucleiStats(i).CellNumber = i;
    nucleiStats(i).FOV = FOV;
    cellStats(i).FormFactor = cellStats(i).Perimeter/cellStats(i).Area;
    cellStats(i).CellNumber = i;
    cytoplasmStats(i).FormFactor = cytoplasmStats(i).Perimeter/cytoplasmStats(i).Area;
    cytoplasmStats(i).CellNumber = i;
    bugsCellStats(i).FormFactor = bugsCellStats(i).Perimeter/bugsCellStats(i).Area;
    bugsCellStats(i).CellNumber = i;
end

% convert morphology data to table and merge
nucleiStatsTab = struct2table(nucleiStats);
nucleiStatsTab.Properties.VariableNames = {'NucArea' 'NucBoundingBox' 'NucCentroid'...
    'NucEccentricity' 'NucMajorAxisLength' 'NucMinorAxisLength'...
    'NucPerimeter' 'NucSolidity' 'NucFormFactor' 'CellNumber' 'FOV'};
cellStatsTab = struct2table(cellStats);
cellStatsTab.Properties.VariableNames = {'CellArea' 'CellBoundingBox' 'CellCentroid'...
    'CellEccentricity' 'CellMajorAxisLength' 'CellMinorAxisLength'...
    'CellPerimeter' 'CellSolidity' 'CellFormFactor' 'CellNumber'};
cytoplasmStatsTab = struct2table(cytoplasmStats);
cytoplasmStatsTab.Properties.VariableNames = {'CytArea' 'CytBoundingBox' 'CytCentroid'...
    'CytEccentricity' 'CytMajorAxisLength' 'CytMinorAxisLength'...
    'CytPerimeter' 'CytSolidity' 'CytFormFactor' 'CellNumber'};
bugsCellStatsTab = struct2table(bugsCellStats);
bugsCellStatsTab.Properties.VariableNames = {'BugArea' 'BugBoundingBox' 'BugCentroid'...
    'BugEccentricity' 'BugMajorAxisLength' 'BugMinorAxisLength'...
    'BugPerimeter' 'BugSolidity' 'BugFormFactor' 'CellNumber'};

morphology = join(nucleiStatsTab, cellStatsTab);
morphology = join(morphology, cytoplasmStatsTab);
morphology = join(morphology, bugsCellStatsTab);

%save morphology data
writetable(morphology, strcat(saveDirectory, 'Morph', name, '.txt'));

% features extracted from every channel
startCount = 1;
endCount = numCells; 
obj = 1; %keeps track of cell number

% cycle through each channel, cycle through each cell (obj) and extract
% all data features
for c = 1:3
    currentChannel = im2double(FOVstack(:,:,c));
    for countCell = startCount:(endCount) %countCell is row number of featureData
        %metadata
        featureData(countCell).experiment = experiment;
        featureData(countCell).timepoint = timepoint;
        featureData(countCell).well = row + column;
        featureData(countCell).row = row;
        featureData(countCell).column = column;
        featureData(countCell).field = field;
        featureData(countCell).round = 1;
        featureData(countCell).object = obj;
        featureData(countCell).channel = channels(c);
        % nuclei data
        featureData(countCell).meanNuclei = mean(currentChannel(nuclei == obj));
        featureData(countCell).modeNuclei = mode(currentChannel(nuclei == obj));
        featureData(countCell).medianNuclei = median(currentChannel(nuclei == obj));
        featureData(countCell).intIntensityNuclei =...
            featureData(countCell).meanNuclei * nucleiStats(obj).Area;
        % cell data
        featureData(countCell).meanCell = mean(currentChannel(cells == obj));
        featureData(countCell).modeCell = mode(currentChannel(cells == obj));
        featureData(countCell).medianCell = median(currentChannel(cells == obj));
        featureData(countCell).intIntensityCell =...
            featureData(countCell).meanCell * cellStats(obj).Area;
        % cytoplasm data
        featureData(countCell).meanCytoplasm = mean(currentChannel(cytoplasm == obj));
        featureData(countCell).modeCytoplasm = mode(currentChannel(cytoplasm == obj));
        featureData(countCell).medianCytoplasm = median(currentChannel(cytoplasm == obj));
        featureData(countCell).intIntensityCytoplasm =...
         featureData(countCell).meanCytoplasm * cytoplasmStats(obj).Area;
        % nuclei/cytoplasm ratio data
        featureData(countCell).meanNC = featureData(countCell).meanNuclei...
            /featureData(countCell).meanCytoplasm;
        featureData(countCell).modeNC = featureData(countCell).modeNuclei...
            /featureData(countCell).modeCytoplasm;
        featureData(countCell).medianNC = featureData(countCell).medianNuclei...
            /featureData(countCell).medianCytoplasm;
        featureData(countCell).intIntensityNC = featureData(countCell).intIntensityNuclei...
            /featureData(countCell).intIntensityCytoplasm;
%         % bug per cell data
%         if c == 2 | c == 3
%         featureData(countCell).meanBug = mean(currentChannel(bugsCellLabel == obj));
%         featureData(countCell).medianBug = median(currentChannel(bugsCellLabel == obj));
%         featureData(countCell).intIntensityBug = ...
%             featureData(countCell).meanBug * bugsCellStats(obj).Area;
%         featureData(countCell).meanRatioBug = mean(bugRatio(bugsCellLabel == obj));
%         featureData(countCell).medianRatioBug = median(bugRatio(bugsCellLabel == obj));
%         featureData(countCell).areaBug = bugsCellStats(obj).Area;
%         else
%             featureData(countCell).meanBug = 0;
%             featureData(countCell).intIntensityBug = 0;
%             featureData(countCell).meanRatioBug = 0;
%         end
        
        
        obj = obj+1; 

    end
    startCount = startCount + numCells;
    endCount = endCount + numCells;
    obj = 1;
end



%save data
writetable(struct2table(featureData), strcat(saveDirectory, name, '.txt'));
