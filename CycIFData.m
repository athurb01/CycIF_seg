function [featureData] = CycIFData(FOVstack, maxCycle, channelNames,...
    nuclei, nucleiShrink, cells, cytoplasm, bugs, bugsCellLabel,...
    saveDirectory, name, FOV, bugGFP, bugmCherry, punctaChannels) 
%% extract data
% background subtraction - current manual value
for c= 2:(4*maxCycle)
    FOVstack(:,:, c) = imsubtract(FOVstack(:,:,c), 110);
end

% slice that is GFP/mCherry
bugRatio = double(FOVstack(:,:,10))./double(FOVstack(:,:,19));

% object property data, once per object

numCells = max(nuclei(:));
properties = {'Area', 'BoundingBox', 'Centroid', 'Eccentricity', 'MajorAxisLength', 'MinorAxisLength', 'Perimeter', 'Solidity'};
nucleiStats = regionprops(nuclei, properties);
cellStats = regionprops(cells, properties);
cytoplasmStats = regionprops(cytoplasm, properties);
bugsCellStats = regionprops(bugsCellLabel, properties);
%pad end of bugsCellStats in case last few cells have 0 bugs
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
for ch = maxCycle+1:(maxCycle*4)
    currentChannel = im2double(FOVstack(:,:,ch));
    % texture filters
    rangeIm = rangefilt(currentChannel);
    stdIm = uint16(stdfilt(currentChannel));
    entropyIm = uint16(entropyfilt(currentChannel));
    % puncta analysis
    if ismember(ch, punctaChannels)
        punctaMask = CycIFFoci(currentChannel, nuclei, cells);
    else
        punctaMask = zeros(2048,2048);
    end
    for countCell = startCount:(endCount) %countCell is row number of featureData
        %metadata
        featureData(countCell).experiment = name(1:4);
        featureData(countCell).timepoint = name(6:8);
        featureData(countCell).well = name(14:16);
        featureData(countCell).row = name(14);
        featureData(countCell).column = name(15:16);
        featureData(countCell).field = name(end);
        featureData(countCell).object = obj;
        featureData(countCell).channel = channelNames(ch);
        % nuclei data
        featureData(countCell).meanNuclei = mean(currentChannel(nuclei == obj));
        featureData(countCell).modeNuclei = mode(currentChannel(nuclei == obj));
        featureData(countCell).medianNuclei = median(currentChannel(nuclei == obj));
        featureData(countCell).intIntensityNuclei =...
            featureData(countCell).meanNuclei * nucleiStats(obj).Area;
        featureData(countCell).meanRangeNuclei = mean(rangeIm(nuclei == obj));
        featureData(countCell).meanStdNuclei = mean(stdIm(nuclei == obj));
        featureData(countCell).meanEntropyNuclei = mean(entropyIm(nuclei == obj));
        % cell data
        featureData(countCell).meanCell = mean(currentChannel(cells == obj));
        featureData(countCell).modeCell = mode(currentChannel(cells == obj));
        featureData(countCell).medianCell = median(currentChannel(cells == obj));
        featureData(countCell).intIntensityCell =...
            featureData(countCell).meanCell * cellStats(obj).Area;
        featureData(countCell).meanRangeCell = mean(rangeIm(cells == obj));
        featureData(countCell).meanStdCell = mean(stdIm(cells == obj));
        featureData(countCell).meanEntropyCell = mean(entropyIm(cells == obj));
        % cytoplasm data
        featureData(countCell).meanCytoplasm = mean(currentChannel(cytoplasm == obj));
        featureData(countCell).modeCytoplasm = mode(currentChannel(cytoplasm == obj));
        featureData(countCell).medianCytoplasm = median(currentChannel(cytoplasm == obj));
        featureData(countCell).intIntensityCytoplasm =...
         featureData(countCell).meanCytoplasm * cytoplasmStats(obj).Area;
        featureData(countCell).meanRangeCytoplasm = mean(rangeIm(cytoplasm == obj));
        featureData(countCell).meanStdCytoplasm = mean(stdIm(cytoplasm == obj));
        featureData(countCell).meanEntropyCytoplasm = mean(entropyIm(cytoplasm == obj));
        % nuclei/cytoplasm ratio data
        featureData(countCell).meanNC = featureData(countCell).meanNuclei...
            /featureData(countCell).meanCytoplasm;
        featureData(countCell).modeNC = featureData(countCell).modeNuclei...
            /featureData(countCell).modeCytoplasm;
        featureData(countCell).medianNC = featureData(countCell).medianNuclei...
            /featureData(countCell).medianCytoplasm;
        featureData(countCell).intIntensityNC = featureData(countCell).intIntensityNuclei...
            /featureData(countCell).intIntensityCytoplasm;
        % bug per cell data
        if ch == bugGFP | ch == bugmCherry
            featureData(countCell).meanBug = mean(currentChannel(bugsCellLabel == obj));
            featureData(countCell).intIntensityBug = ...
                featureData(countCell).meanBug * bugsCellStats(obj).Area;
            featureData(countCell).meanRatioBug = mean(bugRatio(bugsCellLabel == obj));
        else
            featureData(countCell).meanBug = 0;
            featureData(countCell).intIntensityBug = 0;
            featureData(countCell).meanRatioBug = 0;
        end
        %puncta
        featureData(countCell).punctaNuclei = sum(punctaMask(nuclei == obj));
        featureData(countCell).punctaCell = sum(punctaMask(cell == obj));
        obj = obj+1; 

    end
    startCount = startCount + numCells;
    endCount = endCount + numCells;
    obj = 1;
end



%save data
writetable(struct2table(featureData), strcat(saveDirectory, name, '.txt'));
