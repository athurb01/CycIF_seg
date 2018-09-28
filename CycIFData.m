function [morphology, featureData] = CycIFData(FOVstack, maxCycle, channelNames,...
    nuclei, nucleiShrink, cells, cytoplasm, bugs, bugsCellLabel,...
    saveDirectory, name, FOV, bugGFP, bugmCherry, punctaChannels,...
    experiment, timepoint, mag, row, column, field) 
%% extract data
% background subtraction - median of background (non-cell) pixels
for c= 1:(4*maxCycle)
    temp = FOVstack(:,:,c);
    background = double(median(temp(cells == 0)));
    FOVstack(:,:, c) = imsubtract(FOVstack(:,:,c), background);
end

% image that is GFP/mCherry
IndexCh = strfind(channelNames, bugGFP(1));
idx_GFP_Z = find(not(cellfun('isempty', IndexCh)));
IndexCh = strfind(channelNames, bugGFP(2));
idx_GFP_slice = find(not(cellfun('isempty', IndexCh)));
IndexCh = strfind(channelNames, bugmCherry(1));
idx_mCher_Z = find(not(cellfun('isempty', IndexCh)));
IndexCh = strfind(channelNames, bugmCherry(2));
idx_mCher_slice = find(not(cellfun('isempty', IndexCh)));

bugRatio_Z = double(FOVstack(:,:,idx_GFP_Z))./double(FOVstack(:,:,idx_mCher_Z));
bugRatio_slice = double(FOVstack(:,:,idx_GFP_slice))./double(FOVstack(:,:,idx_mCher_slice));

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
    cellStats(i).timepoint = timepoint;
end

% convert morphology data to table and merge
nucleiStatsTab = struct2table(nucleiStats);
nucleiStatsTab.Properties.VariableNames = {'NucArea' 'NucBoundingBox' 'NucCentroid'...
    'NucEccentricity' 'NucMajorAxisLength' 'NucMinorAxisLength'...
    'NucPerimeter' 'NucSolidity' 'NucFormFactor' 'CellNumber' 'FOV'};
cellStatsTab = struct2table(cellStats);
cellStatsTab.Properties.VariableNames = {'CellArea' 'CellBoundingBox' 'CellCentroid'...
    'CellEccentricity' 'CellMajorAxisLength' 'CellMinorAxisLength'...
    'CellPerimeter' 'CellSolidity' 'CellFormFactor' 'CellNumber' 'timepoint'};
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
for ch = 1:(maxCycle*4)
    currentChannel = im2double(FOVstack(:,:,ch));
    currentChannel = (currentChannel.*65535)+1;%convert back to original pixel values and add 1 to get rid of 0s
    
    % texture filters
    rangeIm = rangefilt(currentChannel);
    stdIm = uint16(stdfilt(currentChannel));
    entropyIm = uint16(entropyfilt(currentChannel));
    
    % puncta analysis
    if ismember(channelNames(ch), punctaChannels)
        punctaMask = CycIFFoci(currentChannel, nuclei, cells);
    else
        punctaMask = zeros(2048,2048);
    end %puncta end
    
    for countCell = startCount:(endCount) %countCell is row number of featureData
        %metadata
        featureData(countCell).experiment = experiment;
        featureData(countCell).timepoint = timepoint;
        featureData(countCell).well = row + column;
        featureData(countCell).row = row;
        featureData(countCell).column = column;
        featureData(countCell).field = field;
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
        if ismember(channelNames(ch), bugGFP(1)) | ismember(channelNames(ch), bugmCherry(1))
            featureData(countCell).meanBug = mean(currentChannel(bugsCellLabel == obj));
            featureData(countCell).medianBug = median(currentChannel(bugsCellLabel == obj));
            featureData(countCell).intIntensityBug = ...
                featureData(countCell).meanBug * bugsCellStats(obj).Area;
            featureData(countCell).meanRatioBug = mean(bugRatio_Z(bugsCellLabel == obj));
            featureData(countCell).medianRatioBug = median(bugRatio_Z(bugsCellLabel == obj));
        elseif ismember(channelNames(ch), bugGFP(2)) | ismember(channelNames(ch), bugmCherry(2))
            featureData(countCell).meanBug = mean(currentChannel(bugsCellLabel == obj));
            featureData(countCell).medianBug = median(currentChannel(bugsCellLabel == obj));
            featureData(countCell).intIntensityBug = ...
                featureData(countCell).meanBug * bugsCellStats(obj).Area;
            featureData(countCell).meanRatioBug = mean(bugRatio_slice(bugsCellLabel == obj));
            featureData(countCell).medianRatioBug = median(bugRatio_slice(bugsCellLabel == obj));
        else
            featureData(countCell).meanBug = 0;
            featureData(countCell).medianBug = 0;
            featureData(countCell).intIntensityBug = 0;
            featureData(countCell).meanRatioBug = 0;
            featureData(countCell).medianRatioBug = 0;
        end %bug end
        %puncta
        featureData(countCell).punctaNuclei = sum(punctaMask(nuclei == obj));
        featureData(countCell).punctaCell = sum(punctaMask(cells == obj));
        obj = obj+1; 

    end %countCell end
    startCount = startCount + numCells;
    endCount = endCount + numCells;
    obj = 1;
end %channel end



%save data
writetable(struct2table(featureData), strcat(saveDirectory, name, '.txt'));
