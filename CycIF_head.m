function [morph, fluorescence] = CycIF_head(imageDirectory, saveDirectory,...
    experiment, timepoint, row, column, mag, maxCycle, FOVlimits, channelNames,...
    bugGFP, bugmCherry, punctaChannels)
fields = ["1", "2", "3", "4", "5", "6", "7", "8", "9"];
channels = ["UV - DAPI", "Blue - FITC", "Green - dsRed", "Red - Cy5"];

%% cycle through rows. columns, fields of view, call image input, segment, quantify functions 
for f = FOVlimits(7):FOVlimits(8) %choose fields
    field = fields(f);
    name = char(strcat(experiment, timepoint, mag, '_', row,...
        column, '_fld', field));
    FOV = char(strcat(row, column, '0', field));
    FOVstack = CycIFinputTiffStack(imageDirectory, experiment,...
        timepoint, mag, maxCycle, row, column, field);
    [nuclei, nucleiShrink, nucleiExpand, waterMF] = CycIFNucSeg(FOVstack, maxCycle); 
    if max(nuclei(:)) == 0
        continue
    end
    [cells, cytoplasm] = CycIFCellSeg(FOVstack, nuclei, nucleiExpand, waterMF, channelNames);
    [bugs, bugsCellLabel] = CycIFBugSeg(FOVstack, maxCycle, cells);
    [morph,fluorescence] = CycIFData(FOVstack, maxCycle, channelNames, nuclei,...
        nucleiShrink, cells, cytoplasm, bugs, bugsCellLabel,...
        saveDirectory, name, FOV, bugGFP, bugmCherry, punctaChannels,...
        experiment, timepoint, mag,row, column, field);
     nucleiEdge = edge(nuclei>0);
     cellEdge = edge(cells>0);
     bugEdge = edge(bugs>0);
     cytoEdge = nucleiEdge + cellEdge;
     allEdge = cytoEdge + bugEdge;
     masks = cat(3, nuclei, nucleiShrink, nucleiExpand, cells, cytoplasm, bugs, bugsCellLabel,...
        nucleiEdge, cellEdge, cytoEdge, bugEdge, allEdge);
     edgeMasks = cat(3, nucleiEdge, cellEdge, cytoEdge, bugEdge, allEdge);
            for m=1:length(masks(1, 1, :))
                imwrite(masks(:, :, m), strcat(saveDirectory, name, '_masks.tif'), 'WriteMode', 'append',  'Compression','none');
            end
     clearvars name FOVstack nuclei nucleiShrink nucleiExpand...
                waterMF cells cytoplasm bugs bugsCellLabel featureData...
                masks cellStats cytoplasmStats nucleiStats bugRatio...
                bugCellStats
end %field
