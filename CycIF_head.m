function [morph, fluorescence] = CycIF_head(imageDirectory, saveDirectory,...
    experiment, timepoint, mag, maxCycle, FOVlimits, channelNames,...
    bugGFP, bugmCherry, punctaChannels)

rows =["A", "B", "C", "D", "E", "F", "G", "H"];
columns =["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11"];
fields = ["1", "2", "3", "4", "5", "6", "7", "8", "9"];
channels = ["UV - DAPI", "Blue - FITC", "Green - dsRed", "Red - Cy5"];

%% cycle through rows. columns, fields of view, call image input, segment, quantify functions 
for r = FOVlimits(1):FOVlimits(2); %choose rows
    for c = FOVlimits(3):FOVlimits(4) %choose columns
        for f = FOVlimits(5):FOVlimits(6) %choose fields
            row = rows(r);
            column = columns(c);
            field = fields(f);
            name = char(strcat(experiment, timepoint, mag, row,...
                column, '_fld', field));
            FOV = char(strcat(rows(r), columns(c), '0', fields(f)));
            FOVstack = CycIFinputTiffStack(imageDirectory, experiment,...
                timepoint, mag, maxCycle, row, column, field);
            [nuclei, nucleiShrink, nucleiExpand, waterMF] = CycIFNucSeg(FOVstack, maxCycle); 
            if max(nuclei(:)) == 0
                continue
            end
            [cells, cytoplasm] = CycIFCellSeg(FOVstack, nuclei, nucleiExpand, waterMF);
            [bugs, bugsCellLabel] = CycIFBugSeg(FOVstack, maxCycle, cells);
            featureData = CycIFData(FOVstack, maxCycle, channelNames, nuclei,...
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
%             for m=1:length(masks(1, 1, :))
%                 imwrite(masks(:, :, m), strcat(saveDirectory, name, '_masks.tif'), 'WriteMode', 'append',  'Compression','none');
%             end
            text = strcat('finished file_', name);
            clearvars name FOVstack nuclei nucleiShrink nucleiExpand...
                waterMF cells cytoplasm bugs bugsCellLabel featureData...
                masks cellStats cytoplasmStats nucleiStats bugRatio...
                bugCellStats
        end %field
    end %columns
end %rows