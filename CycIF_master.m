%% master script to input images from Tiff stack, segment, and extract data
%assumes tiff script are in order by channel (not by cycle)

%% inputs

%inputs that change with every batch
imageDirectory = 'F:\FI11\Registered2\';
saveDirectory = 'C:\Users\Amy Thurber\Dropbox (Partners HealthCare)\Experiments\FI11\Matlab analysis\';
experiment = 'FI11_';
timepoint = '3h_';
mag = '20x_';
maxCycle = 5;

%inputs that change with every experiment
channelNames = { 'Hoechst1', 'Hoechst2', 'Hoechst3',...
    'Hoechst4','Hoechst5', 'pH_GFP',...
    'p_TBK1', 'p_cJun', 'Arg-1', 'p_p38', 'mTB_mCherry',...
    'p_ERK', 'p_JNK', 'CD68', 'p_STAT1',...
    'mitoRed', 'p_STAT6', 'p_mTOR', 'RelA', 'p_STAT3'};

%inputs consistent between experiment
rows =["A", "B", "C", "D", "E", "F", "G", "H"];
columns =["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11"];
fields = ["1", "2", "3", "4", "5", "6", "7", "8", "9"];
channels = ["UV - DAPI", "Blue - FITC", "Green - dsRed", "Red - Cy5"];

%% cycle through rows. columns, fields of view, call image input, segment, quantify functions 
for r = 2:2; %choose rows
    for c = 5:5 %choose columns
        for f = 1:3 %choose fields
            % REMOVE 0 IF TIMEPOINT IS ALREADY PADDED
            name = char(strcat(experiment, '0', timepoint, mag, rows(r),...
                columns(c), '_fld', fields(f)));
            FOV = char(strcat(rows(r), columns(c), '0', fields(f)));
            FOVstack = CycIFinputTiffStack(imageDirectory, experiment,...
                timepoint, mag, maxCycle, rows(r), columns(c), fields(f));
            [nuclei, nucleiShrink, nucleiExpand, waterMF] = CycIFNucSeg(FOVstack, maxCycle); 
            if max(nuclei(:)) == 0
                continue
            end
            [cells, cytoplasm] = CycIFCellSeg(FOVstack, nuclei, nucleiExpand, waterMF);
            [bugs, bugsCellLabel] = CycIFBugSeg(FOVstack, maxCycle, cells);
            featureData = CycIFData(FOVstack, maxCycle, channelNames, nuclei,...
                nucleiShrink, cells, cytoplasm, bugs, bugsCellLabel,...
                saveDirectory, name, FOV);
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
