%% master script to input images from Tiff stack, segment, and extract data
%assumes tiff script are in order by channel (not by cycle)

%% inputs

%inputs that change with every batch
imageDirectory = 'F:\FI13\AT_FI13_conf\AT_FI13_conf_rd2_1';
saveDirectory = 'C:\Users\Amy Thurber\Dropbox (Partners HealthCare)\Experiments\FI13_matlab_out\matlab_output\';
experiment = 'FI13_';
timepoint = '03h_';
mag = '20x_';

%inputs that change with every experiment
channelNames = {'Brightfield', 'Hoechst1', 'Hoechst2', 'Hoechst3',...
    'Hoechst4','Hoechst5', 'Hoechst6', 'Hoechst7', 'Hoechst8', 'pH_GFP', 'p_p38',...
    'p_TBK1', 'Bax', 'Nos2', 'p_cJun', 'none', 'none', 'mTB_mCherry', 'p_STAT1',...
    'p_ERK', 'none', 'CD68', 'none', 'actin', 'p50', 'mitoRed', 'p_STAT3',...
    'p_STAT6', 'Bcl2', 'RelA', 'LC_3', 'Ki67', 'none', 'morphology'};

%inputs consistent between experiment
rows =["A", "B", "C", "D", "E", "F", "G", "H"];
columns =["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11"];
fields = ["1", "2", "3", "4", "5", "6", "7", "8", "9"];
channels = ["UV - DAPI", "Blue - FITC", "Green - dsRed", "Red - Cy5"];

%% cycle through rows. columns, fields of view, call image input, segment, quantify functions 
for r = 2:7; %choose rows
    for c = 5:5 %choose columns
        for f = 1:3 %choose fields
            % REMOVE 0 IF TIMEPOINT IS ALREADY PADDED
            name = char(strcat(experiment, timepoint, mag, rows(r),...
                columns(c), '_fld', fields(f), 'rd2'));
            FOV = char(strcat(rows(r), columns(c), '0', fields(f)));
            FOVstack = CycIFinput_1cycle(imageDirectory, channels, rows(r)...
                , columns(c), fields(f));
            [nuclei, nucleiShrink, nucleiExpand, waterMF] = ...
                CycIFNucSeg_1cyc(FOVstack); 
            if max(nuclei(:)) < 5
                clearvars name FOVstack nuclei nucleiShrink nucleiExpand...
                waterMF
                continue
            end
            [cells, cytoplasm] = CycIFCellSeg_1cyc(FOVstack, nuclei,...
                nucleiExpand, waterMF);
            [bugs, bugsCellLabel] = CycIFBugSeg_1cyc(FOVstack, cells);
            featureData = CycIFData_1cyc(FOVstack,...
                nuclei, nucleiShrink, cells, cytoplasm, bugs, bugsCellLabel,...
                saveDirectory, name, FOV, channels);
            nucleiEdge = edge(nuclei>0);
            cellEdge = edge(cells>0);
            bugEdge = edge(bugs>0);
            cytoEdge = nucleiEdge + cellEdge;
            allEdge = cytoEdge + bugEdge;
            masks = cat(3, nuclei, nucleiShrink, nucleiExpand, cells,...
                cytoplasm, bugs, bugsCellLabel,...
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
