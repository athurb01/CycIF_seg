%% master script to input images from Tiff stack, segment, and extract data
%assumes tiff script are in order by channel (not by cycle)

%% inputs

%inputs that change with every batch
imageDirectory = 'D:\fixed_7_24_18\AT_FI10_20x\AT_FI10_20x_1\';
saveDirectory = 'D:\fixed_7_24_18\AT_FI10_20x\';
experiment = 'FIReed_';
timepoint = '30min_';
mag = '20x_';

%inputs that change with every experiment
channelNames = {'Hoechst', 'p50', 'RelA', 'Hoechst3',...
    'Hoechst4','Hoechst5', 'Hoechst6', 'Hoechst7', 'Hoechst8', 'pH_GFP', 'p_p38',...
    'p_TBK1', 'Bax', 'Nos2', 'p_cJun', 'none', 'none', 'mTB_mCherry', 'p_STAT1',...
    'p_ERK', 'none', 'CD68', 'none', 'actin', 'p50', 'mitoRed', 'p_STAT3',...
    'p_STAT6', 'Bcl2', 'RelA', 'LC_3', 'Ki67', 'none', 'morphology'};

%inputs consistent between experiment
rows ={"A", "B", "C", "D", "E", "F", "G", "H"};
columns ={"01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11"};
fields = {"1", "2", "3", "4", "5", "6", "7", "8", "9"};
channels = {"UV - DAPI", "Green - dsRed", "Red - Cy5"};

%% cycle through rows. columns, fields of view, call image input, segment, quantify functions 
for r = 1:4; %choose rows
    for c = 1:11 %choose columns
        for f = 1:4 %choose fields
            % REMOVE 0 IF TIMEPOINT IS ALREADY PADDED
            row = string(rows(r));
            column = string(columns(c));
            field = string(fields(f));
            name = char(strcat(experiment, timepoint, mag, row,...
                column, '_fld', field, '_rd1'));
            FOV = row + column + '0' + field;
            FOVstack = CycIFinput_1cycle(imageDirectory, channels, row...
                , column, field);
            [nuclei, nucleiShrink, nucleiExpand, waterMF] = ...
                CycIFNucSeg_1cyc(FOVstack); 
            if max(nuclei(:)) < 5
                clearvars name FOVstack nuclei nucleiShrink nucleiExpand...
                waterMF
                continue
            end
            [cells, cytoplasm] = CycIFCellSeg_1cyc(FOVstack, nuclei,...
                nucleiExpand, waterMF);
            %nucleiBW = nuclei>=1; %convert final nuclei mask to black/white
            %cells = imdilate(nucleiBW,strel('sphere',100));
            %cells = waterMF.*uint16(nucleiExpand>=1);
            %cytoplasm = cells;
            %cytoplasm(nuclei>0) = 0;
            [bugs, bugsCellLabel] = CycIFBugSeg_1cyc(FOVstack, cells);
            featureData = CycIFData_1cyc(FOVstack,...
                nuclei, nucleiShrink, cells, cytoplasm, bugs, bugsCellLabel,...
                saveDirectory, name, FOV, channels, experiment, timepoint,...
                row, column, field);
%             nucleiEdge = edge(nuclei>0);
%             cellEdge = edge(cells>0);
%             bugEdge = edge(bugs>0);
%             cytoEdge = nucleiEdge + cellEdge;
%             allEdge = cytoEdge + bugEdge;
%             masks = cat(3, nuclei, nucleiShrink, nucleiExpand, cells,...
%                 cytoplasm, bugs, bugsCellLabel,...
%                 nucleiEdge, cellEdge, cytoEdge, bugEdge, allEdge);
%             edgeMasks = cat(3, nucleiEdge, cellEdge, cytoEdge, bugEdge, allEdge);
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
