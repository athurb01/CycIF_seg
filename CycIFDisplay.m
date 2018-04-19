%% random bits of code to display CycIF data and images in useful ways

%% input saved data to featureData (same name as master function)
% used imiport tool to generate code, note that well data loses letters
% not currently working
directoryName = 'C:\Users\Amy Thurber\Dropbox (Partners HealthCare)\Experiments\FI09\CycIF_Test_images_FI09\DataOutput\';
filename = strcat(directoryName, 'FI09_03h_20x_B05_fld2.txt');
delimiter = ',';
startRow = 2;

% Read columns of data as text:
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';

% Open the text file.
fileID = fopen(filename,'r');

% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');

% Close the text file.
fclose(fileID);

% Convert the contents of columns containing numeric text to numbers.
% Replace non-numeric text with NaN.
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

for col=[1,2,3,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26]
    % Converts text in the input cell array to numbers. Replaced non-numeric
    % text with NaN.
    rawData = dataArray{col};
    for row=1:size(rawData, 1)
        % Create a regular expression to detect and remove non-numeric prefixes and
        % suffixes.
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData(row), regexstr, 'names');
            numbers = result.numbers;
            
            % Detected commas in non-thousand locations.
            invalidThousandsSeparator = false;
            if numbers.contains(',')
                thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(numbers, thousandsRegExp, 'once'))
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % Convert numeric text to numbers.
            if ~invalidThousandsSeparator
                numbers = textscan(char(strrep(numbers, ',', '')), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch
            raw{row, col} = rawData{row};
        end
    end
end


% Split data into numeric and string columns.
rawNumericColumns = raw(:, [1,2,3,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26]);
rawStringColumns = string(raw(:, 4));


% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),rawNumericColumns); % Find non-numeric cells
rawNumericColumns(R) = {NaN}; % Replace non-numeric cells

% Make sure any text containing <undefined> is properly converted to an <undefined> categorical
idx = (rawStringColumns(:, 1) == "<undefined>");
rawStringColumns(idx, 1) = "";

% Create output variable
featureData = table;
featureData.experiment = cell2mat(rawNumericColumns(:, 1));
featureData.timepoint = cell2mat(rawNumericColumns(:, 2));
featureData.well = cell2mat(rawNumericColumns(:, 3));
featureData.row = categorical(rawStringColumns(:, 1));
featureData.column = cell2mat(rawNumericColumns(:, 4));
featureData.field = cell2mat(rawNumericColumns(:, 5));
featureData.object = cell2mat(rawNumericColumns(:, 6));
featureData.channel = cell2mat(rawNumericColumns(:, 7));
featureData.meanNuclei = cell2mat(rawNumericColumns(:, 8));
featureData.modeNuclei = cell2mat(rawNumericColumns(:, 9));
featureData.medianNuclei = cell2mat(rawNumericColumns(:, 10));
featureData.intIntensityNuclei = cell2mat(rawNumericColumns(:, 11));
featureData.meanCell = cell2mat(rawNumericColumns(:, 12));
featureData.modeCell = cell2mat(rawNumericColumns(:, 13));
featureData.medianCell = cell2mat(rawNumericColumns(:, 14));
featureData.intIntensityCell = cell2mat(rawNumericColumns(:, 15));
featureData.meanCytoplasm = cell2mat(rawNumericColumns(:, 16));
featureData.modeCytoplasm = cell2mat(rawNumericColumns(:, 17));
featureData.medianCytoplasm = cell2mat(rawNumericColumns(:, 18));
featureData.intIntensityCytoplasm = cell2mat(rawNumericColumns(:, 19));
featureData.meanNC = cell2mat(rawNumericColumns(:, 20));
featureData.modeNC = cell2mat(rawNumericColumns(:, 21));
featureData.medianNC = cell2mat(rawNumericColumns(:, 22));
featureData.intIntensityNC = cell2mat(rawNumericColumns(:, 23));
featureData.meanBug = cell2mat(rawNumericColumns(:, 24));
featureData.intIntensityBug = cell2mat(rawNumericColumns(:, 25));

% Clear temporary variables
clearvars filename delimiter startRow formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp rawNumericColumns rawStringColumns R idx;


%% import data take 2
directoryName = 'C:\Users\Amy Thurber\Dropbox (Partners HealthCare)\Experiments\FI09\CycIF_Test_images_FI09\DataOutput\';
filename = strcat(directoryName, 'FI09_03h_20x_B05_fld2.txt');

table = importdata(filename);


%% select data for one feature one channel
numCells = max(nuclei(:));
channelName = 'RelA';
featureName = 'modeNuclei';
d=1;
data = zeros(numCells, 1);
for i = 1:length(featureData);
    if strcmp(featureData(i).channel, channelName)
        data(d) = featureData(i).(featureName);
        d = d+1;
    end
end

% show data as histogram
figure;
histogram(data, 20)

% scatterplot vs object number
objects = (1:numCells);
figure; scatter(objects, data, 9, objects, 'filled' )

%% select top and bootom 5% expressing cells
N = numCells*.1;
data(isnan(data))=0; %set all NaN to 0
[value, ix] = sort(data, 'descend');
for i = 1:N
    highCells(i) = ix(i);
    lowCells(i) = ix(numCells-(i-1));
end

% display adjusted channel image
slice = find(strcmp(channelName, channelNames));
sliceAdj = imadjust(FOVstack(:,:, slice), [0 .1]);
figure; hold on;
imshowpair(sliceAdj, cellEdge)
% display multiple bounding boxes
for k = 1 : N
    thisBB = cellStats(highCells(k)).BoundingBox;
    rectangle('Position', [thisBB(1),thisBB(2),thisBB(3),thisBB(4)],...
    'EdgeColor',[[1,1,0.4]],'LineWidth',1 )
    thisBB = cellStats(lowCells(k)).BoundingBox;
    rectangle('Position', [thisBB(1),thisBB(2),thisBB(3),thisBB(4)],...
    'EdgeColor',[0,0.5,0.4],'LineWidth',1 )
end

%% display image showing nuclei centroid label
allMaskRGB = cat(3, edgeMasks(:,:,4), edgeMasks(:,:,2), edgeMasks(:,:,1));
figure; imshow(allMaskRGB)
hold on
for k = 1:numCells
    c = nucleiStats(k).Centroid;
    text(c(1), c(2), sprintf('%d', k), ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle');
end
hold off

%% display bug ratio as heatmap
% need to make bug mask more stringent or most of data is really background
B=FOVstack(:,:, 2+(maxCycle*2));
bugsBW=B>4000;
bugs = bwlabel(bugsBW);
bugsCellLabel = cells.*uint16(bugsBW);
bugEdge = edge(bugs>0);

bugRatio = double(FOVstack(:,:,9))./double(FOVstack(:,:,16));
bugRatio(bugs == 0) = 0; % set background to 0
B05f2_24h = bugRatio % save array with unique name

%display as heatmap
clims = [0 0.4];

figure; hold on;
cmap = colormap;
cmap_mod = cmap;
cmap_mod(1,:) = [1 1 1];
colormap(cmap_mod)
imagesc(bugRatio, clims);
colorbar
xlim([0 2048])
ylim([0 2048])
