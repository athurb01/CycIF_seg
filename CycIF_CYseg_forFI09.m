%% import registered images
root = 'C:\Users\Amy Thurber\Dropbox (Partners HealthCare)\Experiments\FI09\CycIF_Test_images_FI09\';
experiment = 'FI09_';
timepoint = '3h_';
mag = '20x_';
well = 'D05';
field = '2';

cd (root);
fname = strcat(experiment, timepoint, mag, well, '_fld', field, '.tif');
info = imfinfo(fname);
imageStack = [];
imageStack = uint16(imageStack);
numberOfImages = length(info);
for k = 1:numberOfImages
    currentImage = imread(fname, k, 'Info', info);
    imageStack(:,:,k) = currentImage;
end

%% nuclei segmentation
% root = 'C:\Users\Amy Thurber\Dropbox (Partners HealthCare)\CycIF_scripts\CycIF_Test_images_FI09\';  
% well = 'B - 05';
% I=imread([root strcat(well, '(fld 2 wv UV - DAPI).tif')]);
maxCycle = 7;
I = FOVstack(:,:, maxCycle+1);
I_close=imclose(I,strel('sphere',4)); %dilation then erosion, blurring effect
%Ith = imtophat(I,strel('disk',15));

%fine tune initial nuclei mask
bw= I_close > 8000; %manual threshold
bw= imclose(bw,strel('sphere',5));
bw = imfill(bw,'holes');
bw = bwareaopen(bw,150); %remove objects with less than x pixels


IdistTF =imcomplement(bwdist(~bw));
dGauss = imhmin(IdistTF,0.5);
dGauss = imgaussfilt3(dGauss,4);
Imax = imregionalmin(dGauss);


imgDist=-bwdist(~bw);
imgDist=imimposemin(imgDist,Imax);
imgDist(~bw)=-inf;
imgWater=watershed(imgDist);
imgMask = zeros(size(imgWater));
imgMask(imgWater>1) = 1;
imgMask = bwareaopen(imgMask,75); %get rid of small nuc made by watershed
imgLabel = bwlabel(imgMask);

%tesselation
markers =bwdist(imgLabel>1);
waterMF=watershed(markers); %assigns every pixel to closest nuclei
waterMF = uint16(waterMF); 
% write if statement to convert waterMF to uint16 if originally uint8
nuclei = waterMF.*uint16(imgLabel>=1); %gives nuclei label from tesselation

%shrink nuclei and expand nuclei to use for measurementsso cyt and nuc do
% not overlap
% **add code to make shrink and expand have the same labels as nuclei
nucleiBW = nuclei>=1; %convert final nuclei mask to black/white
nucleiShrink = imerode(nucleiBW,strel('sphere',3));
nucleiShrink = waterMF.*uint16(nucleiShrink>=1); %gives label # from tesselation
nucleiExpand = imdilate(nucleiBW,strel('sphere',3));
nucleiExpand = waterMF.*uint16(nucleiExpand>=1);

figure,imshowpair(edge(nuclei>0),I)

% % Create masked image.
% maskedImage = I;
% maskedImage(~nuclei) = 0;

%% cytoplasm segmentation
sum = FOVstack(:,:,23) + FOVstack(:,:,27);%change to select markers
% currently mitotracker red, RelA 
 
Ith = imtophat(sum,strel('disk',30));
Ibh = imbothat(sum,strel('disk',30));
I_bg = Ith-Ibh;

thresh=adaptthresh(I_bg,0.9);
bw=imbinarize(I_bg,thresh);
bw=imopen(bw,strel('sphere',2));
bw = imclose(bw, strel('sphere', 8));
bw(nucleiExpand > 0) = 1; %include all of nuclei mask in cell mask
bw=bwareaopen(bw,500);
bw=imfill(bw,'holes');
% can you write if statement to get rid of cells that don't overlap with
% nuclei?
%or can you re-write segmentation so objects with only 1 nuclei are not
%split?

 cells =(waterMF .* uint16(bw)).*uint16(bwareaopen(waterMF .* uint16(bw),200));
 cellEdge=edge((cells>0),'Sobel');

 cytoplasm = cells;
 cytoplasm(nucleiExpand>0) = 0;

nucleiEdge=edge(nuclei>0);

allEdge= cellEdge + nucleiEdge;
figure,imshowpair(allEdge,sum)

%% bug segmentation
%creates a mask of bug location but does not try to separate individual
%bugs
B=FOVstack(:,:, 2+(maxCycle*2));

bugs=B>1000;
bugs = bwlabel(bugs);
figure; imshowpair(edge(bugs>1), B)

%% extract data
%currently no background subtraction or other corrections
%all values extracted from all channels
channelNames = {'Brightfield', 'Hoechst1', 'Hoechst2', 'Hoechst3',...
    'Hoechst4','Hoechst5', 'Hoechst6', 'Hoechst7', 'pH_GFP', 'p_p38',...
    'p_TBK1', 'Bax', 'Nos2', 'p_cJun', 'none', 'mTB_mCherry', 'p_STAT1',...
    'p_ERK', 'none', 'CD68', 'none', 'actin', 'mitoRed', 'p_STAT3',...
    'p_STAT6', 'Bcl2', 'RelA', 'LC_3', 'Ki67'};

properties = {'Area', 'Centroid', 'PixelList'};
nucleiStats = regionprops(nuclei, I, properties);
cellStats = regionprops(cells, I, properties);
cytoplasmStats = regionprops(cytoplasm, I, properties);

numCells = max(nuclei(:));

for ch = maxCycle+1:(maxCycle*4)+1
    currentChannel = FOVstack(:,:,ch);
    for obj = 1:numCells
        meanNuclei(obj).(char(channelNames(ch))) = mean(currentChannel(nuclei == obj));
        modeNuclei(obj).(char(channelNames(ch))) = mode(currentChannel(nuclei == obj));
        medianNuclei(obj).(char(channelNames(ch))) = median(currentChannel(nuclei == obj));
        intIntensityNuclei(obj).(char(channelNames(ch))) =...
            meanNuclei(obj).(char(channelNames(ch))) * nucleiStats(obj).Area;
        meanCell(obj).(char(channelNames(ch))) = mean(currentChannel(cells == obj));
        modeCell(obj).(char(channelNames(ch))) = mode(currentChannel(cells == obj));
        medianCell(obj).(char(channelNames(ch))) = median(currentChannel(cells == obj));
        intIntensityCell(obj).(char(channelNames(ch))) =...
            meanCell(obj).(char(channelNames(ch))) * cellStats(obj).Area;
        meanCytoplasm(obj).(char(channelNames(ch))) = mean(currentChannel(cytoplasm == obj));
        modeCytoplasm(obj).(char(channelNames(ch))) = mode(currentChannel(cytoplasm == obj));
        medianCytoplasm(obj).(char(channelNames(ch))) = median(currentChannel(cytoplasm == obj));
        intIntensityCytoplasm(obj).(char(channelNames(ch))) =...
            meanCytoplasm(obj).(char(channelNames(ch))) * cytoplasmStats(obj).Area;
        meanNC(obj).(char(channelNames(ch))) = meanNuclei(obj).(char(channelNames(ch)))...
            /meanCytoplasm(obj).(char(channelNames(ch)));
        modeNC(obj).(char(channelNames(ch))) = modeNuclei(obj).(char(channelNames(ch)))...
            /modeCytoplasm(obj).(char(channelNames(ch)));
        medianNC(obj).(char(channelNames(ch))) = medianNuclei(obj).(char(channelNames(ch)))...
            /medianCytoplasm(obj).(char(channelNames(ch)));
        intIntensityNC(obj).(char(channelNames(ch))) = intIntensityNuclei(obj).(char(channelNames(ch)))...
            /intIntensityCytoplasm(obj).(char(channelNames(ch)));
    end
end


%% data visualization

% histogram of nuclei pixel intensity from 12 random cells
randomCells = randi([1 numCells],1,12);
edges = linspace(0, 20000, 90);

figure;
for p = 1:12
    subplot(3,4,p)
    hold on
    histogram(cellStats(randomCells(p)).nucleiPixelsCy5, edges)
    plot([cellStats(randomCells(p)).nucleiMean, cellStats(randomCells(p)).nucleiMean],ylim,'r--','LineWidth',2)
    title(strcat('cell ', num2str(randomCells(p))))
end

% histogram of population mean, mode, and median
figure;
hold on
subplot(2,3,1)
histogram(nucleiMeans); title('nuclei means');
subplot(2,3,2)
histogram(nucleiModes); title('nuclei modes');
subplot(2,3,3)
histogram(nucleiMedians); title('nuclei medians');
subplot(2,3,4)
histogram(ncMeansRatio); title('n/c means ratio');
subplot(2,3,5)
histogram(ncModesRatio); title('n/c modes ratio');
subplot(2,3,6)
histogram(ncMediansRatio); title('n/c medians ratio');



    