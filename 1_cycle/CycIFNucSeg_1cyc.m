function [nuclei, nucleiShrink, nucleiExpand, waterMF] = CycIFNucSeg(FOVstack)
%% nuclei segmentation
%based on Clarence Yapp code
%manually set cycle to segment nuclei
I = FOVstack(:,:, 1);
I_close=imclose(I,strel('sphere',4)); %dilation then erosion, blurring effect
%Ith = imtophat(I,strel('disk',15));

%fine tune initial nuclei mask
bw= I_close > 2000; %manual threshold
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

nuclei = waterMF.*uint16(imgLabel>=1); %gives nuclei label from tesselation

%shrink nuclei and expand nuclei to use for measurementsso cyt and nuc do
% not overlap
nucleiBW = nuclei>=1; %convert final nuclei mask to black/white
nucleiShrink = imerode(nucleiBW,strel('sphere',3));
nucleiShrink = waterMF.*uint16(nucleiShrink>=1); %gives label # from tesselation
nucleiExpand = imdilate(nucleiBW,strel('sphere',3));
nucleiExpand = waterMF.*uint16(nucleiExpand>=1);