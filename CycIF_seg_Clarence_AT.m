%% nuclei segmentation
root = 'C:\Users\Amy Thurber\Dropbox (Partners HealthCare)\FI07_for_Clarence\rd3\';  
I=imread([root 'A - 05(fld 2 wv UV - DAPI).tif']);
   
I_close=imclose(I,strel('sphere',2)); %has blurring effect
Ith = imtophat(I,strel('disk',15));

%fine tune initial nuclei mask
bw=Ith>4000;
bw= imclose(bw,strel('sphere',2));
bw = imfill(bw,'holes');
Ierode=bwareaopen(bw,75);

IdistTF =imcomplement(bwdist(~Ierode));
dGauss = imhmin(IdistTF,0.5);
dGauss = imgaussfilt3(dGauss,2);
Imax = imregionalmin(dGauss);

bw=Ierode;
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
waterMF=watershed(markers);
nuclei = waterMF.*uint8(imgLabel>1);

figure,imshowpair(edge(nuclei>0),I)

%% cytoplasm segmentation
FITC=imread([root 'A - 05(fld 2 wv Blue - FITC).tif']);
Cy3=imread([root 'A - 05(fld 2 wv Green - dsRed).tif']);
Cy5=imread([root 'A - 05(fld 2 wv Red - Cy5).tif']);
sum =FITC+Cy3+Cy5;
Ith = imtophat(sum,strel('disk',30));
Ibh = imbothat(sum,strel('disk',30));
I_bg = Ith-Ibh;

thresh=adaptthresh(I_bg,0.9);
bw=imbinarize(I_bg,thresh);
bw=imopen(bw,strel('sphere',2));
bw = imclose(bw, strel('sphere', 8));
bw(nuclei > 0) = 1; %include all of nuclei mask in cell mask
bw=bwareaopen(bw,150);
bw=imfill(bw,'holes');


 cells=(waterMF .* uint8(bw)).*uint8(bwareaopen(waterMF .* uint8(bw),200));
 cellEdge=edge((cells>0),'Sobel');


nucleiEdge=edge(nuclei>0);

allEdge= cellEdge + nucleiEdge;
figure,imshowpair(allEdge,sum)

%% statistics
