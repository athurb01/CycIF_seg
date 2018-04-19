%% nuclei segmentation
root = 'Z:\IDAC\Clarence\LSP\Zenon test 2 layers\';  
I=imread([root 'C - 10(fld 3 wv UV - DAPI).tif']);
   
I_close=imclose(I,strel('sphere',2));
Ith = imtophat(I,strel('disk',15));

%fine tune initial nuclei mask
bw=Ith>4000;
bw= imclose(bw,strel('sphere',2));
bw = imfill(bw,'holes');
Ierode=bwareaopen(bw,100);

IdistTF =imcomplement(bwdist(~Ierode));
dGauss = imhmin(IdistTF,0.5);
dGauss = imgaussfilt3(dGauss,2);
Imax = imregionalmin(dGauss);

bw=Ierode;
imgDist=-bwdist(~bw);
imgDist=imimposemin(imgDist,Imax);
imgDist(~bw)=-inf;
imgLabel=watershed(imgDist);

%tesselation
markers =bwdist(imgLabel>1);
waterMF=watershed(markers);
nuclei = waterMF.*uint8(imgLabel>1);

figure,imshowpair(edge(nuclei>0),I)

%% cytoplasm segmentation
FITC=imread([root 'C - 10(fld 3 wv Blue - FITC).tif']);
Cy3=imread([root 'C - 10(fld 3 wv Green - dsRed).tif']);
Cy5=imread([root 'C - 10(fld 3 wv red - Cy5).tif']);
sum =FITC+Cy3+Cy5;
Ith = imtophat(sum,strel('disk',30));
Ibh = imbothat(sum,strel('disk',30));
I_bg = Ith-Ibh;

thresh=adaptthresh(I_bg,0.9);
bw=imbinarize(I_bg,thresh);
bw=imopen(bw,strel('sphere',2));
bw=bwareaopen(bw,200);
bw=imfill(bw,'holes');


 cells=(waterMF .* uint8(bw)).*uint8(bwareaopen(waterMF .* uint8(bw),200));
 cellEdge=edge(cells,'Sobel');


nucleiEdge=edge(nuclei,'Sobel');

allEdge= cellEdge + nucleiEdge;
imshowpair(allEdge,sqrt(normalize(double(I_bg))))