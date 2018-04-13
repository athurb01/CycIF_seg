function [cells, cytoplasm] = CycIFCellSeg(FOVstack, nuclei, nucleiExpand, waterMF)
%% cytoplasm segmentation
% create sum of multiple channel image to use for cell segmentation
% NEED TO CHANGE SO IT IS DEPENDENT ON CHANNEL NAME NOT POSITION
sum = uint16(zeros(2048,2048));
for c = 1:4
    sum = sum + FOVstack(:,:,c);
end 
 
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
cells = uint16(zeros(2048,2048));
bwLabel = bwlabel(bw);
numObjects = max(bwLabel(:));
for obj = 1:numObjects
    nucleiOverlap = unique(nuclei(bwLabel == obj));
    if length(nucleiOverlap) == 1
        cells(bwLabel == obj) = 0;
    elseif length(nucleiOverlap) == 2
        cells(bwLabel == obj) = max(nucleiOverlap);
    else
        for i = 1:2048
            for j = 1:2048
                if bwLabel(i,j) == obj
                    cells(i,j) = waterMF(i,j);
                end
            end
        end
        
    end
end

bw2 = bwareaopen(cells, 1000); %get rid of small pieces of cells created by waterMF
cells = cells.*uint16(bw2); %give the cells label back
cytoplasm = cells;
cytoplasm(nucleiExpand>0) = 0;

