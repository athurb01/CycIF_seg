% function to adentify foci of set size based on code from Clarence
% filterLoG function must be in path
function [fociMask] = CycIFFoci(image,nuclei, cells)
% convert nuclei and cells mask to logical then double
nuclei_BW = double(nuclei>0);
cells_BW = double(cells>0);
background = double(cells);
background(cells > 0) = 0;
background(cells == 0) = 1;


%apply LoG filter function with std of 1.5
I_log = filterLoG(image,1.5);

%find regional max
Imax= imregionalmax(I_log);

%sample nuclei spots to use as background,convert mean absolute deviation
%to std, and take 20 std above median to set threshold.
backSpots = I_log(background.*double(Imax)>0);
threshold=median(backSpots(:)) + 4*mad(backSpots(:))*1.4826;
fociMask=(I_log.*Imax)>threshold;