function [imageStack] = inputTiffStack(directory, experiment, timepoint, mag, maxCycle, row, column, field)
% function to input all images from a Tiff stack file to a single array

cd (directory);
fname = char(strcat(experiment, timepoint, mag, '_', row, column, '_fld', field, '.tif'));
info = imfinfo(fname);
imageStack = [];
imageStack = uint16(imageStack);
numberOfImages = length(info);
for k = 1:numberOfImages
    currentImage = imread(fname, k, 'Info', info);
    imageStack(:,:,k) = currentImage;
end