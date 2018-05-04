function [imageStack] = inputTiffStack(directory, channels, row, column, field)
% function to input all images from a Tiff stack file to a single array

cd (directory);

for k = 1:4
    a = [row, '-', column];
    b = ['(fld', field, 'wv', channels(k)];
    c = strjoin(a);
    d = strjoin(b);
    imName = strcat(c,d,').tif');
    currentImage = imread(char(imName));
    imageStack(:,:,k) = currentImage;
end