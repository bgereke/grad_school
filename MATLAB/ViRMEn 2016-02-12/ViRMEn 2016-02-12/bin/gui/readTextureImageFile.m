function [img, errorString] = readTextureImageFile(filename)

info = imfinfo(filename);
errorString = '';

switch info.ColorType
    case 'grayscale'
        img = imread(filename);
        img = cat(3,img,img,img);
    case 'truecolor'
        img = imread(filename);
    case 'indexed'
        [img, colorMap] = imread(filename);
        img = ind2rgb(img,colorMap);
    otherwise
        img = [];
        errorString = 'Could not identify file format.';
end

if strcmp(class(img),'uint8')
    img = double(img)/255;
end