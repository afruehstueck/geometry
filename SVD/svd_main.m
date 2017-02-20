% @file     svd_main.m
% @author   anna fruehstueck
% @date     20/02/2017

impath = '../data/bunny.jpg';
img = imread(impath);

if size(img, 3) == 3
    disp('convert image to grayscale');
    img = rgb2gray(img); %convert img to grayscale
end
%imshow(A);

if isa(img, 'integer')
    disp('convert image to double values');
    img = double(img); %convert img to double values
end


[U,S,V] = svd(img); 