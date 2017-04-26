% @file     morphology_main.m
% @author   afruehstueck
% @date     14/03/2017
%
% load image from file
% apply various morphological operations to image

close all
clc;
clear;

%impath = '../data/img/beaker.jpg';
%impath = '../data/img/squares.jpg';
%impath = '../data/img/beta.jpg';
%impath = '../data/img/bunnybw.jpg';
impath = '../data/img/cat.jpg';
%impath = '../data/img/turing.jpg';
%impath = '../data/img/zebra.jpg';

shapepath_3 = '../data/img/shape_3.jpg';
shapepath_8 = '../data/img/shape_8.jpg';

img = imread(impath);
%convert image to black and white
BWimg = im2bw(img, 0.5);

[img_y, img_x] = size(img(:, :, 1));
se = [];
se_idx = 1;
se_sz = floor(img_x / 20);

%define a series of structure elements
se = cat(3, se, strel('line', 2*se_sz, 90));
se = cat(3, se, strel('square', 2*se_sz));
se = cat(3, se, strel('disk', se_sz));
se = cat(3, se, strel('diamond', se_sz));
se = cat(3, se, strel('rectangle', [se_sz, 2*se_sz]));

shape = zeros(2*se_sz,2*se_sz);
shape(1:5, 1:5) = ones(5, 5);
shape(1:5, end-4:end) = ones(5, 5);
shape(end-4:end, se_sz-2:se_sz+2) = ones(5, 5);
se = cat(3, se, strel('arbitrary', shape));

%load strel from image
shape_8 = im2bw(imread(shapepath_8), 0.5);
se = cat(3, se, strel('arbitrary', shape_8));

shape_3 = im2bw(imread(shapepath_3), 0.5);
se = cat(3, se, strel('arbitrary', shape_3));

n_se = length(se);
w = 6;

scr = get(0, 'ScreenSize'); 
morph_fig = figure('Name', 'Morphological Operations', 'NumberTitle', 'off', 'Position', scr);
%set dark background color
set(gcf,'Color', [0.1  0.1  0.1])

%show original image to the left
subplot(n_se, w, 1)
imshow(BWimg)
title('ORIGINAL', 'Color', 'w')
for i=1:n_se    
    cur_se = se(:, :, i);
    
    %show strel
    ne = subplot(n_se, w, (i - 1) * w + 2);
    pos = get(ne, 'position');
    imshow(cur_se.Neighborhood) 
    [se_y, se_x] = size(cur_se.Neighborhood);
    strsz = (se_x / img_x);
    %rescale strel figure to correspond to comparative 'real' size
    set(ne, 'position', (pos.*[1 1 strsz strsz])-[pos(3)*(strsz - 1)/2 pos(4)*(strsz - 1)/2 0 0])
    
    %erode
    erodedBW = imerode(BWimg, cur_se);
    subplot(n_se, w, (i - 1) * w + 3)
    imshow(erodedBW)
    
    %dilate
    dilatedBW = imdilate(BWimg, cur_se);
    subplot(n_se, w, (i - 1) * w + 4)
    imshow(dilatedBW)
       
    %open
    subplot(n_se, w, (i - 1) * w + 5)
    imshow(imdilate(erodedBW, cur_se))
    
    %close
    subplot(n_se, w, (i - 1) * w + 6)
    imshow(imerode(dilatedBW, cur_se))
    
end

%set figure titles in top row
subplot(n_se, w, 3)
title('ERODE', 'Color', 'w')
subplot(n_se, w, 4)
title('DILATE', 'Color', 'w')
subplot(n_se, w, 5)
title('OPEN', 'Color', 'w')
subplot(n_se, w, 6)
title('CLOSE', 'Color', 'w')

subplotsqueeze(morph_fig, 1.35);