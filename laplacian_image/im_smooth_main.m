% @file     im_smooth_main.m
% @author   anna fruehstueck
% @date     01/03/2017
close all
clc;
clear;

%impath = '../data/img/bunny.jpg';
impath = '../data/img/sloth.jpg';
%impath = '../data/img/flower.jpg';

%impath = '../data/img/stripes.png';
%impath = '../data/img/zebra.jpg';
%impath = '../data/img/stripes_diagonal.jpg';
%impath = '../data/img/loops.jpg';
img = imread(impath);

laplace1  = [0 1 0; 1 -4 1; 0 1 0];
laplace2  = [1 1 1; 1 -8 1; 1 1 1];

doRGB = true; %filter in rgb or grayscale
lambda = 0.25; %limit to break the filtering at 0.25 | <0.01, change visually not really noticeable
iterations = 100;
filtertype = string('filter2'); %'matrix', 'filter2' or 'circshift'
filter = laplace1;

if size(img, 3) == 3 && ~doRGB
    disp('convert image to grayscale...');
    img = rgb2gray(img); %convert img to grayscale
end

if isa(img, 'integer')
    disp('convert image to double values...');
    img = im2double(img);
end

scr = get(0, 'ScreenSize'); 

[sy, sx] = size(img(:, :, 1));
num_px = sy * sx;
L = laplacian_matrix(sx, sy);

symm = issymmetric(L);
if symm
    msg = 'symmetric';
else
    msg = 'not symmetric';
end

disp(['Laplacian matrix is ', msg]);

L_eval = laplacian_matrix(30, 30);
condition_number = condest(L_eval);
disp(['Condition number for ', num2str(30), 'x', num2str(30), ' laplacian matrix: ', num2str(condition_number)]);
determinant = det(L_eval);
%disp(['Determinant: ', num2str(determinant)]);
[V, D] = eig(full(L_eval));
plot(1:length(D), diag(D))
title('900x900 laplacian matrix')
xlabel('index')
ylabel('eigenvalues')

figure;
spy(L_eval);
title('sparse 900x900 laplacian matrix');

condition_number = condest(L);
disp(['Condition number for ', num2str(num_px), 'x', num2str(num_px), ' laplacian matrix: ', num2str(condition_number)]);

figure('Name', 'Laplacian Smoothing', 'NumberTitle', 'off', 'Position', [scr(3)/4 50 scr(3)/2 scr(3)/4]);

subplot(1,2,1)
imshow(img)
title('Original')
tic;

for t = 1:iterations
    img_filtered = [];
    for c = 1:size(img, 3)
        channel = img(:, :, c);
        switch filtertype
            case 'filter2'
                tmp_img = padarray(channel, [1 1], 'replicate');
                channel_filtered = filter2(filter, tmp_img, 'valid');    
            case 'circshift'
                channel_filtered = zeros(size(channel));
                d = ceil(size(filter)./2);
                for fy = 1:size(filter, 1)
                    for fx = 1:size(filter, 2)
                        f = filter(fy, fx);
                        if f ~= 0
                            channel_filtered = channel_filtered + f * circshift(channel, [fy - d(1), fx - d(2)]);
                        end
                    end
                end
            case 'matrix'    
                img_vec = reshape(channel, [num_px, 1]);
                img_vec = L * img_vec;
                channel_filtered = reshape(img_vec, [sy, sx]);
        end       
        img_filtered = cat(3, img_filtered, channel_filtered);
    end
    subplot(1,2,2)
    img = img + lambda * img_filtered;
    imshow(img)
    title(['\lambda= ', num2str(lambda), ' iteration ' , num2str(t), '/', num2str(iterations)])
    drawnow;
end
toc;