% @file     svd_main.m
% @author   anna fruehstueck
% @date     20/02/2017
%
% load image from file, decompose using SVD and reconstructs image using
% varying numbers of singular values. evaluate the error of low-rank
% approximation

close all
clc;
clear;

impath = '../data/img/bunny.jpg';
%impath = '../data/img/stripes.png';
%impath = '../data/img/zebra.jpg';
%impath = '../data/img/flower.jpg';
%impath = '../data/img/carpet.jpg';
%impath = '../data/img/stripes_diagonal.jpg';
%impath = '../data/img/loops.jpg';
img = imread(impath);

doRGB = 1;

if size(img, 3) == 3 && ~doRGB
    disp('convert image to grayscale');
    img = rgb2gray(img); %convert img to grayscale
end
%imshow(A);

if isa(img, 'integer')
    disp('convert image to double values');
    img = im2double(img);
end

U = []; S = []; V = [];

%do SVD for all current channels (either grayscale or RGB separately
for channel = 1:size(img, 3)
    [u,s,v] = svd(img(:,:,channel)); 
    U = cat(3, U, u);
    S = cat(3, S, s);
    V = cat(3, V, v);
end

figure('Name', 'low-rank approximation for N singular values', 'NumberTitle', 'off', 'Position', get(0, 'ScreenSize'));

t = [1, 3, 5, 10, 25, 50, 100, 200]; %numbers of used singular values 
w = 4; h = ceil(length(t) / 4); %number of subplots
meanErrs = [];
nVals = [];
[ha, pos] = tight_subplot(h*2, w, .01, .05, .02); %using tight_subplot because subfigures leave so much unused space
for idx = 1:length(t)
    reconstruction = [];
    N = t(idx);
    for channel = 1:size(img, 3)
        C=S(:, :, channel); 
        C(N+1:end, N+1:end) = 0; %set all diagonal values above N to zero
        
        D = U(:, :, channel) * C * V(:, :, channel)'; %reconstruct low-rank approximation
        reconstruction = cat(3, reconstruction, D);  
    end
    
    row = floor((idx - 1) / w); %find out current row in plot
    
    axes(ha(idx + row * w)); %position in subplot
    imshow(reconstruction);
    title(sprintf('N = %d', N));  
    
    %calculate mean squared error
    error=sum(sum(sum((img - reconstruction).^2)));
    meanErrs = [meanErrs; error];
    nVals = [nVals; N];
    
    axes(ha(idx + (row + 1) * w)); %position in subplot
    
    %imshow(img - reconstruction);
    imshow(1 - (img - reconstruction)); %showing the inverse since it is much easier to see (my subjective opinion)
end

% dislay the error graph
figure('Name', 'Error graph', 'NumberTitle', 'off'); 
plot(nVals, meanErrs);
title('Mean squared error in low-rank approximation');
grid on
xlabel('N');
ylabel('error');
