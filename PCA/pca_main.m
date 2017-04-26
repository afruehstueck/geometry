% @file     pca_main.m
% @author   afruehstueck
% @date     15/03/2017
%
% load some 2D and 3D pointcloud data from file and/or autogenerate random
% sample data
% apply some rotations and translations to the data to make it more
% interesting
% do principal component analysis on the data and visualize the results
% generate an oriented bounding box of the data in 2D and 3D on the basis 
% of the principal components of the data and visualize the bounding box

close all
clc;
clear;

%% PARAMETERS: select 2D or 3D data
%load 3D data from file
pcpath = 'teapot.ply';
pcpath = '../data/pointcloud/face.ply';
%pcpath = '../data/pointcloud/woman.ply';
%pcpath = '../data/pointcloud/bunny.ply';
%pcpath = '../data/pointcloud/lobster.ply';
%pcpath = '../data/pointcloud/pepper.ply';

dims = 3; %set to 2 or 3
randomRotation = true; %apply random rotation to data for more interesting results
loadFromFile = true; %true = load file, false = generate random data
N = 5000; %number of randomly generated data points

evaluatePCA(dims, randomRotation, loadFromFile, pcpath, N);

%load 2D data from file
csvpath = '../data/csv/freshman_kgs.csv';
csvpath = '../data/csv/heightweight_7500.csv';
csvpath = '../data/csv/snakes_count_1000.csv';

dims = 2; %set to 2 or 3
randomRotation = true; %apply random rotation to data for more interesting results
loadFromFile = false; %true = load file, false = generate random data
N = 5000; %number of randomly generated data points
evaluatePCA(dims, randomRotation, loadFromFile, csvpath, N);

function evaluatePCA(dims, randomRotation, loadFromFile, path, N)

%% load or generate data
if dims == 2 && loadFromFile
    %load data from file
    data = csvread(path);
elseif dims == 2
    %create randomized 2D data
    data = [random(makedist('Normal'), N, 1) random(makedist('Logistic'), N, 1)];
elseif dims == 3 && loadFromFile
    %load data from file
    pcloud = pcread(path);
    data = pcloud.Location;
    
    %subsample 3D data for big files > 15000pts
    % s = max(floor(length(data) / 15000), 1);
    % data = data(1:s:length(data), :);
elseif dims == 3
    %create randomized 3D data
    data = [random(makedist('Normal'), N, 1) random(makedist('Logistic'), N, 1) random(makedist('Weibull'), N, 1)];
else
    disp('Invalid number of dims');
    return;
end

[M, dims] = size(data); 

%rotate data points by random angle (to make things a bit more interesting)
if randomRotation && dims == 2
    randAngle = rand * (pi/2);
    disp(['Rotate by random angle ', num2str(randAngle)]);
    rot = matrix(2, 'rotation', randAngle);
    data = transformPoints(data', rot)';
elseif randomRotation && dims == 3 
    randAngleX = rand * pi;
    randAngleY = rand * pi;
    randAngleZ = rand * pi;
    disp(['Rotate around x by random angle ', num2str(randAngleX)]);
    disp(['Rotate around y by random angle ', num2str(randAngleY)]);
    disp(['Rotate around z by random angle ', num2str(randAngleZ)]);
    rot = matrix(3, 'rotation', 'x', randAngleX);
    data = transformPoints(data', rot)';
    rot = matrix(3, 'rotation', 'y', randAngleY);
    data = transformPoints(data', rot)';
    rot = matrix(3, 'rotation', 'z', randAngleZ);
    data = transformPoints(data', rot)';
end

disp([num2str(M), ' data values, ', num2str(dims), ' dimensions']);
% calculate the mean per dimension = column
data_mean = mean(data, 1);
disp(['Mean: (', num2str(data_mean), ')']);
%subtract the mean from all data points
data_centered = data - repmat(data_mean, M, 1);

%% PCA using COVARIANCE MATRIX
% calculate the covariance matrix 
covariance = (data_centered' * data_centered) / (M - 1);
% find the eigenvectors and eigenvalues 
[eVecs, eVals] = eig(covariance);
% extract eigenvalues along diagonal of matrix 
eVals = diag(eVals); 
% sort the variances in decreasing order 
[eVals, order] = sort(eVals, 'descend');
eVecs = eVecs(:, order);

% project the centralized original data set to new basis
data_PCAbasis = data_centered * eVecs;

figure1 = figure('Name', 'PCA Individual Steps', 'NumberTitle', 'off', 'Position', get(0, 'ScreenSize'));
if dims == 2 
    subplot(1,3,1);
    plot(data(:, 1), data(:, 2), '.')
    title('Original data');
    
    subplot(1,3,2);
    plot(data_centered(:, 1), data_centered(:, 2), '.')
    title('Centered data');
    
    subplot(1,3,3);
    plot(data_PCAbasis(:, 1), data_PCAbasis(:, 2), '.')
    title('Data in PCA basis');
elseif dims == 3
    subplot(1,3,1);
    pcshow(data);
    title('Original data');
    
    subplot(1,3,2);
    pcshow(data_centered);
    title('Centered data');
    
    subplot(1,3,3);
    pcshow(data_PCAbasis);
    title('Data in PCA basis');
end

%find minima, maxima and ranges
minima = min(data_PCAbasis, [], 1);
maxima = max(data_PCAbasis, [], 1);
ranges = maxima - minima;

%scale range by eigenvalue size, with biggest eigenvector corresponding to 
% biggest eigenvalue spanning the largest data range
max_range = max(ranges ./ 2);
ranges = repmat(max_range, dims, 1) .* ( eVals ./ max(eVals));

rangeVecs = ranges' .* eVecs;

%position labels close to arrows
textPos = (ranges' + 0.05 * max_range) .* eVecs;

plts = zeros(0, 2);
figure('Name', 'PCA', 'NumberTitle', 'off', 'Position', get(0, 'ScreenSize'));

plts(1) = subplot(1, 2, 1);
hold on;
if dims == 3
    %plot data points
    pcshow(data_centered);
    
    %display arrows for eigenvectors
    quiver3([0 0 0], [0 0 0], [0 0 0], rangeVecs(1, :), rangeVecs(2, :), rangeVecs(3, :), 'LineWidth', 2.5, 'Color', [0.27 0 0.58], 'AutoScale', 'off');

    %add text labels to eigenvectors
    text(textPos(1, :), textPos(2, :), textPos(3, :), {'E1', 'E2', 'E3'}, 'FontSize', 13)
    
    %get bounding box coordinates in PCA basis
    box_PCAbasis = zeros(8, 3);
    box_PCAbasis(1, :) = [ minima(1), minima(2), minima(3) ];
    box_PCAbasis(2, :) = [ maxima(1), minima(2), minima(3) ];
    box_PCAbasis(3, :) = [ maxima(1), maxima(2), minima(3) ];
    box_PCAbasis(4, :) = [ minima(1), maxima(2), minima(3) ];
    box_PCAbasis(5, :) = [ minima(1), minima(2), maxima(3) ];
    box_PCAbasis(6, :) = [ maxima(1), minima(2), maxima(3) ];
    box_PCAbasis(7, :) = [ maxima(1), maxima(2), maxima(3) ];
    box_PCAbasis(8, :) = [ minima(1), maxima(2), maxima(3) ];

    %transform bounding box coordinates back to original basis
    box = box_PCAbasis * eVecs';

    %face list for bounding box
    faces = [ 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8; 1 2 3 4; 5 6 7 8 ];
    
    %display semitransparent bounding box
    patch('Vertices', box, 'Faces', faces, 'FaceVertexCData', (1:6)', 'FaceColor', 'flat', 'FaceAlpha', 0.08);
elseif dims == 2 
    %plot data points
    plot(data_centered(:, 1), data_centered(:, 2), '.')
    
    %display arrows for eigenvectors
    quiver([0 0], [0 0], rangeVecs(1, :), rangeVecs(2, :), 'LineWidth', 2.5, 'Color', [0.27 0 0.58], 'AutoScale', 'off');
    
    %add text labels to eigenvectors
    text(textPos(1, :), textPos(2, :), {'E1', 'E2'}, 'FontSize', 13)
    
    %get bounding rectangle coordinates in PCA basis
    rect_PCAbasis = zeros(4, 2);
    rect_PCAbasis(1, :) = [ minima(1), minima(2) ];
    rect_PCAbasis(2, :) = [ maxima(1), minima(2) ];
    rect_PCAbasis(3, :) = [ maxima(1), maxima(2) ];
    rect_PCAbasis(4, :) = [ minima(1), maxima(2) ];
    
    %transform rectangle coordinates back to original basis
    rect = rect_PCAbasis * eVecs';
    
    %draw rectangle
    fill(rect(:,1), rect(:,2), [0.27 0 0.58], 'FaceAlpha', 0.08);
end
axis equal
title('Covariance');

%% PCA using SVD
[u, eVals, eVecs] = svd(data_centered, 0);
eVals = diag(eVals);
[eVals, order] = sort(eVals, 'descend');
eVecs = eVecs(:, order);
plts(2) = subplot(1, 2, 2);
hold on;
if dims == 3
    pcshow(data_centered);
    
    %display arrows for eigenvectors
    quiver3([0 0 0], [0 0 0], [0 0 0], eVecs(1, :), eVecs(2, :), eVecs(3, :), 'LineWidth', 2.5, 'Color', [0.27 0 0.58], 'AutoScale', 'off');
elseif dims == 2 
    plot(data_centered(:, 1), data_centered(:, 2), '.')
    
    %display arrows for eigenvectors
    quiver([0 0], [0 0], eVecs(1, :), eVecs(2, :), 'LineWidth', 2.5, 'Color', [0.27 0 0.58], 'AutoScale', 'off');
end
axis equal
title('SVD');

%link rotation of plts
link = linkprop(plts, {'CameraPosition', 'CameraUpVector'} );
rotate3d on

end