% @file     pca_main.m
% @author   afruehstueck
% @date     15/03/2017

close all
clc;
clear;

%load 3D data from file
pcpath = 'teapot.ply';
%pcpath = '../data/pointcloud/face.ply';
pcpath = '../data/pointcloud/woman.ply';
%pcpath = '../data/pointcloud/bunny.ply';
%pcpath = '../data/pointcloud/lobster.ply';
%pcpath = '../data/pointcloud/pepper.ply';
pcloud = pcread(pcpath);
data3D = pcloud.Location;
% max_pts = max(data3D, [], 1);
% min_pts = min(data3D, [], 1);
% range_pts = max_pts - min_pts;
% 
% %normalize vertex values to [-1, 1] range for any plot
% for dim = 1:3
%     data3D(:, dim) = 2*(data3D(:, dim) - min_pts(dim)) / range_pts(dim) - 1;
% end

%subsample 3D data for big files > 5000pts
s = max(floor(length(data3D) / 15000), 1);
data3D = data3D(1:s:length(data3D), :);

%load 2D data from file
csvpath = '../data/csv/freshman_kgs.csv';
csvpath = '../data/csv/heightweight_200.csv';
csvpath = '../data/csv/heightweight_7500.csv';
csvpath = '../data/csv/snakes_count_1000.csv';
data2D = csvread(csvpath);

%create random 2D data
data2D = [random(makedist('Normal'), 1000, 1) random(makedist('Logistic'), 1000, 1)];

%% select 2D or 3D data
data = data2D;
[M, dims] = size(data); 

%rotate data points by random angle (to make it more interesting)
if dims == 2
    randAngle = rand * (pi/2);
    disp(['Rotate by random angle ', num2str(randAngle)]);
    rot = matrix(2, 'rotation', randAngle);
    data = transformPoints(data', rot)';
elseif dims == 3 
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
data_centered = data - repmat(data_mean, M, 1);%repmat(data_mean, M, 1); 

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
% project the original data set 

data_PCAbasis = (data_centered * eVecs);

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

figure('Name', 'PCA', 'NumberTitle', 'off', 'Position', get(0, 'ScreenSize'));
subplot(1, 2, 1);
hold on;
if dims == 3
    pcshow(data_centered);
    labels = {'E1', 'E2', 'E3'};
    
    %display arrows for eigenvectors
    quiver3([0 0 0], [0 0 0], [0 0 0], rangeVecs(1, :), rangeVecs(2, :), rangeVecs(3, :), 'LineWidth', 2.5, 'Color', [0.27 0 0.58], 'AutoScale', 'off');

    %add text labels to eigenvectors
    text(textPos(1, :), textPos(2, :), textPos(3, :), labels, 'FontSize', 13)
    
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

    %transform back to original basis
    box = box_PCAbasis * eVecs';

    %face list for bounding box
    faces = [ 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8; 1 2 3 4; 5 6 7 8 ];
    %display transparent bounding box
    patch('Vertices', box, 'Faces', faces, 'FaceVertexCData', (1:6)', 'FaceColor', 'flat', 'FaceAlpha', 0.08);

elseif dims == 2 
    plot(data_centered(:, 1), data_centered(:, 2), '.')
    labels = {'E1', 'E2'};
    
    %display arrows for eigenvectors
    quiver([0 0], [0 0], rangeVecs(1, :), rangeVecs(2, :), 'LineWidth', 2.5, 'Color', [0.27 0 0.58], 'AutoScale', 'off');
    
    %add text labels to eigenvectors
    text(textPos(1, :), textPos(2, :), labels, 'FontSize', 13)
    
    %get bounding box coordinates in PCA basis
    rect_PCAbasis = zeros(4, 2);
    rect_PCAbasis(1, :) = [ minima(1), minima(2) ];
    rect_PCAbasis(2, :) = [ maxima(1), minima(2) ];
    rect_PCAbasis(3, :) = [ maxima(1), maxima(2) ];
    rect_PCAbasis(4, :) = [ minima(1), maxima(2) ];
    
    %transform back to original basis
    rect = rect_PCAbasis * eVecs';
    
    fill(rect(:,1), rect(:,2), [0.27 0 0.58], 'FaceAlpha', 0.08);
    
end
axis equal
title('Covariance');



%% PCA using SVD
[u, eVals, eVecs] = svd(data_centered, 0);
eVals = diag(eVals);
[eVals, order] = sort(eVals, 'descend');
eVecs = eVecs(:, order);
subplot(1, 2, 2);
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
%size(eVecs)
%biplot(eVecs, 'LineWidth', 1.5, 'Color', [0.27 0 0.58], 'Marker', 'o', 'MarkerSize', 3); %'VarLabels', labels, 

% %matlab princomp
% subplot(1,2,2);
% hold on;
% if dims == 3
%     pcshow(data_centered);
%     labels = {'X1' 'X2' 'X3'};
% else
%     plot(data_centered(:, 1), data_centered(:, 2), 'g.')
%     labels = {'X1' 'X2'};
% end
% [pc,score,latent,tsquare] = princomp(data);
% biplot(pc(1:dims,1:dims), 'VarLabels', labels)
