% @file     ICP_main.m
% @author   afruehstueck
% @date     17/04/2017

close all
clc;
clear;

%% settings
%%%% BUNNY %%%%
path = '../data/pointcloud/bunny.ply';
plane_A = [ 0 0 -1; 0.3 0.1 0; -1 0 0.5 ];
plane_B = [-0.008 0.075 -1.005; 0.292 0.175 -0.005; -1.008 0.075  0.495];
threshold = 0.005;

%path = '../data/pointcloud/lobster.ply';
%path = '../data/pointcloud/pepper.ply';

%%%% TEAPOT %%%%
path = 'teapot.ply';
plane_A = [ 0.5 0.6 -1; 0.3 0.1 0; -1 0 0.5 ];
plane_B = [0.56 1.08 -0.748; 0.36 0.58 0.252; -0.94 0.48 0.752];
threshold = 0.1;


% %%%% FACE %%%%
% path = '../data/pointcloud/face.ply';
% plane_A = [ 0 0 -1; 0.3 0.1 0; -1 0 0.5 ];
% plane_B = [-0.008 0.075 -1.005; 0.292 0.175 -0.005; -1.008 0.075  0.495];
% threshold = 0.1;

pcloud = pcread(path);
data = pcloud.Location;
[M, ~] = size(data);
max_num = 10000;
if M > max_num
    sample = randsample(M, max_num);
    M = max_num;
    data = data(sample, :);
end

minima = min(data, [], 1);
maxima = max(data, [], 1);
ranges = maxima - minima;
   
cutPlane = true;            %apply cut in data interesting results
randomRotation = true;      %apply random rotation to data for more interesting results
randomTranslation = true;   %apply random translation to data for more interesting results
addNoise = true;            %apply noise to data poinfor more interesting results

max_noise = max(ranges) / 300;
max_rotation = pi/25;
max_translation = ranges / 25;

%%
data_A = data; 
data_B = data; 
if cutPlane
    side = -1;
    f = 0.4; %distance between planes
    
    plane_A_V1 = plane_A(2, :) - plane_A(1, :);
    plane_A_V2 = plane_A(3, :) - plane_A(1, :);
    N = cross(plane_A_V1, plane_A_V2);
    
%     plane_B(1, :) = plane_A(1, :) + side * f * N
%     plane_B(2, :) = plane_A(2, :) + side * f * N
%     plane_B(3, :) = plane_A(3, :) + side * f * N
    
    plane_B_V1 = plane_B(2, :) - plane_B(1, :);
    plane_B_V2 = plane_B(3, :) - plane_B(1, :);
    
    plane_A_V3s = data_A - plane_A(1, :);
    plane_B_V3s = data_B - plane_B(1, :);
    fprintf('size before cut plane: %d\n', M);
    new_data_A = [];
    new_data_B = [];
    for p = 1:M
        plane_A_V3 = plane_A_V3s(p, :);
        plane_B_V3 = plane_B_V3s(p, :);
        
        M_plane_A = [ plane_A_V1; plane_A_V2; plane_A_V3 ];
        D_plane_A = det(M_plane_A);
        
        M_plane_B = [ plane_B_V1; plane_B_V2; plane_B_V3 ];
        D_plane_B = det(M_plane_B);
        
        if sign(D_plane_A) == side
            new_data_A = [new_data_A; data_A(p, :)];
        end
        
        if sign(D_plane_B) ~= side
            new_data_B = [new_data_B; data_B(p, :)];
        end
    end
    [M, ~] = size(new_data_A);
    fprintf('size after cut plane 1: %d\n', M);
    fprintf('size after cut plane 2: %d\n', size(new_data_B, 1));
    data_A = new_data_A;
    data_B = new_data_B;
end

if addNoise
   data_A = data_A + rand(size(data_A)) * 2 * max_noise - max_noise; 
   data_B = data_B + rand(size(data_B)) * 2 * max_noise - max_noise;
end

if randomRotation 
    rotX = matrix(3, 'rotation', 'x', rand * 2*max_rotation - max_rotation);
    rotY = matrix(3, 'rotation', 'y', rand * 2*max_rotation - max_rotation);
    rotZ = matrix(3, 'rotation', 'z', rand * 2*max_rotation - max_rotation);
    r = rotX * rotY * rotZ;
    data_A = transformPoints(data_A', r)';
end

if randomTranslation
    randDX = rand * 2*max_translation(1) - max_translation(1);
    randDY = rand * 2*max_translation(2) - max_translation(2);
    randDZ = rand * 2*max_translation(3) - max_translation(3);
    fprintf('Translate by random vector (%f %f %f)\n', randDX, randDY, randDZ);
    data_A = transformPoints(data_A', matrix(4, 'translation', randDX, randDY, randDZ))';
end


%figure stuff
color1 = [0.27 0.00 0.58];
color2 = [1.00 0.80 0.00];
%color2 = [125 207 182]/255;

scr = get(0, 'ScreenSize'); 
figsize = [80 120 5*scr(4)/3 2*scr(4)/3];
figure('Name', 'ICP', 'NumberTitle', 'off', 'Position', figsize);  

plot1 = subplot(1,2,1);

hold on;
pcshow(data_A, color1);
pcshow(data_B, color2);
axis equal;
title('Data');

for iteration = 1:100    
    
    %%%% KNN SEARCH %%%
    modelEX = ExhaustiveSearcher(data_B);
    %kd tree is about 4x faster
    modelKD = KDTreeSearcher(data_B);

    %finds the nearest neighbor in data_B for each point in data_A -> data_B(idx(i),:) corresponds to data_A(i,:)
    %[idx, d] = knnsearch(modelEX, data_A); 
    [idx, d] = knnsearch(modelKD, data_A); 
    
    %find all data points in both data sets that are closer than a specific threshold
    similar_A = data_A(d < threshold, :);
    similar_B = data_B(idx(d < threshold), :);
    
    fprintf('%d of %d below threshold\n', size(similar_A, 1), size(data_A, 1));
    
    %subtract the mean from all data points
    centroid_A = mean(similar_A);
    centroid_B = mean(similar_B);
    
    data_centered_A = similar_A - centroid_A;
    data_centered_B = similar_B - centroid_B;
    
    C = data_centered_A' * data_centered_B;
    
    [U, S, V] = svd(C);
    
    R_opt = V * U';
    
    if det(R_opt) < 0
        %diag([ones(M-1, 1); -1])
        %R_opt = U * diag([1; 1; -1]) * V';
        disp('inverse rotation')
        R_opt(:, 3) = R_opt(:, 3) * -1;
    end
    
    t_opt = centroid_B' - R_opt * centroid_A';

    rotated = R_opt * data_A';
    translated = rotated + repmat(t_opt, 1, M);
    data_A = translated';
 
    %data_A = mean(data_B)' + R * data_centered_A'

    plot2 = subplot(1, 2, 2);
    cla;
    hold on;
    pcshow(data_A, color1);
    pcshow(data_B, color2);
    hold off;
    
    link = linkprop([plot1, plot2], {'CameraPosition', 'CameraUpVector'} );
    rotate3d on

    drawnow;
    
    error = similar_A - similar_B;
    error = error .* error;
    meansquarederr = sqrt(sum(error(:))/ size(similar_A, 1));
    fprintf('RMSE: %f\n', meansquarederr);
    if meansquarederr < 1e-3 
        fprintf('converged after %d iterations\n', iteration);
        return;
    end
end