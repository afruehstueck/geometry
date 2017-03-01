% @file     mesh_viewer_main.m
% @author   anna fruehstueck
% @date     10/02/2017

close all;
clc;
clear;
%[V, F] = read_obj('../data/mesh/shuttle.obj');          %  ~300 vertices
[V, F] = read_obj('../data/mesh/iris.obj');              % ~1100 vertices
%[V, F] = read_obj('../data/mesh/teddy.obj');            % ~1500 vertices
%[V, F] = read_obj('../data/mesh/moredata/doghead.obj');   
%[V, F] = read_obj('../data/mesh/head.obj');             % ~5100 vertices
%[V, F] = read_obj('../data/mesh/sphere.obj');           % ~2500 vertices
%[V, F] = read_obj('../data/mesh/sphere_distorted.obj'); % ~2500 vertices
%[V, F] = read_obj('../data/mesh/bunny.obj');            % ~2500 vertices
%[V, F] = read_obj('../data/mesh/teapot.obj');           % ~3600 vertices
%[V, F] = read_obj('../data/mesh/lamp.obj');             % ~4400 vertices
%[V, F] = read_obj('../data/mesh/cow.obj');              % ~4500 vertices
%[V, F] = read_obj('../data/mesh/pumpkin.obj');          % ~5000 vertices
%[V, F] = read_obj('../data/mesh/suzanne.obj');          % ~7800 vertices
%[V, F] = read_obj('../data/mesh/trumpet.obj');          %~11900 vertices
%[V, F] = read_obj('../data/mesh/minicooper.obj');       %big!  ~44000 vertices
%[V, F] = read_obj('../data/mesh/tyrannosaurus.obj');	%BIG! ~100000 vertices
%[V, F] = read_obj('../data/mesh/armadillo.obj');        %BIG! ~106200 vertices

[VxV, FxV] = calc_adjacent_vertices(V, F);
FxF = calc_adjacent_faces(F);

% figure('Name', 'Vertex Adjacency')
% spy(S_v)
% figure('Name', 'Face Adjacency')
% spy(S_f)

% N = calc_normal(V, F, FxV);

%find data range for axes
minel = min(V, [], 1);
maxel = max(V, [], 1);

type = 'cotan';

scr = get(0, 'ScreenSize'); 
plt = zeros(1, 2);
original = figure('Name', 'Laplacian', 'NumberTitle', 'off', 'Position', [scr(3)/4 50 scr(3)/2 scr(3)/4]);

plt(1) = subplot(1,3,1);

%colorbar

[Luni, Wuni, K] = calc_laplacian('uniform', V, F, FxV, VxV);
[Lcot, Wcot, K] = calc_laplacian('cotan', V, F, FxV, VxV);
[Vuni, D] = applyMatrices(V, Luni, Wuni);
[Vcot, D] = applyMatrices(V, Lcot, Wcot);

%D = L;
hold on;
trisurf(F, V(:, 1), V(:, 2), V(:, 3));
quiver3(V(:, 1), V(:, 2), V(:, 3), D(:, 1), D(:, 2), D(:, 3), 2);

% light('Position', [0, 60.0, 30.0], 'Style', 'infinite');
% light('Position', [0, -60.0, -30.0], 'Style', 'infinite');
% lighting phong;

%V = V + D;

%set fixed axes for both plots
axis([minel(1) maxel(1) minel(2) maxel(2)]);
plt(2) = subplot(1,3,2);
axis([minel(1) maxel(1) minel(2) maxel(2)]);
plt(3) = subplot(1,3,3);
axis([minel(1) maxel(1) minel(2) maxel(2)]);

hlink = linkprop(plt, {'CameraPosition','CameraUpVector'} );
rotate3d on

for i = 1:25%80
    subplot(1,3,2)
    cla(gca);
    hold on; 
    trisurf(F, Vuni(:, 1), Vuni(:, 2), Vuni(:, 3));
%     light('Position', [0, 60.0, 30.0], 'Style', 'infinite');
%     light('Position', [0, -60.0, -30.0], 'Style', 'infinite');
%     lighting phong;
    hold off;
    subplot(1,3,3)
    cla(gca);
    hold on;
    trisurf(F, Vcot(:, 1), Vcot(:, 2), Vcot(:, 3));
%     light('Position', [0, 60.0, 30.0], 'Style', 'infinite');
%     light('Position', [0, -60.0, -30.0], 'Style', 'infinite');
%     lighting phong;
    hold off;
    
%     [Luni, Wuni, K] = calc_laplacian('uniform', Vuni, F, FxV, VxV);
%     [Vuni, D] = applyMatrices(Vuni, Luni, Wuni);
    
    [Lcot, Wcot, K] = calc_laplacian('cotan', Vcot, F, FxV, VxV);
    [Vcot, D] = applyMatrices(Vcot, Lcot, Wcot);
    %D = L;
    %V = V + D;
    
    drawnow;
end


function [V, D] = applyMatrices(V, L, W)
    M = W * L;
    D = M * V;
    V = V + D;
end

%quiver3(V(:, 1), V(:, 2), V(:, 3), N(:, 1), N(:, 2), N(:, 3))
%quiver3(V(:, 1), V(:, 2), V(:, 3), L(:, 1), L(:, 2), L(:, 3));

%add light sources
