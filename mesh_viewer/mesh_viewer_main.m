% @file     mesh_viewer_main.m
% @author   anna fruehstueck
% @date     10/02/2017

close all;
clc;
clear;
%[V, F] = read_obj('../data/shuttle.obj');          %  ~300 vertices
%[V, F] = read_obj('../data/magnolia.obj');         %  ~800 vertices
%[V, F] = read_obj('../data/teddy.obj');            % ~1500 vertices
%[V, F] = read_obj('../data/bunny.obj');            % ~2500 vertices
%[V, F] = read_obj('../data/teapot.obj');           % ~3600 vertices
%[V, F] = read_obj('../data/lamp.obj');             % ~4400 vertices
%[V, F] = read_obj('../data/cow.obj');              % ~4500 vertices
[V, F] = read_obj('../data/pumpkin.obj');          % ~5000 vertices
%[V, F] = read_obj('../data/suzanne.obj');          % ~7800 vertices
%[V, F] = read_obj('../data/trumpet.obj');          %~11900 vertices
%[V, F] = read_obj('../data/minicooper.obj');       %big!  ~44000 vertices
%[V, F] = read_obj('../data/tyrannosaurus.obj');	%BIG! ~100000 vertices
%[V, F] = read_obj('../data/armadillo.obj');        %BIG! ~106200 vertices

[VxV, FxV] = calc_adjacent_vertices(V, F);
FxF = calc_adjacent_faces(F);


% figure('Name', 'Vertex Adjacency')
% spy(S_v)
% figure('Name', 'Face Adjacency')
% spy(S_f)

N = calc_normal(V, F, FxV);
L = calc_laplacian('cotan', V, F, FxV, VxV);

figure;
trisurf(F, V(:, 1), V(:, 2), V(:, 3));

scr = get(0, 'ScreenSize'); 
fig = figure('Name', '2D Viewer', 'NumberTitle', 'off', 'Position', [scr(3)/2 50 scr(3)/3 scr(3)/3]);
hold on;
%axis vis3d
V = V + L
trisurf(F, V(:, 1), V(:, 2), V(:, 3));%,'FaceColor',[0.26,0.33,1.0 ]);
%quiver3(V(:, 1), V(:, 2), V(:, 3), N(:, 1), N(:, 2), N(:, 3))
%quiver3(V(:, 1), V(:, 2), V(:, 3), L(:, 1), L(:, 2), L(:, 3));

%add light sources
%light('Position', [0, 60.0, 30.0], 'Style', 'infinite');
%light('Position', [0, -60.0, -30.0], 'Style', 'infinite');

lighting phong;