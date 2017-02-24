% @file     mesh_viewer_main.m
% @author   anna fruehstueck
% @date     10/02/2017

close all;
clc;
clear;
%[V, F] = read_obj('../data/shuttle.obj');          %  ~300 vertices
%[V, F] = read_obj('../data/magnolia.obj');         %  ~800 vertices
[V, F] = read_obj('../data/teddy.obj');            % ~1500 vertices
%[V, F] = read_obj('../data/bunny.obj');            % ~2500 vertices
%[V, F] = read_obj('../data/teapot.obj');           % ~3600 vertices
%[V, F] = read_obj('../data/lamp.obj');             % ~4400 vertices
%[V, F] = read_obj('../data/cow.obj');              % ~4500 vertices
%[V, F] = read_obj('../data/pumpkin.obj');          % ~5000 vertices
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

% N = calc_normal(V, F, FxV);
L = calc_laplacian('cotan', V, F, FxV, VxV);
scr = get(0, 'ScreenSize'); 
plt = zeros(1, 2);
original = figure('Name', 'Laplacian', 'NumberTitle', 'off', 'Position', [scr(3)/4 50 scr(3)/2 scr(3)/4]);
plt(1) = subplot(1,2,1);
hold on;
%trisurf(F, V(:, 1), V(:, 2), V(:, 3));
quiver3(V(:, 1), V(:, 2), V(:, 3), L(:, 1), L(:, 2), L(:, 3));
view(2)
plt(2) = subplot(1,2,2);
axis(plt, 'vis3d');
grid off;
for i = 1:10
    cla(gca);
    V = V + 0.2*L;
    hold on;
    trisurf(F, V(:, 1), V(:, 2), V(:, 3));%,'FaceColor',[0.26,0.33,1.0 ]);
    hold off;
    view(2)
    L = calc_laplacian('cotan', V, F, FxV, VxV);
    drawnow;
    %pause(0.1)
end


%quiver3(V(:, 1), V(:, 2), V(:, 3), N(:, 1), N(:, 2), N(:, 3))
%quiver3(V(:, 1), V(:, 2), V(:, 3), L(:, 1), L(:, 2), L(:, 3));

%add light sources
%light('Position', [0, 60.0, 30.0], 'Style', 'infinite');
%light('Position', [0, -60.0, -30.0], 'Style', 'infinite');
%lighting phong;