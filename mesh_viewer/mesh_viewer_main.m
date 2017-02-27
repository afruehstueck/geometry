% @file     mesh_viewer_main.m
% @author   anna fruehstueck
% @date     10/02/2017

close all;
clc;
clear;
%[V, F] = read_obj('../data/shuttle.obj');          %  ~300 vertices
%[V, F] = read_obj('../data/magnolia.obj');         %  ~800 vertices
[V, F] = read_obj('../data/teddy.obj');            % ~1500 vertices
%[V, F] = read_obj('../data/moredata/doghead.obj');   
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

minel = min(V, [], 1)
maxel = max(V, [], 1)
type = 'cotan';
scr = get(0, 'ScreenSize'); 
plt = zeros(1, 2);
original = figure('Name', 'Laplacian', 'NumberTitle', 'off', 'Position', [scr(3)/4 50 scr(3)/2 scr(3)/4]);
plt(1) = subplot(1,2,1);
hold on;
trisurf(F, V(:, 1), V(:, 2), V(:, 3));
%colorbar

[L, W, K] = calc_laplacian(type, V, F, FxV, VxV);
% size(L)
% size(W)
% L(1,:)
% W(1,:)
[Vn, D] = applyMatrices(V, L, W);
%D = L;
quiver3(V(:, 1), V(:, 2), V(:, 3), D(:, 1), D(:, 2), D(:, 3), 2);
%V = V + D;

V = Vn;

axis([minel(1) maxel(1) minel(2) maxel(2)])
plt(2) = subplot(1,2,2);
axis([minel(1) maxel(1) minel(2) maxel(2)])

hlink = linkprop(plt,{'CameraPosition','CameraUpVector'});
rotate3d on

for i = 1:10%80
    cla(gca);
    hold on;    
    trisurf(F, V(:, 1), V(:, 2), V(:, 3));%,'FaceColor',[0.26,0.33,1.0 ]);
    hold off;
    
    [L, W, K] = calc_laplacian(type, V, F, FxV, VxV);
    [V, D] = applyMatrices(V, L, W);
    
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
%light('Position', [0, 60.0, 30.0], 'Style', 'infinite');
%light('Position', [0, -60.0, -30.0], 'Style', 'infinite');
%lighting phong;