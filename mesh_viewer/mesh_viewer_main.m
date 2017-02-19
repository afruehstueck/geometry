clc;
clear;
[V, F] = read_obj('../data/shuttle.obj');          %  ~300 vertices
%[V, F] = read_obj('../data/magnolia.obj');         %  ~800 vertices
%[V, F] = read_obj('../data/teddy.obj');            % ~1500 vertices
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

scr = get(0, 'ScreenSize'); 

fig = figure('Name', '2D Viewer', 'NumberTitle', 'off', 'Position', [scr(3)/2 50 scr(3)/3 scr(3)/3]);
hold on;
axis vis3d;
trisurf(F, V(:, 1), V(:, 2), V(:, 3));%,'FaceColor',[0.26,0.33,1.0 ]);

%add light sources
light('Position', [0, 60.0, 30.0], 'Style', 'infinite');
light('Position', [0, -60.0, -30.0], 'Style', 'infinite');

lighting phong;