% @file     mesh_viewer_main.m
% @author   anna fruehstueck
% @date     10/02/2017

close all;
clc;
clear;

set(0,'DefaultFigureColormap', viridis)
[V, F] = read_obj('../data/mesh/simple_bunny.obj'); lambdauni = 0.1; lambdacot = 0.02;  % ~500 vertices
%[V, F] = read_obj('../data/mesh/sphere_distorted.obj'); % ~2500 vertices
%[V, F] = read_obj('../data/mesh/ankylosaurus.obj'); lambdauni = 0.5; lambdacot = 0.05;   % ~5100 vertices
%[V, F] = read_obj('../data/mesh/pumpkin.obj');          % ~5000 vertices


%[V, F] = read_obj('../data/mesh/shuttle.obj');          % ~300 vertices
%[V, F] = read_obj('../data/mesh/simple_bunny.obj');     % ~500 vertices
%%[V, F] = read_obj('../data/mesh/duck.obj');            % ~900 vertices
%[V, F] = read_obj('../data/mesh/LEGO_man.obj');    
%[V, F] = read_obj('../data/mesh/iris.obj');             % ~1100 vertices
%[V, F] = read_obj('../data/mesh/teddy.obj');            % ~1500 vertices
%[V, F] = read_obj('../data/mesh/gecko.obj');            % ~2200 vertices
%[V, F] = read_obj('../data/mesh/sphere.obj');           % ~2500 vertices
%[V, F] = read_obj('../data/mesh/sphere_distorted.obj'); % ~2500 vertices // holes
%!%[V, F] = read_obj('../data/mesh/camel.obj');            % ~3700 vertices
%%[V, F] = read_obj('../data/mesh/bear_pnd.obj');            % ~3900 vertices
%[V, F] = read_obj('../data/mesh/octopus.obj');          % ~4000 vertices
%[V, F] = read_obj('../data/mesh/lamp.obj');             % ~4400 vertices // holes
%[V, F] = read_obj('../data/mesh/cow.obj');              % ~4500 vertices
%[V, F] = read_obj('../data/mesh/atenea.obj');           % ~4700 vertices
%[V, F] = read_obj('../data/mesh/head.obj');             % ~5100 vertices
%[V, F] = read_obj('../data/mesh/sabtooth.obj');            % ~5800 vertices
%!%[V, F] = read_obj('../data/mesh/mammoth.obj');         % ~7400 vertices
%[V, F] = read_obj('../data/mesh/suzanne.obj');          % ~7800 vertices

max_vertices = max(V, [], 1);
min_vertices = min(V, [], 1);
range_vertices = max_vertices - min_vertices;

%normalize vertex values to [-5, 5] range for any plot
for dim = 1:3
    V(:, dim) = 10*(V(:, dim) - min_vertices(dim)) / range_vertices(dim) - 5;
end

[VxV, FxV] = calc_adjacent_vertices(V, F);

%inspect vertex and face adjacency matrices
% figure('Name', 'Vertex Adjacency')
% spy(S_v)
% figure('Name', 'Face Adjacency')
% spy(S_f)

scr = get(0, 'ScreenSize'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% normal vectors
N = calc_normal(V, F, FxV);

%show normal vectors on separate plot
normals = figure('Name', 'Normal vectors', 'NumberTitle', 'off', 'Position', [scr(3)/4 30 scr(3)/4 scr(3)/4]);
hold on;
trisurf(F, V(:, 1), V(:, 2), V(:, 3));
quiver3(V(:, 1), V(:, 2), V(:, 3), N(:, 1), N(:, 2), N(:, 3));
title('Vertex normals');
camlight
lighting gouraud;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% settings for laplacian

type1 = 'uniform';
type2 = 'cotan';
showoriginal = true;
update1 = true;
update2 = true;

numplots = showoriginal + update1 + update2;
plot1 = showoriginal + update1;
plot2 = numplots;

numIter = 25;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% update laplacian 
plt = zeros(1, 2);
laplacian_fig = figure('Name', 'Laplacian', 'NumberTitle', 'off', 'Position', scr);
%colorbar

%find data range for axes
minel = min(V, [], 1);
maxel = max(V, [], 1);

Vcot = V;
Vuni = V;

[Mcot, Dcot, Kcot] = calc_laplacian(type2, V, F, FxV, VxV);
%Vcot = applyCotanMatrices(Vcot, Mcot, Dcot, lambdacot);

[~, LapVcot, ~] = applyMatrices(V, Mcot, Dcot, lambdacot);
Klims = [-3, 3];

if showoriginal 
    plt(1) = subplot(1, numplots, 1);
    hold on;
    surface = trisurf(F, V(:, 1), V(:, 2), V(:, 3), Kcot);
    caxis(Klims)
    colorbar
    shading interp
    set(surface', 'edgecolor', 'k');
    quiver3(V(:, 1), V(:, 2), V(:, 3), LapVcot(:, 1), LapVcot(:, 2), LapVcot(:, 3), 2);
    camproj perspective
end

%set fixed axes for both plots
axis([minel(1) maxel(1) minel(2) maxel(2)]);
plt(2) = subplot(1,numplots,numplots-1);
title(type1)
axis([minel(1) maxel(1) minel(2) maxel(2)]);
plt(3) = subplot(1,numplots,numplots);
title(type2)   
axis([minel(1) maxel(1) minel(2) maxel(2)]);

%link rotation for all subplots
hlink = linkprop(plt, {'CameraPosition','CameraUpVector'} );
rotate3d on

for i = 1:numIter
    disp(['***** Iteration ', num2str(i), '/', num2str(numIter), ' *****']);
    if update1
        %redraw plot1
        subplot(1, numplots, plot1);
        cla(gca);
        hold on; 
        trisurf(F, Vuni(:, 1), Vuni(:, 2), Vuni(:, 3));  
        hold off;
        %recalculate laplacian1
        [Muni, Duni, Kuni] = calc_laplacian(type1, Vuni, F, FxV, VxV);
        [Vuni, LapVuni, Luni] = applyMatrices(Vuni, Muni, Duni, lambdauni);
        
    end
    
    if update2
        %redraw plot2
        subplot(1, numplots, plot2);
        cla(gca);
        
        hold on;
        surface = trisurf(F, Vcot(:, 1), Vcot(:, 2), Vcot(:, 3), Kcot);
        caxis(Klims);
        %caxis(Klims)
        shading interp     
        set(surface', 'edgecolor', 'k');
        hold off;

        %recalculate laplacian2
        [Mcot, Dcot, Kcot] = calc_laplacian(type2, Vcot, F, FxV, VxV);
        Vcot = applyCotanMatrices(Vcot, Mcot, Dcot, lambdacot);
        %Hcot = sqrt(sum(LapVcot.^2, 2)) / 2;
    end
    drawnow;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% eigen decomposition
[eVecs, eVals] = eig(full(Luni));
evecs_fig = figure('Name', 'Eigenvectors projected onto mesh', 'NumberTitle', 'off', 'Position', scr);

%sort eigenvectors according to eigenvalue sizes
eVals = diag(eVals); 
% sort the variances in decreasing order 
[~, order] = sort(eVals, 'descend');
eVecs = eVecs(:, order);

%show eigenvectors color coded on mesh
eig_plts = zeros(0, 2);
selected_eigenvectors = [1, 2, 3, 5, 10, length(eVecs) - 10, length(eVecs) - 5, length(eVecs )- 2, length(eVecs) - 1, length(eVecs)];
for i = 1:length(selected_eigenvectors)
    eig_plts(i) = subplot(2, length(selected_eigenvectors)/2, i);
    c = eVecs(:, selected_eigenvectors(i));
    surface = trisurf(F, V(:, 1), V(:, 2), V(:, 3), c);
    title(['Eigenvector #', num2str(selected_eigenvectors(i))])
    shading interp
    lighting gouraud
    camlight
    colormap jet
end

%link rotation for all subplots
hlink2 = linkprop(eig_plts, {'CameraPosition', 'CameraUpVector'} );
rotate3d on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% spectral mesh compression
figure('Name', 'Low-pass filter', 'NumberTitle', 'off', 'Position', scr);
coeff_plts = zeros(0, 2);
to_keep = [0.01, 0.02, 0.05, 0.15]; % percents of coefficients kept
for i = 1:length(to_keep)
    coeff_plts(i) = subplot(1, length(to_keep), i);
    
    pf = (eVecs' * V); % projection of the vector in the laplacian basis
    q = round(to_keep(i) * length(V)); % number of zeros
    %set all coefficients above percentage to zero
    pf(q+1:end, :) = 0;
    
    Vrec = (eVecs * pf);
    
    trisurf(F, Vrec(:, 1), Vrec(:, 2), Vrec(:, 3));
    title([num2str(to_keep(i) * 100), '% of coefficients'])
end

%link rotation for all subplots
hlink3 = linkprop(coeff_plts, {'CameraPosition', 'CameraUpVector'} );
rotate3d on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% enhanced eigenvectors
[eVecs, ~] = eig(full(Muni));

pf = (eVecs' * V); % projection of the vector in the laplacian basis
for i =length(pf)-2:length(pf)
    pf(i, :) = (1 + i * 0.005) * pf(i, :); % number of zeros
end
Vrec = (eVecs * pf);
comp_fig = figure('Name', 'Enhanced eigenvectors', 'NumberTitle', 'off');
trisurf(F, Vrec(:, 1), Vrec(:, 2), Vrec(:, 3));

%helper function to apply Laplacian and mass matrix to data
%returns new vertex positions and discrete average of laplace-beltrami
function [V, Dvec, L] = applyMatrices(V, M, D, lambda)
    L = D * M;
    Dvec = lambda * L * V;
    V = V + Dvec;
end

%helper function to apply Laplacian and mass matrix to data
%returns new vertex positions and discrete average of laplace-beltrami
function V = applyCotanMatrices(V, M, D, lambda)
    b = D * V;
    A = (D + lambda * D * M);
    V = A \ b;
end
