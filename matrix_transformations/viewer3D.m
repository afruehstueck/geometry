% @file     viewer3D.m
% @author   afruehstueck
% @date     30/01/2017

%create a 3D view and show a 3D object in it
function [viewer, points] = viewer3D(dimX, dimY, dimZ)
    if (~exist('dimX', 'var'))
        dimX = [-4, 4];
    end

    if (~exist( 'dimY', 'var'))
        dimY = [-4, 4];
    end
    
    if (~exist( 'dimZ', 'var'))
        dimZ = [-4, 4];
    end

    scr = get(groot, 'ScreenSize');  
    viewer = figure('Name', '3D Viewer', 'NumberTitle', 'off', 'Position', [scr(3)/2 50 scr(3)/3 scr(3)/3]);
    hold on;
    grid on;
    
    axis([dimX dimY dimZ]);
    view(45, 15);

    %[x, y, z] = sphere()
    t = 0 : pi/15 : pi;
    [x, y, z] = cylinder(sin(t) - pi/3 - 0.2);
    
    points(:, :, 1) = x;
    points(:, :, 2) = y;
    points(:, :, 3) = 2*z - 1;
    surf(points(:, :, 1),points(:, :, 2),points(:, :, 3));
    
    drawnow;
end