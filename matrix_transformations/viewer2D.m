% @file     viewer2D.m
% @author   afruehstueck
% @date     30/01/2017

%create a 2D view and request a 2D usergenerated object to plot
function [viewer, points] = viewer2D(dimX, dimY, closedLine)
    clear;
    clc;

    if (~exist('dimX', 'var'))
        dimX = [-5, 5];
    end

    if (~exist( 'dimY', 'var'))
        dimY = [-5, 5];
    end

    if (~exist('closedLine', 'var'))
        closedLine = 0;
    end

    scr = get(groot, 'ScreenSize');  
    viewer = figure('Name', '2D Viewer', 'NumberTitle', 'off', 'Position', [scr(3)/2 50 scr(3)/3 scr(3)/3]);
    title('Draw an arbitrary object...');
    hold on;
    grid on;
    axis([dimX dimY]);

    %input points manually
    %[ points( :, 1 ), points( :, 2 ) ] = getpts( view );
    %[ x, y ] = getpts( view );
    %plot( x, y );

    fH = imfreehand('closed', closedLine);
    points = getPosition(fH);
    points = points';
    
    cla(viewer);
    title('');
    
    plot(points(1, :), points(2, :), 'Color', [0.956, 0.258, 0.258]);
    drawnow;
end