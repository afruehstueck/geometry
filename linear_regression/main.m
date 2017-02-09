% @file     main.m
% @author   afruehstueck
% @date     07/02/2017

function [] = main() 
    clear;
    clc;

    dimX = [-10, 10];
    dimY = [-5, 5];
    scr = get(0, 'ScreenSize');  
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   FIVE dataset (94 data points)             %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %READ DATA FROM FILE
    data = dlmread('data_five.txt'); % read points from ASCII file
    
    x = data(1:size(data,1), 1);
    y = data(1:size(data,1), 2);
    
    sizez = size(x, 1);
    lin = linspace(0, 6*pi, sizez);
    %choose z-coords along sine curve
    z = -100 + 200 * sin(lin);
    
    fig_five = figure('Name', 'Five Viewer', 'NumberTitle', 'off', 'Position', [50 150 scr(3)/2 scr(3)/3]);
    hold on;
    axis auto;
    plot3(x, y, z, 'o');
    cubicSplines(x, y, z, 1);    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   BUNNY dataset (73 data points)            %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %READ DATA FROM FILE
    data = dlmread('data_bunny.txt'); % read points from ASCII file
    x = data(1:size(data,1), 1);
    y = data(1:size(data,1), 2);
    
    sizez = size(x, 1);
    lin = linspace(0, 8*pi, sizez);
    %choose z-coords along sine curve
    z = -2*pi + 4*pi * sin(lin);
    
    fig_bunny = figure('Name', 'Bunny Viewer', 'NumberTitle', 'off', 'Position', [50 150 scr(3)/2 scr(3)/3]);
    hold on;
    axis auto;
    plot3(x, y, z, 'o');
    cubicSplines(x, y, z, 1);    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   SPIRAL (generated)                        %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    b = 0.3; %tightness
    a = 0; %angle
    N = 7; %num curls
    t = linspace(0, N*2*pi, 70);
    x = t'.* (b * cos(t(:) + a));
    y = t'.* (b * sin(t(:) + a));
    z = -t; %spiral upwards
    
    fig_spiral = figure('Name', 'Spiral Viewer', 'NumberTitle', 'off', 'Position', [50 150 scr(3)/2 scr(3)/3]);
    hold on;
    axis auto;
    plot3(x, y, z, 'o');
    cubicSplines(x, y, z, 0); 
    
    %test points
    %x = [-2, -1.4, -1.0, -0.8, 0.8, 1.3, 1.5];
    %y = [-1.1, -3.4, -2., -0.3, 1.3, 1.8, 3.5];
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   2D GENERATED NOISY FUNCTION               %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    x = linspace(-7*pi,7*pi,12*pi);
    abscos = @(x)cos(x).*abs(x);
    [foo, y] = noisyFunction(x, 6., abscos);
    %[foo, y] = noisyFunction(x, 500., @polynomial, [2.8 -0.3 1.4 3.2]);
    fig_poly = figure('Name', '2D Function', 'NumberTitle', 'off', 'Position', [50 150 scr(3)/2 scr(3)/3]);
    hold on;
    axis auto;
    plot(x, y, 'o');
    cubicSplines(x, y); 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   MOUSE PICKING                             %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    fig_picking = figure('Name', 'Pick points', 'NumberTitle', 'off', 'Position', [50 150 scr(3)/2 scr(3)/3]);
    hold on;
    axis([dimX dimY]);  
    
    mousePicking(fig_picking, @cubicSplines);
end

function [x, y] = mousePicking(fig, fun) 
    % Give instructions
    disp('Left mouse button picks points.')
    disp('Right mouse button picks last point.')
    % Stores integer value for which button is pressed
    % 1 for left, 2 for middle, 3 for right
    button = 1; 
    pts = 0;
    
    while button == 1
        pts = pts+1;
        try %catch figure deletion error
            [x(pts), y(pts), button] = ginput(1); %Gets mouse click input
            cla(fig);
            plot(x, y, 'o');

            if pts == 1 %calculate for >=2 points 
                continue
            end

            fun(x, y);
        catch
            button = 0;
            fprintf('Figure closed.');
        end
    end
end