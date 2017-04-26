% @file     splines_main.m
% @author   afruehstueck
% @date     07/02/2017
%
% generate several different examples of natural splines in 2D and 3D
% first viewer will request 2D user input and generate splines connecting
% these data points. Click right button to stop input or close figure
% window.
%
% Several other 2D and 3D examples are generated in various figures.

function [] = splines_main() 
    clear;
    clc;

    dimX = [-10, 10];
    dimY = [-5, 5];
    scr = get(0, 'ScreenSize');  
          
%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % LOAD AND CONNECT 3D VERTICES fun, but makes only limited sense ;) %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     [V, ~] = read_obj('../data/mesh/cow.obj');
%     x = V(:, 1);
%     y = V(:, 2);
%     z = V(:, 3);
%     fig_teddy = figure('Name', 'Teddy Viewer', 'NumberTitle', 'off', 'Position', [50 150 scr(3)/2 scr(3)/3]);
%     hold on;
%     axis auto;
%     %plot3(x, y, z, 'o');
%     naturalSplines(x, y, z, 1);   
    
%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   BUNNY dataset (73 data points)            %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    %READ DATA FROM FILE
    data = dlmread('../data/points/data_bunny.txt'); % read points from ASCII file
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
    naturalSplines(x, y, z, 1);    
    
    rotate3d on;
    
%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   FIVE dataset (94 data points)             %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %READ DATA FROM FILE
    data = dlmread('../data/points/data_five.txt'); % read points from ASCII file
    
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
    naturalSplines(x, y, z, 1);    
    
    rotate3d on;
    
    %%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   LOXODROME   (generated)                   %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    a = 0.13; %'curliness'
    t = linspace(-10*pi, 10*pi, 100);
    c = atan(0.13*t);
    x = cos(t).*cos(c);
    y = sin(t).*cos(c);
    z = -sin(c);
    
    fig_spirals = figure('Name', 'Loxodrome Viewer', 'NumberTitle', 'off', 'Position', [50 150 scr(3)/2 scr(3)/3]);
    hold on;
    axis auto;
    plot3(x, y, z, 'o');
    naturalSplines(x, y, z, 0); 
    
    rotate3d on;

%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   TWO SPIRALS   (generated)                 %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    b = 0.3; %tightness
    a = 0; %angle
    N = 5; %num curls
    t = linspace(-N*2*pi, N*2*pi, 100);
    t2 = [linspace(0, N*2*pi, 50) linspace(N*2*pi, 0, 50)];
    x = t'.* (b * sin(t2(:) + a));
    y = t'.* (b * cos(t2(:) + a));
    z = t; %spiral upwards
    
    fig_spirals = figure('Name', 'Two Spiral Viewer', 'NumberTitle', 'off', 'Position', [50 150 scr(3)/2 scr(3)/3]);
    hold on;
    axis auto;
    plot3(x, y, z, 'o');
    naturalSplines(x, y, z, 0); 
    
    rotate3d on;
    
%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   SPIRAL (generated)                        %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%    %commented because two spirals looks even cooler
%     b = 0.3; %tightness
%     a = 0; %angle
%     N = 7; %num curls
%     t = linspace(0, N*2*pi, 70);
%     x = t'.* (b * cos(t(:) + a));
%     y = t'.* (b * sin(t(:) + a));
%     z = -t; %spiral upwards
%     
%     fig_spiral = figure('Name', 'Spiral Viewer', 'NumberTitle', 'off', 'Position', [50 150 scr(3)/2 scr(3)/3]);
%     hold on;
%     axis auto;
%     plot3(x, y, z, 'o');
%     naturalSplines(x, y, z, 0); 
    
    %test points
    %x = [-2, -1.4, -1.0, -0.8, 0.8, 1.3, 1.5];
    %y = [-1.1, -3.4, -2., -0.3, 1.3, 1.8, 3.5];
   
%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    naturalSplines(x, y); 
    
%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   MOUSE PICKING                             %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    fig_picking = figure('Name', 'Pick points', 'NumberTitle', 'off', 'Position', [50 150 scr(3)/2 scr(3)/3]);
    hold on;
    axis([dimX dimY]);  
    
    mousePicking(fig_picking, @naturalSplines);
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