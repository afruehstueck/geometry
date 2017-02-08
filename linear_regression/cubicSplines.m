% @file     cubicSplines.m
% @author   afruehstueck
% @date     07/02/2017

%create a 2D view and plot a function
function [viewer, points] = cubicSplines()
    clear;
    clc;

    dimX = [-10, 10];
    dimY = [-5, 5];

    scr = get(0, 'ScreenSize');  
    fig = figure('Name', '2D Viewer', 'NumberTitle', 'off', 'Position', [scr(3)/2 50 scr(3)/3 scr(3)/3]);
    hold on;
    axis([dimX dimY]);    
    
    x = [];
    y = [];
    
    %READ DATA FROM FILE
    data = dlmread('data_bunny.txt'); % read from ASCII file
    %data = dlmread('data_five.txt'); % read from ASCII file
    
    x = data(1:size(data,1), 1);
    y = data(1:size(data,1), 2);
    
    sizez = size(x, 1);
    lin = linspace(0, 4*pi, sizez);
    %choose z-coords along sine curve
    z = -100 + 200 * sin(lin);
    
    
    %SPIRAL
%     b = 0.1;
%     a = 0;
%     N = 4;
%     t = linspace(0, N*2*pi, 100)
%     cost = b *cos(t(:) + a)
%     sint = b *sin(t(:) + a)
%     x = t.* cost
%     y = t.* sint
    %choose z-coords randomly
    %z = -25 + 50 * rand(sizez, 1);
    
    %test points
    %x = [-2, -1.4, -1.0, -0.8, 0.8, 1.3, 1.5];
    %y = [-1.1, -3.4, -2., -0.3, 1.3, 1.8, 3.5];
   
    %x = -10:10;
    %y = cos(x).*abs(x);
    %[~, y] = noisyFunction(x, .7, @sin);
    %[~, y] = noisyFunction(x, 1000., @polynomial, [2.8 -0.3 1.4 3.2]);
    
    if length(x) == 0 && length(y) == 0 %do mousepicking
        mousePicking(fig);
    else %take predefined points
        axis auto;
        plot3(x, y, z, 'o');
        calculateSpline(x, y, z, 1);
        %calculateSpline(x, y, 2);
    end
end

function [] = mousePicking(fig) 
    % Give instructions
    disp('Left mouse button picks points.')
    disp('Right mouse button picks last point.')
    % Stores integer value for which button is pressed
    % 1 for left, 2 for middle, 3 for right
    
    button = 1; 
    pts = 0;
    
    while button == 1
        pts = pts+1;
        [x(pts), y(pts), button] = ginput(1); %Gets mouse click input
        
        cla(fig);
        plot(x, y,'o');
        
        if pts == 1 %calculate spline for >=2 points 
            continue
        end
        
        calculateSpline(x, y);
    end
end

function [] = calculateSpline(x, y, z, closed)
    if (~exist('closed', 'var'))
        closed = 0;
    end
    
    abcd_x = evaluate1DSpline(x, closed);
    abcd_y = evaluate1DSpline(y, closed);
    abcd_z = evaluate1DSpline(z, closed);
    spl = size(abcd_x, 1) %number of splines
    
    u = linspace(0, 1);
    u3 = u.^3;
    u2 = u.^2;
    
    for i = 1:spl %step through segments
        c_x = abcd_x(i,:);
        c_y = abcd_y(i,:);
        c_z = abcd_z(i,:);
        
        p_x = c_x(4).*u3 + c_x(3).*u2 + c_x(2).*u + c_x(1);
        p_y = c_y(4).*u3 + c_y(3).*u2 + c_y(2).*u + c_y(1);
        p_z = c_z(4).*u3 + c_z(3).*u2 + c_z(2).*u + c_z(1);

        %plot spline segment
        %plot(p_x, p_y,'-');
        plot3(p_x, p_y, p_z,'-')
    end
end

function abcd = evaluate1DSpline( v, closed ) 
    pts = length(v);
    spl = pts - 1; %number of splines

    M = zeros(pts);
    b = zeros(pts, 1);

%           M       .   D     =         b   
%     [ 2 1        ] [ D(1) ]   [ 3(x(2) - x(1))   ]
%     | 1 4 1      | | D(2) |   | 3(x(3) - x(1))   |
%     |   1 4 1    | |  .   | = |       .          |
%     |     .....  | |  .   |   |       .          |
%     |      1 4 1 | |  .   |   | 3(x(n) - x(n-2)) |
%     [        1 2 ] [ D(n) ]   [ 3(x(n) - x(n-1)) ]
       
    for i = 2:pts-1
        M(i, i-1:i+1) = [1 4 1];
        b(i) = 3*(v(i+1) - v(i-1));
    end

    closed 
    %first row
    if closed == 1
        M(1, 1:2) = [4 1];
        M(1, pts) = 1;
        b(1) = 3*(v(2) - v(pts));
    else
        M(1, 1:2) = [2 1];
        b(1) = 3*(v(2) - v(1));
    end

    %last row
    if closed == 1
        M(pts, 1) = 1;
        M(pts, pts-1:pts) = [1 4];
        b(pts) = 3*(v(1) - v(pts-1));
    else
        M(pts, pts-1:pts) = [1 2];
        b(pts) = 3*(v(pts) - v(pts-1));
    end
    
    %solve system
    D = M \ b;
    
    %store coefficients
    abcd = zeros(spl, 4);
    for i = 1:spl %step through segments
        a = v(i);
        b = D(i);
        c = 3*(v(i+1) - v(i)) - 2*D(i) - D(i+1);
        d = 2*(v(i) - v(i+1)) + D(i) + D(i+1);
        abcd(i, :) = [a b c d];
    end   
    
    if closed == 1 %add another segment for closing the curve
        a = v(pts);
        b = D(pts);
        c = 3*(v(1) - v(pts)) - 2*D(pts) - D(1);
        d = 2*(v(pts) - v(1)) + D(pts) + D(1);
        
        abcd(spl+1, :) = [a b c d];
    end
end
