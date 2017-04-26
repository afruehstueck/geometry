% @file     piecewiseSplines.m
% @author   afruehstueck
% @date     07/02/2017
%
% piecewise splines for points drawn in 2D
% not a parametric function

function [viewer, points] = piecewiseSplines()
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
        calculateSpline(x, y);
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
        
        if pts == 1 %calculate spline for >=2 points 
            continue
        end
        
        cla(fig);
        calculateSpline(x, y);
    end
end

function [] = calculateSpline( x, y, boundary )
    if (~exist('boundary', 'var'))
        boundary = 1;
    end
    pts = length(x);
    spl = pts - 1; %number of splines

    msz = 4*spl; %matrix size
    A = zeros(msz);
    b = zeros(msz, 1);

    for i = 1:spl %fill matrix 'cells' for each segment
        idx = (i-1)*4+1; %upper left matrix position of current segment
        n = i;      %current point index
        m = i+1;    %next point index (for readability)

        A(idx,   idx:idx+3) = [x(n)^3	x(n)^2	x(n) 1]; %p(x_i)   = y_i
        A(idx+1, idx:idx+3) = [x(m)^3   x(m)^2	x(m) 1]; %p(x_i+1) = y_i+1
        b(idx)      = y(i);
        b(idx+1)	= y(i+1);

        if i<spl %relations to next segment
            A(idx+2, idx:idx+7)	= [3*x(m)^2	2*x(m) 1 0	-3*x(m)^2	-2*x(m)	-1 0]; %first derivative has to be equal
            A(idx+3, idx:idx+7)	= [6*x(m)	2	   0 0	-6*x(m)     -2     	 0 0]; %second derivative has to be equal
        end
    end

    if boundary == 1    
        A(msz-1, 1:4)       = [3*x(1)^2     2*x(1)	 1 0]; %first derivative = 0 condition for first point
        A(msz,   msz-3:msz)	= [3*x(pts)^2	2*x(pts) 1 0]; %first derivative = 0 condition for last point
    elseif boundary == 2
        A(msz-1, 1:4)       = [6*x(1)	2 0 0]; %second derivative = 0 condition for first point
        A(msz,   msz-3:msz)	= [6*x(pts) 2 0 0];  %second derivative = 0 condition for last point
    end

    X = A \ b;
    
    plot(x, y, 'o'); %plot control points as red circles
    
    for i = 1:spl %step through segments
        idx = (i-1)*4+1;
        u = linspace(x(i),x(i+1));
        pu = X(idx).*u.^3 + X(idx+1).*u.^2 + X(idx+2).*u + X(idx+3);
        
        %plot spline segment
        plot(u,pu,'-');
    end
    
%     xx=linspace(x(1),x(pts));
%     pp = csape(x,y,'clamped',[0 0]);
%     yy0 = ppval(pp,xx);
%     plot(xx,yy0,'.r')
    %plot(x,y,'go')
    
    %MATLAB spline function (for comparison)
    %it has a different boundary condition
%     xx=linspace(x(1),x(pts));
%     zz=spline(x,y,xx);
%     plot(xx,zz,'g--');
    %xx=linspace
end        