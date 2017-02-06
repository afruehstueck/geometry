% @file     regression.m
% @author   afruehstueck
% @date     02/02/2017

%create a 2D view and request a 2D usergenerated object to plot
function [viewer, points] = regression(dimX, dimY)
    clear;
    clc;
    
    if (~exist('dimX', 'var'))
        dimX = [-10, 10];
    end

    if (~exist( 'dimY', 'var'))
        dimY = [-10, 10];
    end

    scr = get(groot, 'ScreenSize');  
    fig = figure('Name', '2D Viewer', 'NumberTitle', 'off', 'Position', [scr(3)/2 50 scr(3)/3 scr(3)/3]);
    hold on;
    grid on;
    axis([dimX dimY]);

    x = zeros(1,1);
    y = zeros(1,1);
    
    num = 1;
    r=linspace(-10,10);
    while ishandle(fig)
        [ x(num), y(num) ] = ginput(1);

        A = ones(num);
         
        for col=1:num-1 
            A(:, col) = x'.^(num - col);
        end
        
        b = y';
        X = linsolve(A, b);
        
        cla(fig);
        
%         polynomial = '';
%         for row=num:-1:1
%             coeff = num2str(X(row), 3);
%             %exp = num2str(row-1);
%             exp = '';
%             if row > 2 
%                 exp = [ '*x^' num2str(row-1) ];
%             elseif row == 2
%                 exp = 'x';
%             end
%             polynomial = [polynomial coeff exp];
%             if row>1
%                 polynomial = [polynomial '+'];
%             end
%         end
%         title(polynomial);
        plot(x, y, '*');
        plot(r, polyval(X, r));
        num = num + 1;
    end
end