% @file     regression_inputClick.m
% @author   afruehstueck
% @date     02/02/2017

%create a 2D view and plot a function
function [viewer, points] = regression_inputClick()
    clear;
    clc;
    
    dimX = [-5, 5];
    dimY = [-5, 5];

    scr = get(groot, 'ScreenSize');  
    fig = figure('Name', '2D Viewer', 'NumberTitle', 'off', 'Position', [scr(3)/2 50 scr(3)/3 scr(3)/3]);
    hold on;
    axis([dimX, dimY]);    
  
    x = zeros(1,1);
    y = zeros(1,1);
    r = linspace(-5, 5);
    num = 1;
    while ishandle(fig)
        [ x(num), y(num) ] = ginput(1);

%         prompt = 'Input degree of output polynomial: ';
%         str = input(prompt, 's');
%         if isempty(str)
%             deg = num-1
%         else
%             deg = str2num(str)
%         end
    
        A = ones(num);
         
        for col=1:num-1 
            A(:, col) = x'.^(num - col);
        end
        
        b = y';
        %X = linsolve(A, b);
        X = A \ b;
        
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