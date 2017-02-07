% @file     regression_inputClick.m
% @author   afruehstueck
% @date     02/02/2017

%create a 2D view and plot a function
function [] = regression_inputClick()
    allowInput = 1;
    clear;
    clc;
    
    dimX = [-5, 5];
    dimY = [-5, 5];

    scr = get(0, 'ScreenSize');  
    fig = figure('Name', '2D Viewer', 'NumberTitle', 'off', 'Position', [scr(3)/2 50 scr(3)/3 scr(3)/3]);
    hold on;
    axis([dimX, dimY]);    
  
    x = zeros(1,1);
    y = zeros(1,1);
    r = linspace(-5, 5);
    num = 1;
    
    % Give instructions
    disp('Left mouse button picks points.')
    disp('Right mouse button picks last point.')
    % Stores integer value for which button is pressed
    % 1 for left, 2 for middle, 3 for right
    
    button = 1; 
    while button == 1
        [x(num), y(num), button] = ginput(1);
    
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
    
    if allowInput == 1
        while ishandle(fig)    
            prompt = 'Input degree of output polynomial: ';
            str = input(prompt, 's');
            if isempty(str), break, end

            deg = str2num(str);

            num = length(x);
            A = ones(num, deg + 1);

            for col=1:deg 
                A(:, col) = x'.^(deg + 1 - col);
            end

            b = y';
            %X = linsolve(A, b);
            X = A \ b;

            cla(fig);

            plot(x, y, '*');
            plot(r, polyval(X, r));
        end
    end
end