% @file     regression_inputFn.m
% @author   afruehstueck
% @date     06/02/2017

%create a 2D view and plot a function
function [viewer, points] = regression_inputFn()
    clear;
    clc;
    
    if (~exist('dimX', 'var'))
        dimX = [-5, 5];
    end

    if (~exist( 'dimY', 'var'))
        dimY = [-5, 5];
    end

    scr = get(groot, 'ScreenSize');  
    fig = figure('Name', '2D Viewer', 'NumberTitle', 'off', 'Position', [scr(3)/2 50 scr(3)/3 scr(3)/3]);
    hold on;
    axis auto;    
    
    x = linspace(-2, 2);
    
    %[y, y_noise] = noisyFunction(x, 0.4, @sin);
    %[y, y_noise] = noisyFunction(x, 0.4, @cos);
    %[y, y_noise] = noisyFunction(x, 2., @polynomial, [0.8 -0.3 0.4 0.5 -1.2]);
    [y, y_noise] = noisyFunction(x, 20., @polynomial, [2.8 -0.3 1.4 3.2]);
    
    num = length(x);  
    
    plot(x, y_noise, '*');
    drawnow;
    
    prompt = 'Input degree of output polynomial: (press enter to exit) ';
    while ishandle(fig)    
        str = input(prompt, 's');
        if isempty(str), break, end
        
        deg = str2num(str);
        
        A = ones(num, deg + 1);

        for col=1:deg 
            A(:, col) = x'.^(deg + 1 - col);
        end

        b = y_noise';
        X = A \ b;

        cla(fig);

        plot(x, y_noise, '*');
        plot(x, polyval(X, x));
        drawnow;
    end
end