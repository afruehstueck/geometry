% @file     regression_main.m
% @author   afruehstueck
% @date     11/02/2017
%
% execute several instances of regression with different input functions
% note: there is user input requested in most of these functions
% abort input request by pressing 'Enter'

clear;
clc;
% input data points by mouse click in figure window
regression_inputClick(1, 0);

% specify degree of polynomial by input in console window (requests more
% input until input is aborted)


%leave some values out one x-axis
xl = linspace(-4, -1, 25);
xr = linspace(1.5, 3.5, 25);
x = [xl xr];
[y, y_noise] = noisyFunction(x, 0.8, @sin);
descr = 'f(x) = sin(x)';
regression_inputFn(x, y_noise, descr);

x = linspace(-2, 2, 30);
descr = 'polynomial f(x) = 3x^4 + x^3 - x^2 + x + 3';
[y, y_noise] = noisyFunction(x, 20., @polynomial, [3 1 -1 1 3]);
regression_inputFn(x, y_noise, descr);

%remove return for some more examples
return;

x = linspace(-2, 2, 50);
[y, y_noise] = noisyFunction(x, 0.4, @cos);
descr = 'f(x) = cos(x)';
regression_inputFn(x, y_noise, descr); 

x = linspace(-2, 2, 60);
descr = 'polynomial f(x) = -4x^5 + 2x^4 + 3x^3 - 3x^2 + 2x - 1';
[y, y_noise] = noisyFunction(x, 50., @polynomial, [-1 2 -3 3 2 -4]);
regression_inputFn(x, y_noise, descr);