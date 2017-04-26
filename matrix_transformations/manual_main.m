% @file     manual_main.m
% @author   afruehstueck
% @date     30/01/2017
%
% investigate 2D/3D viewers and matrices manually
% start a 2D or 3D viewer (depending on user input), then request user to
% input properties for transformation matrices and apply these matrices to
% the input data
clear;
clc;

prompt = 'Display 2D or 3D? [2/3]: ';
str = input(prompt, 's');
if isempty(str)
    str = '2';
end

if strcmp(str, '2') || strcmp(str, '2D')
    [view, points] = viewer2D();
elseif strcmp(str, '3') || strcmp(str, '3D')
    [view, points] = viewer3D();
end

fprintf(2, '***************************************************************************\n');
fprintf(1, 'Please input a matrix in the following format:\n');
fprintf(2, 'm = matrix( dimensions, matrix type, additional parameters )\n\n');
fprintf(1, 'dimensions = { 2, 3, 4 (3D + homogeneous coordinates) }\n');
fprintf(1, 'matrix type = { ''rotation'', ''shear'', ''reflection'', ''scale'', ''translation'' }\n\n');
fprintf(1, 'apply the matrix using function ');
fprintf(2, 'points = transformPoints(view, points, m);\n');
fprintf(2, '***************************************************************************\n');

    