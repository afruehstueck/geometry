% @file     auto2D_main.m
% @author   afruehstueck
% @date     02/02/2017
%
% asks for user input of 2D continuous line drawing and applies several
% transformations to the input data
% automatically generate a number of 2D transformation matrices
% and show animations of the applied transformations

clear;
clc;

[view, points] = viewer2D();
title('Draw an arbitrary object...');

%rotate
disp('***** ROTATION *****');
m = matrix( 2, 'rotation', 3*pi/512 );
analyzeMatrix(m);

title('Rotate 270° around origin');
for i = 1:256
    points = transformPoints(view, points, m);
end

%invert rotation
minv = inv(m)
title('Rotate back');
for i = 1:256
    points = transformPoints(view, points, minv);
end

%scale
disp('***** SCALE *****');
m = matrix( 2, 'scale', 1.002, 1.001 );
analyzeMatrix(m);

title('Anisotropic scaling');
for i = 1:256
    points = transformPoints(view, points, m);
end

%invert scaling
title('Inverse scaling');
minv = inv(m)
for i = 1:256
    points = transformPoints(view, points, minv);
end

%shearx
disp('***** SHEAR *****');
m = matrix( 2, 'shear', 'x', 0.005 )
analyzeMatrix(m);

title('Shear along x axis');
for i = 1:128
    points = transformPoints(view, points, m);
end

%invert shear
minv = inv(m)
title('Invert shear');
for i = 1:128
    points = transformPoints(view, points, minv);
end

%sheary
m = matrix( 2, 'shear', 'y', 0.006 )

title('Shear along y axis');
for i = 1:128
    points = transformPoints(view, points, m);
end

%invert shear
minv = inv(m)
title('Invert shear');
for i = 1:128
    points = transformPoints(view, points, minv);
end

cla(view);

title('Now draw something off-center...');
fH = imfreehand('closed', 0);
points = getPosition(fH);
points = points';

center = sum(points,2)/length(points);

t = matrix(3,'translation', -center(1)/200, -center(2)/200)
tinv = inv(t);
title('Translate to origin');
for i = 1:200
    points = transformPoints(view, points, t);
end

title('Rotate 90° around origin');
r = matrix(2, 'rotation', (pi/2)/200);
for i = 1:200
    points = transformPoints(view, points, r);
end

s = matrix(2, 'scale', 1.003, 1.003);
title('Scale');
for i = 1:128
    points = transformPoints(view, points, s);
end

title('Translate back');
for i = 1:200
    points = transformPoints(view, points, tinv);
end