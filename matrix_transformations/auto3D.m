% @file     auto3D.m
% @author   afruehstueck
% @date     02/02/2017

%automatically generate a number of 3D transformation matrices
%and show animations of the applied transformations
clear;
clc;

[view, points] = viewer3D();

%rotate
disp('***** ROTATION *****');
m1 = matrix(3, 'rotation', 'x', 0.01);
m2 = matrix(3, 'rotation', 'y', 0.01);
m3 = matrix(3, 'rotation', 'z', 0.01);
analyzeMatrix(m1);

%invert rotation
minv1 = inv(m1)
minv2 = inv(m2)
minv3 = inv(m3)

title('Rotate around x axis');
for i = 1:128
    points = transformPoints(view, points, m1);
end

title('Rotate around y axis');
for i = 1:64
    points = transformPoints(view, points, m2);
end

title('Rotate around z axis');
for i = 1:128
    points = transformPoints(view, points, m3);
end

title('Rotate backwards around z axis');
for i = 1:128
    points = transformPoints(view, points, minv3);
end

title('Rotate backwards around y axis');
for i = 1:64
    points = transformPoints(view, points, minv2);
end

title('Rotate backwards around x axis');
for i = 1:128
    points = transformPoints(view, points, minv1);
end

title('Rotate around all axes');
for i = 1:64
    points = transformPoints(view, points, m1);
    points = transformPoints(view, points, m2);
    points = transformPoints(view, points, m3);
end

%scale
disp('***** SCALE *****');
m = matrix( 3, 'scale', 1.006, 1.004, 1.003 );
analyzeMatrix(m);
minv = inv(m)

title('Anisotropic scaling along all axes');
for i = 1:128
    points = transformPoints(view, points, m);
end

%invert scaling
title('Scale back');
for i = 1:128
    points = transformPoints(view, points, minv);
end

%shear
disp('***** SHEAR *****');
m1 = matrix( 3, 'shear', 'x', 0.008, 0.003 )
m2 = matrix( 3, 'shear', 'y', 0.004, -0.005 )
m3 = matrix( 3, 'shear', 'z', -0.002, 0.007 )

analyzeMatrix(m1);

minv1 = inv(m1)
minv2 = inv(m2)
minv3 = inv(m3)

title('Shear along x axis');
for i = 1:100
    points = transformPoints(view, points, m1);
end

title('Shear along y axis');
for i = 1:100
    points = transformPoints(view, points, m2);
end

% title('Shear along z axis');
% for i = 1:128
%     points = transformPoints(view, points, m3);
% end
% 
% title('Invert shear along z axis');
% for i = 1:128     
%     points = transformPoints(view, points, minv3);
% end

title('Invert shear along y axis');
for i = 1:100
    points = transformPoints(view, points, minv2);
end

title('Invert shear along x axis');
for i = 1:100
    points = transformPoints(view, points, minv1);
end

% title('Iteratively shear along all axes');
% for i = 1:64
%     points = transformPoints(view, points, m1);
%     points = transformPoints(view, points, m2);
%     points = transformPoints(view, points, m3);
% end
% 
% title('Revert previous shearing operations');
% for i = 1:64
%     points = transformPoints(view, points, minv3);
%     points = transformPoints(view, points, minv2);
%     points = transformPoints(view, points, minv1);
% end

%translate





disp('***** TRANSLATE *****');
m = matrix(4,'translation', 0.008, 0.007, 0.004)
minv = inv(m)
analyzeMatrix(m)

m2 = matrix(4,'translation', 0.0, 0.0, 0.007)
minv2 = inv(m2)
title('Translate');
for i = 1:128
    points = transformPoints(view, points, m);
end
title('Translate only one direction');
for i = 1:128
    points = transformPoints(view, points, m2);
end

r1 = matrix(3, 'rotation', 'x', 0.01);
r2 = matrix(3, 'rotation', 'y', 0.01);
r3 = matrix(3, 'rotation', 'z', 0.01);
r1inv = inv(r1);
r2inv = inv(r2);
r3inv = inv(r3);
title('Rotate around y axis');
for i = 1:200
    points = transformPoints(view, points, r2);
end

title('Rotate back');
for i = 1:200
    points = transformPoints(view, points, r2inv);
end

title('Translate back');
for i = 1:128
    points = transformPoints(view, points, minv2);
end
for i = 1:128
    points = transformPoints(view, points, minv);
end

m = matrix(4,'translation', 0.0, 0.0, 0.007)
title('Translate');
for i = 1:128
    points = transformPoints(view, points, m);
end
