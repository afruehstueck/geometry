% @file     transformPoints.m
% @author   afruehstueck
% @date     30/01/2017

%applies matrix to points, updates viewer and returns transformed points
function transformed = transformPoints(viewer, points, matrix)

%clear window content
cla(viewer);
    
if size(points, 1) == 2 %2D matrix transformations
    if length(matrix) == 2
        transformed = (matrix * points);
    elseif length(matrix) == 3 %homogeneous transformations
        homogeneousPoints = [points; ones(1, length(points))];
        transformedHomogeneous = (matrix * homogeneousPoints);
        transformed = transformedHomogeneous(1:2, :);
    end
    %plot result of transformation in viewer
    plot(transformed(1, :), transformed(2, :), 'Color', [0.956, 0.258, 0.258]);
else %3D matrix transformations
    transformed = zeros(size(points));
    point = ones(length(matrix), 1);
    for x = 1:size(points, 1)
        for y = 1:size(points, 2)
            %extract point from n*m*3 matrix
            point(1:3) = squeeze(points(x, y, :));
            %transform point
            transformedPoint = matrix * point;
            if length(matrix) == 4
                transformedPoint = transformedPoint./transformedPoint(4);
            end
            transformed(x, y, :) = transformedPoint(1:3);
        end
    end
    %plot result of transformation in viewer
    surf(transformed(:, :, 1),transformed(:, :, 2),transformed(:, :, 3));
end

%animate
drawnow;
end