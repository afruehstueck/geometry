% @file     transformPoints.m
% @author   afruehstueck
% @date     30/01/2017

%applies matrix to points, updates viewer and returns transformed points
function transformed = transformPoints(points, matrix)
    if size(points, 1) == 2 %2D matrix transformations
        if length(matrix) == 2
            transformed = (matrix * points);
        elseif length(matrix) == 3 %homogeneous transformations
            homogeneousPoints = [points; ones(1, length(points))];
            transformedHomogeneous = (matrix * homogeneousPoints);
            transformed = transformedHomogeneous(1:2, :);
        end
    else %3D matrix transformations
        transformed = zeros(size(points));
        point = ones(length(matrix), 1);
        for i = 1:size(points, 2)
                %extract point from n*m*3 matrix
                point(1:3) = points(:, i);
                %transform point
                transformedPoint = matrix * point;
                if length(matrix) == 4
                    transformedPoint = transformedPoint./transformedPoint(4);
                end
                transformed(:, i) = transformedPoint(1:3);
        end
    end
end