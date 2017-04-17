% @file     drawRectangles.m
% @author   afruehstueck
% @date     08/04/2017
% helper function that takes a list of rectangles in format [x y w h]
% where (x y) = coordinates of center point, (w h) = width&height of rectangle
% and draws the rectangles as well as their center points
function drawRectangles(rectangles, face_colors, edge_colors, x_bounds, y_bounds)
    hold on;
    %calculate bottom left of rectangle to use with MATLAB rectangle drawing function
    bottom_left_rectangles = rectangles;
    bottom_left_rectangles(:, 1) = rectangles(:, 1) - rectangles(:, 3)/2;
    bottom_left_rectangles(:, 2) = rectangles(:, 2) - rectangles(:, 4)/2;
    
    for rect = 1:size(rectangles, 1)
        %plot rectangle
        rectangle('Position', bottom_left_rectangles(rect, :), 'FaceColor', face_colors(rect, :), 'EdgeColor', edge_colors(rect, :));
        %plot rectangle center
        plot(rectangles(rect, 1), rectangles(rect, 2), '+', 'Color', edge_colors(rect, :));
    end
    
    xlim(x_bounds);
    ylim(y_bounds);
end
