% @file     rectangle_main.m
% @author   afruehstueck
% @date     03/04/2017

function rectangle_main(rows, cols)

    close all
    clc;
    clear;
    %% settings
    if nargin < 1
        rows = 5;
        cols = 6;
        
        spacing_x = 20;
        spacing_y = 30;
        noise_sx = 30;
        noise_sy = 50;
        
        rect_w = 30;
        rect_h = 50;
        noise_rw = 10; 
        noise_rh = 20;
        
    end
    
    num_rects = rows*cols;
    
    %static color for all rects
    %edge_color = [0.27 0.00 0.58];
    %face_color = [0.96 0.95 0.98];
    
    %colors = jet(10);
    %edge_colors = colors(randperm(rows*cols), :);
    %edge_colors = colors(randi(10, rows*cols, 1), :);
    edge_colors = winter(num_rects);
    face_colors = edge_colors + (1 - edge_colors) * 0.8;
    
    %store [x y w h] of rectangles (lower left point plus width and height)
    %
    rectangles = zeros(num_rects, 4);
    %create grid of rectangles with noise
    for r = 1:rows
        for c = 1:cols
            n_w = randNoise(noise_rw);
            n_h = randNoise(noise_rh);
            w = rect_w + n_w; %create widths with noise
            h = rect_h + n_h; %create heights with noise
            
            x = (rect_w + spacing_x) * (c - 1) + randNoise(noise_sx);% - n_w/2;
            y = (rect_h + spacing_y) * (r - 1) + randNoise(noise_sy);% - n_h/2;
            disp([num2str(x), ', ', num2str(y)]);
            
            rectangles((r - 1) * cols + c, :) = [x y w h];
        end
    end
    
    %calculate constrained rectangles
    constrained_rectangles_objective = rectangle_constrain(rectangles, rows, cols, true);
    constrained_rectangles_constrained = rectangle_constrain(rectangles, rows, cols, false);
    
    %% output results
    scr = get(0, 'ScreenSize'); 
    figsize = [100 150 2*scr(4) 3*scr(4)/4];
    figure('Name', 'Rectangle Constraints', 'NumberTitle', 'off', 'Position', figsize);
    %draw original rectangles
    subplot(1, 3, 1);
    hold on;
    
    for rect = 1:num_rects
        rectangle('Position', rectangles(rect, :), 'FaceColor', face_colors(rect, :), 'EdgeColor', edge_colors(rect, :));
    end
    xlim([min(rectangles(:, 1)) - spacing_x, max(rectangles(:, 1)) + rect_w + spacing_x]);
    ylim([min(rectangles(:, 2)) - spacing_y, max(rectangles(:, 2)) + rect_h + spacing_y]);
    
    %draw constrained rectangles
    subplot(1, 3, 2);
    hold on;
    
    for rect = 1:num_rects
        rectangle('Position', constrained_rectangles_objective(rect, :), 'FaceColor', face_colors(rect, :), 'EdgeColor', edge_colors(rect, :));
    end
    
    xlim([min(constrained_rectangles_objective(:, 1)) - spacing_x, max(constrained_rectangles_objective(:, 1)) + rect_w + spacing_x]);
    ylim([min(constrained_rectangles_objective(:, 2)) - spacing_y, max(constrained_rectangles_objective(:, 2)) + rect_h + spacing_y]);
    
    %draw constrained rectangles
    subplot(1, 3, 3);
    hold on;
    
    for rect = 1:num_rects
        rectangle('Position', constrained_rectangles_constrained(rect, :), 'FaceColor', face_colors(rect, :), 'EdgeColor', edge_colors(rect, :));
    end
    
    xlim([min(constrained_rectangles_constrained(:, 1)) - spacing_x, max(constrained_rectangles_constrained(:, 1)) + rect_w + spacing_x]);
    ylim([min(constrained_rectangles_constrained(:, 2)) - spacing_y, max(constrained_rectangles_constrained(:, 2)) + rect_h + spacing_y]);
end

% helper function that returns a random value between (-noise/2, noise/2)
function randVal = randNoise(noise)
    randVal = (rand * noise) - (noise / 2);
end