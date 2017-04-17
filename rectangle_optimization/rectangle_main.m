% @file     rectangle_main.m
% @author   afruehstueck
% @date     03/04/2017

function rectangle_main(rows, cols,rect_w, rect_h, spacing_x, spacing_y, noise_rw, noise_rh, noise_sx, noise_sy, minimize_variables, alignment)
%     close all
%     clc;
%     clear;
    %% rectangle grid settings
    if nargin < 1
        %number of grid rows and columns
        rows = 5;       cols = 8;   
        %base width and heights of rectangles
        rect_w = 45;    rect_h = 50;       
        %base spacing between rectangles
        spacing_x = 25; spacing_y = 20;      
        %range of noise added to width and height [from -noise/2 to noise/2]        
        noise_rw = 30;  noise_rh = 30;       
        %range of noise added to spacing [from -noise/2 to noise/2]        
        noise_sx = 30;  noise_sy = 50;
    
        %% string of characters describing the rectangle constraints
        %c: column alignment
        %r: row alignment
        %w: width alignment 
        %h: height alignment
        %s: equal row spacing
        %t: equal column spacing
        minimize_variables = 'cr';
        %align everything:
        %minimize_variables = 'crwhst';

        %% string of two (!) characters describing the alignment in y and x direction
        %y-direction: c = center / t = top  / b = bottom alignment
        %x-direction: c = center / l = left / r = right alignment
        alignment = 'bl';     
    end
    
    num_rects = rows*cols;
    
    % use static color for all rects
%    edge_color = [0.27 0.00 0.58];
%    face_color = [0.96 0.95 0.98];
    
    % define color mapping for rectangles 
    % (darker color for edges, lighter color for faces)
    edge_colors = winter(num_rects);
    face_colors = edge_colors + (1 - edge_colors) * 0.8;
     
    %% create grid of rectangles with noise   
    %store [x y w h] of rectangles (center point, width, height)
    noisy_rectangles = zeros(num_rects, 4);
    for r = 1:rows
        for c = 1:cols
            n_w = randNoise(noise_rw);
            n_h = randNoise(noise_rh);
            w = rect_w + n_w; %create widths with noise
            h = rect_h + n_h; %create heights with noise
            
            x = (rect_w + spacing_x) * (c - 1) + randNoise(noise_sx);
            y = (rect_h + spacing_y) * (r - 1) + randNoise(noise_sy);

            noisy_rectangles((r - 1) * cols + c, :) = [x y w h];
        end
    end
    
    %calculate constrained rectangles    
    constrained_rectangles_objective = rectangle_optim_quadprog(noisy_rectangles, rows, cols, minimize_variables, alignment, true);
    constrained_rectangles_constrained = rectangle_optim_quadprog(noisy_rectangles, rows, cols, minimize_variables, alignment, false);
    constrained_rectangles_lsq = rectangle_optim_lsq(noisy_rectangles, rows, cols, minimize_variables, alignment);
    
    %% output results
    scr = get(0, 'ScreenSize'); 
    figsize = [80 120 5*scr(4)/3 2*scr(4)/3];
    figure('Name', 'Rectangle Constraints', 'NumberTitle', 'off', 'Position', figsize);
    
    %define common bounds for all plots
    outer_spacing = 2 * max(spacing_x, spacing_y);
    x_bounds = [-rect_w / 2 - outer_spacing, (cols - 1) * (rect_w + spacing_x) + rect_w / 2 + outer_spacing];
    y_bounds = [-rect_h / 2 - outer_spacing, (rows - 1) * (rect_h + spacing_y) + rect_h / 2 + outer_spacing];
    
    %draw original rectangles
    subplot(1, 4, 1);
    drawRectangles(noisy_rectangles, face_colors, edge_colors, x_bounds, y_bounds);
    title('rectangle grid with noise');
    
    %draw constrained rectangles
    subplot(1, 4, 2);
    drawRectangles(constrained_rectangles_objective, face_colors, edge_colors, x_bounds, y_bounds);
    title('constraints formulated as part of objective function');
    
    %draw constrained rectangles
    subplot(1, 4, 3);
    drawRectangles(constrained_rectangles_constrained, face_colors, edge_colors, x_bounds, y_bounds);
    title('quadprog with equality constraints');
    
    %draw constrained rectangles
    subplot(1, 4, 4);
    drawRectangles(constrained_rectangles_lsq, face_colors, edge_colors, x_bounds, y_bounds);
    title('lsqlin with equality constraints');
    
    %% output info/evaluation string
    constraintstr = '';
    for m=1:numel(minimize_variables)
        if m > 1  
            constraintstr = [constraintstr, ', '];
        end
        
        switch minimize_variables(m) 
            case 'c',    constraintstr = [constraintstr, 'column alignment'];
            case 'r',    constraintstr = [constraintstr, 'row alignment'];
            case 'w',    constraintstr = [constraintstr, 'equal widths'];
            case 'h',    constraintstr = [constraintstr, 'equal heights'];
            case 's',    constraintstr = [constraintstr, 'equal row spacing'];
            case 't',    constraintstr = [constraintstr, 'equal column spacing'];
        end
    end
    
    switch alignment(1) 
        case 't',    alignstr = 'top';
        case 'b',    alignstr = 'bottom';
        otherwise,   alignstr = 'center';
    end
    
    switch alignment(2)  
        case 'l',    alignstr = [alignstr, '/left']; 
        case 'r',    alignstr = [alignstr, '/right'];
        otherwise,   alignstr = [alignstr, '/center'];
    end
    fprintf('*********************************************************\n');
    fprintf(['Aligned %dx%d rectangle grid with ', alignstr, ' alignment.\nApplied constraints: ', constraintstr, '\n'], rows,  cols);
    fprintf('*********************************************************\n');
end


