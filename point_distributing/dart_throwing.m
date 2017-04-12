% @file     dart_throwing.m
% @author   afruehstueck
% @date     24/03/2017
function dart_throwing(radius, grid_resolution, use_grid, pause_val)  
    %% settings
    %radius of disks
    if nargin < 1
        radius = 0.08;
        %use grid to enhance dart throwing (true/false)
        use_grid = true;
        %resolution of grid
        grid_resolution = 50;
        %pause after each throw for 'better' visualization (length of pause)
        pause_val = 0.0;
    end
    %%
    
    %colors used in plot
    purple = [0.27 0.00 0.58];
    yellow = [1.00 0.80 0.00];
    
    %radius size in grid coordinates
    radius_raster = radius * grid_resolution;
    %cell size in data coordinates
    d = 1 / grid_resolution;
    
    %create grid of zeroes for empty cells
    cell_grid = zeros(grid_resolution, grid_resolution);
    [xx, yy] = meshgrid(1:grid_resolution);
    centers = zeros(0, 2);

    % create figure
    scr = get(0, 'ScreenSize'); 
    figure('Name', 'Dart Throwing', 'NumberTitle', 'off', 'Position', [100 150 4*scr(4)/3 2*scr(4)/3]);
    
    %plot of empty/non-empty grid
    subplot(1, 2, 2);
    hold on;
    xlim([0.5 grid_resolution+.5]);
    ylim([0.5 grid_resolution+.5]);
    [xshow, yshow] = meshgrid(0.5:1:grid_resolution+0.5);
    matrix_plot = pcolor(xshow, yshow, padarray(cell_grid, [1 1], 'post'));
    %create color map for grid plot
    map = [1 1 1; purple .* (1:-0.2:0.6)'];
    colormap(map)
    caxis([0, 4])
    
    current_grid_plot = plot([0 0 0 0 0], [0 0 0 0 0], 'Color', yellow, 'LineWidth', 2);
    
    %plot of actual data
    subplot(1, 2, 1);
    xlim([0 1]);
    ylim([0 1]);
    hold on;
    current_cell_plot = plot([0 0 0 0 0], [0 0 0 0 0], 'Color', yellow, 'LineWidth', 2);
    centers_plot = scatter(centers(:, 1), centers(:, 2), 150, '.');
           
    t = 0:.1:2*pi;
    num_darts_thrown = 0;
    num_valid_darts = 0;
    all_cells_filled = false;
    while num_darts_thrown < 100000 && ~all_cells_filled
        valid_dart = false;
        while ~valid_dart && num_darts_thrown < 100000
            if use_grid
                %find all empty grid cells
                [empty_y, empty_x] = find(~cell_grid);
                if isempty(empty_x)
                    all_cells_filled = true;
                    break;
                end
                %pick a random cell by index
                random_idx = randi(length(empty_x));
                
                %grid coordinates of randomly selected cell
                x_center = empty_x(random_idx); %coordinate lies at center of grid cell
                y_center = empty_y(random_idx);
                fprintf('picked empty cell at (%d, %d)\n', x_center, y_center);
                
                %rectangle around picked cell in grid
                x_l = x_center - 0.5; x_r = x_center + 0.5;
                y_b = y_center - 0.5; y_t = y_center + 0.5;
                rect_grid = [x_l y_b; x_r y_b; x_r y_t; x_l y_t; x_l y_b];
                
                %visualize position of current cell in grid as rectangle
                set(current_grid_plot, 'XData', rect_grid(:, 1), 'YData', rect_grid(:, 2));
               
                %rectangle in original data space
                rect_data = (rect_grid - 0.5) / grid_resolution;
                %visualize position of current cell in data space
                set(current_cell_plot, 'XData', rect_data(:, 1), 'YData', rect_data(:, 2));
                
                %cell coordinates to original data space
                trnsf_x = (x_center - 1) / grid_resolution;
                trnsf_y = (y_center - 1) / grid_resolution;
                %pick a random dart in current cell region
                rand_dart_x = trnsf_x + rand * d;
                rand_dart_y = trnsf_y + rand * d;
                
                dart = [rand_dart_x rand_dart_y];
                pause(pause_val)
            else
                dart = [rand rand];
            end
            
            num_darts_thrown = num_darts_thrown + 1;
            %calculate distances of all circle centers from current dart position
            distances = norm_per_row(centers - repmat(dart,  num_valid_darts, 1));
            if all(distances >= radius)
                %accept dart position
                valid_dart = true;
                num_valid_darts = num_valid_darts + 1;
            else
                fprintf('invalid dart at [%1.2f, %1.2f]\n', dart(1), dart(2));
            end
        end
        
        if ~valid_dart
            continue;
        end
        
        %update current cell
        centers = [centers; dart];

        x = radius * cos(t) + dart(1);
        y = radius * sin(t) + dart(2);
        %plot new circle (as patch)
        patch('XData', x, 'YData', y, 'FaceColor', purple, 'FaceAlpha', 0.05, 'EdgeColor', purple);
        
        %transform dart coordinates back to raster to rasterize circle
        dart_raster = ([dart(1) dart(2)] .* grid_resolution) + 1;
        circle = sqrt((xx - dart_raster(1)).^2 + (yy - dart_raster(2)).^2) <= radius_raster;

        %add circle (ones in variable circle) to cell grid
        cell_grid = cell_grid + circle;
        
        %update matrix plot contents
        set(matrix_plot, 'CData', padarray(cell_grid, [1 1], 'post'));
        %update circle centers
        set(centers_plot, 'XData', centers(:, 1), 'YData', centers(:, 2));
          
        if isempty(find(~cell_grid, 1))
            all_cells_filled = true;
            break;
        end
        drawnow;
        pause(pause_val)
    end
    
    %remove yellow indicators for current cell when done
    set(current_cell_plot, 'XData', [-1 -1 -1 -1 -1], 'YData', [-1 -1 -1 -1 -1]);
    set(current_grid_plot, 'XData', [-1 -1 -1 -1 -1], 'YData', [-1 -1 -1 -1 -1]);
               
    fprintf('*****************************\n %d darts thrown\n %d darts chosen\n', num_darts_thrown, num_valid_darts);
end


% helper function that calculates row-wise norm
function N = norm_per_row(M) 
    N = sqrt(sum(M.^2, 2));
end

% helper function that clamps values to range
function C = clamp(M, minval, maxval) 
    if nargin == 1
        minval = 0;
        maxval = 1;
    end
    one_vec = ones(size(M));
    C = min(max(M, minval * one_vec), maxval * one_vec);
end
