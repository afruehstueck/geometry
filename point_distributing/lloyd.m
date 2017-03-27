% @file     lloyd.m
% @author   afruehstueck
% @date     23/03/2017

function lloyd(N, max_iterations, max_error, seed_from_center)
    %% settings
    % number of seed points
    if nargin < 1
        N = 128;
        seed_from_center = false;

        % iterate until error gets small
        max_error = 1e-4;
        % additional constraint to avoid iterating forever
        max_iterations = 500;
    end
    
    %%
    %initialize random point positions    
    if seed_from_center
        % seed from center
        pts_x = 0.4 + rand(N, 1) * 0.2;
        pts_y = 0.4 + rand(N, 1) * 0.2;
    else
        % seed in entire plot
        pts_x = rand(N, 1);
        pts_y = rand(N, 1);
    end
        
    cts_x = ones(N, 1);
    cts_y = ones(N, 1);
    
    % bounding box
    bounds = [ 0, 0; 0, 1; 1, 1; 1, 0 ];
    
    % create figure
    scr = get(0, 'ScreenSize'); 
    figure('Name', 'Lloyd''s Algorithm', 'NumberTitle', 'off', 'Position', [100 150 2*scr(4)/3 2*scr(4)/3]);
    hold on;
    
    % assign random colors from colormap to cells
    colors = viridis(N);
    
    % polygon plots for each cell
    cells_plots = zeros(N, 1);
    for idx = 1:N
        cells_plots(idx) = patch(pts_x(idx), pts_y(idx), colors(idx, :));
    end
    
    % initialize error measure
    errors = ones(N, 1);
    
    % plot containing the current points - colormap shows the current error
    % for each point: the more yellow it gets, the smaller the error
    pts_plot = scatter(pts_x, pts_y, 300, 'r.');
    colormap(flipud(autumn))
    caxis([1e-5 1e-3])
    
    % plot containing the polygon centroids
    cts_plot = plot(pts_x, pts_y, 'wo', 'MarkerSize', 8, 'LineWidth', 2);
    
    cur_iteration = 0;
    while(max(errors) > max_error && cur_iteration < max_iterations)   
        cur_iteration = cur_iteration + 1;
        
        % use helper function to calculate bounded voronoi regions
        [V, C] = voronoiBounded(pts_x, pts_y, bounds);

        for idx = 1:numel(C)
            verts_x = V(C{idx}, 1);
            verts_y = V(C{idx}, 2);
            [cts_x(idx), cts_y(idx)] = centroid(verts_x, verts_y);
        
            % update the voronoi cells in plot
            set(cells_plots(idx), 'XData', verts_x, 'YData', verts_y);
        end
        
        % calculate euclidean distances between point positions and new centroid positions
        errors = norm_per_row([cts_x cts_y] - [pts_x pts_y]);
        disp(['current maximum error: ', num2str(max(errors))]);
        
        % move all points to new centroid positions (unless these are nan)
        pts_x = ~isnan(cts_x) .* cts_x + isnan(cts_x) .* pts_x;
        pts_y = ~isnan(cts_y) .* cts_y + isnan(cts_y) .* pts_y;
        
        % update centroid positions in plot
        set(cts_plot, 'XData', cts_x, 'YData', cts_y);
        
        % redraw all plot contents during animation
        drawnow;
        
        % update point positions in plot
        set(pts_plot, 'XData', pts_x, 'YData', pts_y, 'CData', errors);
    end
    
    % output evaluation
    if cur_iteration == max_iterations
        disp(['suspended algorithm after ', num2str(max_iterations), ' iterations']);
    else
        disp(['algorithm converged after ', num2str(cur_iteration), ' iterations']);
    end
end

% calculate coordinates of centroid of polygon from polygon vertices vX/
% according to https://en.wikipedia.org/wiki/Centroid#Centroid_of_polygon
function [cX, cY] = centroid(vX, vY)
    vXp1 = vX([2:end 1]); %x coordinates with shifted indices (plus 1)
    vYp1 = vY([2:end 1]); %y coordinates with shifted indices (plus 1)

    %term used multiple times in calculation
    multiplier = vX.*vYp1 - vXp1.*vY;

    A = sum(multiplier) / 2; %signed area of the polygon

    cX = sum((vX + vXp1) .* multiplier) / (6 * A);
    cY = sum((vY + vYp1) .* multiplier) / (6 * A);
end

% http://www.mathworks.com/matlabcentral/fileexchange/34428-voronoilimit
function [V, C] = voronoiBounded(x, y, bounds)
    % add 4 additional edges
    xA = [x; 0;  0; -2; 2];
    yA = [y; -2; 2;  0; 0];

    [vi,ci] = voronoin([xA, yA]);

    % remove the last 4 cells
    C = ci(1:end-4);
    V = vi;
    
    for idx=1:length(C)
            % convert the contour coordinate to clockwise order
            [X2, Y2] = poly2cw(V(C{idx}, 1), V(C{idx}, 2));
            % polybool restricts polygons to domain
            [vert_x, vert_y] = polybool('intersection', bounds(:, 1), bounds(:, 2), X2, Y2);
            inds = nan(1, length(vert_x));
            for vert_i = 1:length(vert_x)
                if any(V(:, 1) == vert_x(vert_i)) && any(V(:, 2) == vert_y(vert_i))
                    i_x = find(V(:, 1) == vert_x(vert_i));
                    i_y = find(V(:, 2) == vert_y(vert_i));
                    for ib=1:length(i_x)
                        if any(i_x(ib) == i_y)
                            inds(vert_i) = i_x(ib);
                        end
                    end
                    if isnan(inds(vert_i)) == 1
                        lv = length(V);
                        V(lv+1,1) = vert_x(vert_i);
                        V(lv+1,2) = vert_y(vert_i);
                        inds(vert_i) = lv + 1;
                    end
                else
                    lv = length(V);
                    V(lv+1, 1) = vert_x(vert_i);
                    V(lv+1, 2) = vert_y(vert_i);
                    inds(vert_i) = lv + 1;
                end
            end
            C{idx} = inds;
    end
end

% helper function that calculates row-wise norm
function N = norm_per_row(M) 
    N = sqrt(sum(M.^2, 2));
end