% @file     ccvt.m
% @author   afruehstueck
% @date     27/03/2017
%
% implementation of Capacity Constrained Voronoi Tesselation (Balzer et al.)
% use image or custom function as distribution function for the points
% show iterative updates
% note: first iteration can be rather slow - following iterations are increasingly fast

function ccvt(N_sites, M_points, max_iterations, show_voronoi_regions, plot_in_two_figures, impath)    
    %% settings
    % number of seed points
    if nargin < 1
        N_sites = 4096;
        M_points = 64 * N_sites;
        
        % additional constraint to avoid iterating forever
        max_iterations = 50;
        
        %display the voronoi regions corresponding to the sites
        show_voronoi_regions = false;
        plot_in_two_figures = true;
        
        use_simple_distribution = false;
        use_exp_distribution = false;
        use_image = true;
        impath = '../data/img/face.jpg';
        impath = '../data/img/tiger.jpg';
        %impath = '../data/img/sloth_drawing.jpg';
    else
        use_simple_distribution = false;
        use_exp_distribution = false;
        use_image = true;
    end
    
    color_pts = [0.3, 0.3, 0.55];
    color_centroids = [0.27 0.00 0.58];
    color_lines = [0.4 0.4 0.6];
    line_width = 0.7;
    
    %% create random points with beta distribution
    if use_simple_distribution
        points = [random('Beta', 4, 4, M_points, 1) random('Beta', 4, 4, M_points, 1)];
    end
        
    %% create random points with distribution from paper
    if use_exp_distribution
        %create a grid of values
        sz_x = floor(M_points / 100); sz_y = sz_x;
        [xx, yy] = meshgrid(linspace(-1, 1, sz_x), linspace(-1, 1, sz_x));
        %evaluate function for each grid point
        pp = exp(1).^(-20.0 .* xx .* xx - 20.0 .* yy .* yy) + 0.2 .* sin(pi .* xx) .* sin(pi .* xx) .* sin(pi .* yy) .* sin(pi .* yy);
        %sample by rejecting points larger than rand
        rr = rand(size(xx));
        select = find(pp > rr);
    end
    
    %% create random points from image
    if use_image
        %read image
        img = imread(impath);
        if size(img, 3) == 3
            disp('convert image to grayscale...');
            img = rgb2gray(img); %convert img to grayscale
        end
        %convert back to double values
        img = double(img);

        sz_x = size(img, 2);
        sz_y = size(img, 1);

        %divide image by max value to get values between [0, 1]
        img = img / max(max(img));%mean(mean(img));

        %sample by rejecting points smaller than rand
        r = rand(sz_y, sz_x);
        select = find(img < r);    
    end
    
    %% subsample for too many selected points
    if use_image || use_exp_distribution
        if length(select) > M_points
            %select M_points from indices
            select = select(randsample(1:length(select), M_points));
        else
            modulo = mod(length(select), N_sites);
            M_points = length(select) - modulo;
            select = select(randsample(1:length(select), M_points));
        end
        %find coordinates of selected points
        [y, x] = ind2sub([sz_y, sz_x], select);
        %move coordinates to [0, 1] range    
        points = [1 - (y ./ sz_y), x ./ sz_x];
    end
    
    %%
    % display properties
    c_star = M_points / N_sites;
    fprintf('%d sites\n%d points\nc* = %d\n',N_sites, M_points, c_star);
    
    % seed random sites in plot
    sites = [rand(N_sites, 1) rand(N_sites, 1)];

    % create figure
    scr = get(0, 'ScreenSize'); 
    if plot_in_two_figures
        figsize = [100 150 scr(4) scr(4)/2];
    else
        figsize = [100 150 scr(4)/2 scr(4)/2];
    end
    
    figure('Name', 'CCVT', 'NumberTitle', 'off', 'Position', figsize);
    hold on;
    
    if plot_in_two_figures
        subplot(1, 2, 1)
    end
    
    %plot points extracted from distribution or image
    scatter(points(:, 2), points(:, 1), 4, 'MarkerFaceColor', color_pts, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .15)

    if plot_in_two_figures   
        %hide axes labels
        set(gca, 'XTick', []);
        set(gca, 'YTick', []);
        subplot(1, 2, 2)
        hold on;
    end
    
    if show_voronoi_regions
        [vx, vy] = voronoi(sites(:, 2), sites(:, 1));
        voronoi_plot = plot(sites(:, 2), sites(:, 1), 'ko', vx, vy, 'Color', color_lines, 'LineWidth', line_width, 'MarkerFaceColor', color_pts, 'MarkerEdgeColor', 'none');
    end
    % plot containing the sites
    sites_plot = plot(sites(:, 2), sites(:, 1), 'o', 'Color', color_centroids, 'MarkerSize', 2.0, 'LineWidth', 1.5);
    
    xlim([0 1]);
    ylim([0 1]);
    %hide axes labels
    set(gca, 'XTick', []);
    set(gca, 'YTick', []);
    drawnow;
    
    %initialize assignment of points to sites
    A = initialize_A(sites, points, c_star);
    
    %create lists containing all points per site
    [points_in_sites, radii_sites] = create_point_lists(sites, points, A);
    
    disp('Create vector of all possible combinations...');
    %use combinator function from mathworks to find all possible combinations of indices
    %which is way faster than combnk
    combinations = combinator(N_sites, 2, 'c', 'no');
    
    cur_iteration = 0;
    while(cur_iteration < max_iterations)  
        fprintf('************* CCVT iteration #%d *************\n', cur_iteration);
        [A, points_in_sites, radii_sites, num_iterations] = capacity_constrained_voronoi(sites, points, A, points_in_sites, radii_sites, combinations);

        if num_iterations == 0  
            fprintf('Algorithm converged after %d iterations.\n', cur_iteration);
            break;
        end

        for s = 1:N_sites
         points_in_current_site = points(points_in_sites{s}, :);
         sites(s, :) = mean(points_in_current_site, 1);
        end
        %update positions of sites to center of mass of contained points
        set(sites_plot, 'XData', sites(:, 2), 'YData', sites(:, 1));

        if show_voronoi_regions
         %update voronoi plot
         delete(voronoi_plot);
         [vx, vy] = voronoi(sites(:, 2), sites(:, 1));
         voronoi_plot = plot(sites(:, 2), sites(:, 1), 'ko', vx, vy, 'Color', color_lines, 'LineWidth', line_width, 'MarkerFaceColor', color_pts, 'MarkerEdgeColor', 'none');
        end

        drawnow;
        cur_iteration = cur_iteration + 1;
     end
end

%do a random initialization of the assignments of points to sites by assigning 
%points to the nearest neighboring site that still has available capacity
function A = initialize_A(sites, points, c_star)
    disp('Assign points to sites randomly as a start...');
    %do an initial assignment
    A = zeros(length(points), 1);
    current_capacities = zeros(length(sites), 1);

    all_assigned = false;
    iteration = 0;
    while ~all_assigned
        available = find(current_capacities < c_star);
        unassigned_points = find(A == 0);
        closest_sites = knnsearch(sites(available, :), points(unassigned_points, :));

        for pt=1:length(unassigned_points)
            current_point = unassigned_points(pt);
            current_closest_site = available(closest_sites(pt));
            current_capacity = current_capacities(current_closest_site);
            if(current_capacity < c_star)
                current_capacities(current_closest_site) = current_capacity + 1;
                A(current_point) = current_closest_site;
            end
        end
        iteration = iteration + 1;
        if all(A)
            all_assigned = true;
        end
    end
end

function [points_in_sites, radii_sites] = create_point_lists(sites, points, A)
    disp('Create point lists...');
    N_sites = length(sites);
    %create list of all points associated with each site
    points_in_sites = cell(N_sites, 1);
    radii_sites = zeros(N_sites, 1);
    for i=1:N_sites
        current_points = find(A == i);
        points_in_sites{i} = current_points;
        
        distances = distance(sites(i, :), points(current_points, :));
        radii_sites(i) = max(distances);
    end
end

function [A, points_in_sites, radii_sites, num_iterations] = capacity_constrained_voronoi(sites, points, A, points_in_sites, radii_sites, combinations)
    N_sites = length(sites);
    
    sites_stable = false(N_sites, 1);
    
    %do skipping calculations once and vectorized for speedup
    distances = distance(sites(combinations(:, 2), :), sites(combinations(:, 1), :));
    added_radii = radii_sites(combinations(:, 1)) + radii_sites(combinations(:, 2));
    skip_because_too_distant = distances > added_radii;
    
    iteration = 0;
    reverseStr = '';
   
    while ~all(sites_stable) 
        num_switches = zeros(N_sites, 1);
        
        %skip iteration if both of the sites are stable
        skip_because_both_stable = sites_stable(combinations(:, 1)) & sites_stable(combinations(:, 2));
        
        skip = skip_because_both_stable | skip_because_too_distant;
        current_combinations = combinations(~skip, :);
        for i=1:length(current_combinations)
            site_si = current_combinations(i, 1);
            site_sj = current_combinations(i, 2);
            

            si = sites(site_si, :);
            sj = sites(site_sj, :);
            
%            %%skipping criteria are now vectorized, which is less legible, but faster
%            %%skip iteration if sites lie far apart
%            if sqrt(sum((sj - si).^2)) > radii_sites(site_si) + radii_sites(site_sj)
%            %%skip iteration if both of the sites are stable           
%            if (sites_stable(site_si) == true && sites_stable(site_sj) == true) || skip_because_too_distant(i)
            
            points_xi = points_in_sites{site_si};
            xi = points(points_xi, :);
            
            points_xj = points_in_sites{site_sj};
            xj = points(points_xj, :);
            
            max_sq_radius = max(radii_sites(site_si)^2, radii_sites(site_sj)^2);
            
            dist_xi_sj = sq_distance(xi, sj);
            dist_xj_si = sq_distance(xj, si);
            
            %skipping criterion
            if all(dist_xi_sj > max_sq_radius) || all(dist_xj_si > max_sq_radius)
                continue;
            end
            
            %calculate distance values
            Hi = sq_distance(xi, si) - dist_xi_sj;
            Hj = sq_distance(xj, sj) - dist_xj_si;
            
            %sort distance values by size, sort point indices accordingly
            [Hi_sorted, sort_i] = sort(Hi);
            [Hj_sorted, sort_j] = sort(Hj);
            
            %sum distance values
            Hi_plus_Hj = Hi_sorted + Hj_sorted;
            
            %find indices of summed values that are greater than zero
            greater_zero = find(Hi_plus_Hj > 0);
            if numel(greater_zero) > 0 
                points_xi_sorted = points_xi(sort_i);
                points_xj_sorted = points_xj(sort_j);
            
                %switch assignments of sites for these points in A
                A(points_xi_sorted(greater_zero)) = site_sj;
                A(points_xj_sorted(greater_zero)) = site_si;
                
                remove_from_xi = points_xi_sorted(greater_zero);
                remove_from_xj = points_xj_sorted(greater_zero);
                
                %switch indices in point arrays
                points_xi_sorted(greater_zero) = remove_from_xj;
                points_xj_sorted(greater_zero) = remove_from_xi;
                
                %replace point arrays for site 
                points_in_sites{site_si} = points_xi_sorted;
                points_in_sites{site_sj} = points_xj_sorted;
                
                %update radii
                distances_si = distance(si, points(points_xi_sorted, :));
                radii_sites(site_si) = max(distances_si);
                
                distances_sj = distance(sj, points(points_xj_sorted, :));
                radii_sites(site_sj) = max(distances_sj);
                
                %count number of switches per sites
                num_switches(site_si) = num_switches(site_si) + 1;
                num_switches(site_sj) = num_switches(site_sj) + 1;
                
            end 
        end
        %site is counted as stable if there were no switches
        sites_stable(num_switches == 0) = true;
        
        %output progress, always printing to same line
        percent_stable = round(100 * numel(find(sites_stable == true)) / N_sites);
        msg = sprintf('Iteration %d: %3.0f percent of sites stable.', iteration, percent_stable); %Don't forget this semicolon
        fprintf([reverseStr, msg]);
        %clever trick: \b is the backspace character
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
   
        iteration = iteration + 1;
    end
    fprintf('\n');
    num_iterations = iteration - 1;
end

% helper function that calculates row-wise distance
function D = distance(M, N) 
    D = sqrt(sum((M - N).^2, 2));
end

% helper function that calculates squared row-wise distance
function D = sq_distance(M, N) 
    D = sum((M - N).^2, 2);
end