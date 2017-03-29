% @file     ccvt.m
% @author   afruehstueck
% @date     27/03/2017

function ccvt(N_sites, M_points, max_iterations, max_error)
    close all
    clc;
    clear;
    
    %% settings
    % number of seed points
    if nargin < 1
        N_sites = 128;
        M_points = 256 * N_sites;
        
        % iterate until error gets small
        max_error = 1e-8;
        % additional constraint to avoid iterating forever
        max_iterations = 30;
    end
    
    %% create random points with beta distribution
    %points = [random('Beta', 4, 4, M_points, 1) random('Beta', 4, 4, M_points, 1)];
       
    %% create random points from image
    impath = '../data/img/bunny_drawing.jpg';
    impath = '../data/img/sloth_drawing.jpg';
    %impath = '../data/img/lena.jpg';
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

    %divide image by mean value to get values between [0, 1]
    img = img / mean(mean(img));
    
    %sample by rejecting points smaller than rand
    r = rand(sz_y, sz_x);
    select = find(img < r);    

    %% create random points with distribution from paper
    %create a grid of values
    sz_x = floor(M_points / 10); sz_y = sz_x;
    [xx, yy] = meshgrid(linspace(-1, 1, sz_x), linspace(-1, 1, sz_x));
    %evaluate function for each grid point
    pp = exp(1).^(-20.0 .* xx .* xx - 20.0 .* yy .* yy) + 0.2 .* sin(pi .* xx) .* sin(pi .* xx) .* sin(pi .* yy) .* sin(pi .* yy);
    %sample by rejecting points larger than rand
    rr = rand(size(xx));
    select = find(pp > rr);
    
    %% subsample for too many selected points
    if length(select) > M_points
        %select M_points from indices
        select = select(randsample(1:length(select), M_points));
    else
        modulo = mod(length(select), N_sites);
        select(end-modulo+1:end, :) = [];
        M_points = length(select);
    end
    %find coordinates of selected points
    [y, x] = ind2sub([sz_y, sz_x], select);
    %move coordinates to [0, 1] range    
    points = [1 - (y ./ sz_y), x ./ sz_x];
    
    %%
    % seed random sites in plot
    sites = [rand(N_sites, 1) rand(N_sites, 1)];

    c_star = M_points / N_sites;
    disp(['c* = ', num2str(c_star)]);
    
    % create figure
    scr = get(0, 'ScreenSize'); 
    figure('Name', 'CCVT', 'NumberTitle', 'off', 'Position', [100 150 2*scr(4)/3 2*scr(4)/3]);
    hold on;
    
    % polygon plots for each cell
    cells_plots = zeros(N_sites, 1);
    for idx = 1:N_sites
        cells_plots(idx) = patch(sites(idx, 2), sites(idx, 1), 'red', 'FaceAlpha', 1.0);
    end
    
    hold on;
    color_pts = [0.4, 0.3, 0.95];
    color_centroids = [0.27 0.00 0.58];
    color_lines = [0.4 0.4 0.6];
    line_width = 0.7;
    
    scatter(points(:, 2), points(:, 1), 4, 'MarkerFaceColor', color_pts, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .2)
    %pts_plot = scatter(sites(:, 2), sites(:, 1), 300, 'r.');
    
    [vx, vy] = voronoi(sites(:, 2), sites(:, 1));
    voronoi_plot = plot(sites(:, 2), sites(:, 1), 'ko', vx, vy, 'Color', color_lines, 'LineWidth', line_width, 'MarkerFaceColor', color_pts, 'MarkerEdgeColor', 'none');

    % plot containing the polygon centroids
    sites_plot = plot(sites(:, 2), sites(:, 1), 'o', 'Color', color_centroids, 'MarkerSize', 5, 'LineWidth', 1.5);
    
    xlim([0 1]);
    ylim([0 1]);
    set(gca, 'XTick', []);
    set(gca, 'YTick', []);
    drawnow;
    
    
%     for pt = 1:M_points
%         available = find(current_capacities < c_star);
%         idx = knnsearch(sites(available, :), points(pt, :));
%         A(pt) = available(idx);
%         current_capacities(available(idx)) = current_capacities(available(idx)) + 1;
%     end
    
    A = initialize_A(sites, points, c_star);
    
    cur_iteration = 0;
    while(cur_iteration < max_iterations)  
         [A, points_in_sites] = capacity_constrained_voronoi(sites, points, A);
         
         for s = 1:N_sites
             points_in_current_site = points(points_in_sites{s}, :);
             sites(s, :) = mean(points_in_current_site, 1);
         end
         set(sites_plot, 'XData', sites(:, 2), 'YData', sites(:, 1));
         %voronoi(sites(:, 2), sites(:, 1))
         
         delete(voronoi_plot)
         [vx, vy] = voronoi(sites(:, 2), sites(:, 1));
         voronoi_plot = plot(sites(:, 2), sites(:, 1), 'ko', vx, vy, 'Color', color_lines, 'LineWidth', line_width, 'MarkerFaceColor', color_pts, 'MarkerEdgeColor', 'none');
         drawnow;
         
         cur_iteration = cur_iteration + 1;
     end
end

%do a random initialization of the assignments of points to sites by assigning 
%points to the nearest neighboring site that still has available capacity
function A = initialize_A(sites, points, c_star)
    disp('Assign points to sites randomly as a start');
    %do an initial assignment
    tic;

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
    toc;
end

function [A, points_in_sites] = capacity_constrained_voronoi(sites, points, A)
    N_sites = length(sites);
    stable = false;
    
    %create list of all points associated with each site
    points_in_sites = cell(N_sites, 1);
    sites_stable = zeros(N_sites, 1);
    for i=1:N_sites
        points_in_sites{i} = find(A == i);
    end
        
    
    combinations = combnk(1:N_sites, 2);
    iteration = 0;
    while ~all(sites_stable) 
        %sites_stable = ones(N_sites, 1);
        num_switches = zeros(N_sites, 1);
        
        for i=1:length(combinations)
            site_si = combinations(i, 1);
            site_sj = combinations(i, 2);
            si = sites(site_si, :);
            sj = sites(site_sj, :);
            
            if sites_stable(site_si) == true && sites_stable(site_sj) == true
                continue;
            end
            
            points_xi = points_in_sites{site_si};
            xi = points(points_xi, :);
            
            points_xj = points_in_sites{site_sj};
            xj = points(points_xj, :);
            
            %calculate distance values
            Hi = sum((xi - si).^2, 2) - sum((xi - sj).^2, 2);
            Hj = sum((xj - sj).^2, 2) - sum((xj - si).^2, 2);
            
            %sort distance values by size, sort point indices accordingly
            [Hi_sorted, sorti] = sort(Hi);
            points_xi_sorted = points_xi(sorti);
            
            [Hj_sorted, sortj] = sort(Hj);
            points_xj_sorted = points_xj(sortj);
            
            %sum distance values
            Hi_plus_Hj = Hi_sorted + Hj_sorted;
            
            %find indices of summed values that are greater than zero
            greater_zero = find(Hi_plus_Hj > 0);
            if numel(greater_zero) > 0 
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
                
                %set stability for sites that had switches
                %sites_stable(site_si) = false;
                %sites_stable(site_sj) = false;
                num_switches(site_si) = num_switches(site_si) + 1;
                num_switches(site_sj) = num_switches(site_sj) + 1;
            end 
        end
        
        sites_stable(num_switches == 0) = true;
        iteration = iteration + 1;
        disp(['Iteration ', num2str(iteration), ': ', num2str(numel(find(sites_stable > 0))), '/' num2str(N_sites), ' sites stable.']);
    end
end