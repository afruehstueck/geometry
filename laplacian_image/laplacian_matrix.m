% @file     laplacian_matrix.m
% @author   anna fruehstueck
% @date     01/03/2017
%
% generate a laplacian matrix for an image with dimensions sx and sy

function L = laplacian_matrix(sx, sy) 
    num_px = sy * sx;

    %main diagonal is all -4s
    diag0 = ones(num_px, 1)*(-4);
    
    diag1l = repmat([ones(sy-1, 1); 0], sx, 1);
    diag1r = repmat([0; ones(sy-1, 1)], sx, 1);

    diag2 = ones(num_px, 1);
    L = spdiags([diag2 diag1l diag0 diag1r diag2], [-sy; -1; 0; 1; sy], num_px, num_px);
end