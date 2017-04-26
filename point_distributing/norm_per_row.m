% @file     norm_per_row.m
% @author   afruehstueck
% @date     20/02/2017
%
% calculate row-wise norm
function N = norm_per_row(M) 
    N = sqrt(sum(M.^2, 2));
end