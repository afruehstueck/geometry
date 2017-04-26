% @file     clamp.m
% @author   afruehstueck
% @date     20/02/2017
%
% clamp values in M to range specified by minval and maxval
function C = clamp(M, minval, maxval) 
    if nargin == 1
        minval = 0;
        maxval = 1;
    end
    one_vec = ones(size(M));
    C = min(max(M, minval * one_vec), maxval * one_vec);
end