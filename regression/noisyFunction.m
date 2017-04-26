% @file     noisyFunction.m
% @author   afruehstueck
% @date     06/02/2017
%
% noisyFunction evaluates the function in function handle fun (must be
% passed using '@' operator, e.g. @sin) and adds some noise to the output
% values in y (range of the noise is determined by range

function [y, points] = noisyFunction(x, range, fun, coeff)
    if exist('coeff', 'var') 
        y = fun(x, coeff); %if variable coeff exists, pass coeff to fun
    else
        y = fun(x);
    end

    % scale random values to noise range e.g. [-0.1, +0.1] for range 0.2
    noise = -(range / 2) + range * rand(size(y)); 
    %points(:) = y(:) + noise(:);
    points = y + noise;
end
