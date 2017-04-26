% @file     polynomial.m
% @author   afruehstueck
% @date     06/02/2017
%
% polynomial solves f(x) = sum_{i+1:n)(coeff(i)*x^(i-1)) for all x
% degree of the polynomial is determined by the number of elements in coeff

function y = polynomial(x, coeff)
    y = zeros(1,length(x));
    for c=1:length(coeff) 
        y(:) = y(:) + coeff(c) * x(:).^(c-1);
    end
end