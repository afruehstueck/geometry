% @file     applyMatrices.m
% @author   afruehstueck
% @date     10/03/2017
%
% apply Laplacian and mass matrix to data
% return new vertex positions, laplacian matrix and vector of update distances of laplace-beltrami
function [V, Dvec, L] = applyMatrices(V, M, D, lambda)
    L = D * M;
    Dvec = lambda * L * V;
    V = V + Dvec;
end