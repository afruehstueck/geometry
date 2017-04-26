% @file     drawSpectralGraph.m
% @author   afruehstueck
% @date     20/03/2017
%
% spectral graph drawing of graph with adjacency matrix A
% following Jean H. Gallier, CIS 515, Spectral Graph Drawing

function [V, Dvec, L] = drawSpectralGraph(A, labels)    
    %degree matrix (#edges for each vertex)
    D = diag(sum(A));
    %laplacian matrix (degree matrix - adjacency matrix)
    L = D - A;
    %find eigenvectors and eigenvalues of laplacian matrix
    [eVecs, ~] = eig(L);

    hold on
    %fix limits to make all limits equal and make room for matrix in plot
    xlim([-1.05, 1.05])
    ylim([-1.05, 1.05])
    %draw vertices
    gplot(A, eVecs(:, [3 2]), 'o')
    %draw edges
    gplot(A, eVecs(:, [3 2]))

    text(eVecs(:, 3) .* 1.1, eVecs(:, 2) .* 1.1, labels(1:size(A, 1)));
end

    