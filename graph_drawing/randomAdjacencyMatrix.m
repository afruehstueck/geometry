% @file     randomAdjacencyMatrix.m
% @author   afruehstueck
% @date     20/03/2017
%
% generate a random adjacency matrix of size NxN
% matrix should be symmetric and have zeros on the main diagonal
% generates adjacency matrix that has at least degree 2 for each vertex

function [A] = randomAdjacencyMatrix(N)    
    %generate NxN matrix of random 0s and 1s
    G = round(rand(N));
    %take only upper triangular part and replicate on lower triangular
    %while leaving the main diagonal as zeros (no loops)
    randA = triu(G, 1) + triu(G, 1)';

    %check if each vertex has at least one edge
    if any(sum(randA) < 2)
        disp('<2 edges for at least one vertex... generate new matrix');
        A = randomAdjacencyMatrix(N);  
    else
        A = randA;
    end
end
