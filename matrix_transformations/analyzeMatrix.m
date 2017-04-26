% @file     analyzeMatrix.m
% @author   afruehstueck
% @date     02/02/2017
%
% calculate several interesting matrix characteristics

function d = analyzeMatrix(matrix)
    d = det(matrix)
    [U, V] = eig(matrix)
    [V, J] = jordan(matrix)
    C = cond(matrix)
end