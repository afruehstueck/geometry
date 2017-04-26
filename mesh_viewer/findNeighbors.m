% @file     findNeighbors.m
% @author   anna fruehstueck
% @date     25/02/2017
%
% generate data structure containing indices of all neighbors of vertex V

function [neighbors] = findNeighbors(V, F)
    neighbors = cell(1, size(V, 1));

    for i=1:length(F)
        neighbors{F(i,1)} = [neighbors{F(i,1)} [F(i,2) F(i,3)]];
        neighbors{F(i,2)} = [neighbors{F(i,2)} [F(i,3) F(i,1)]];
        neighbors{F(i,3)} = [neighbors{F(i,3)} [F(i,1) F(i,2)]];
    end

    for i=1:size(V,1)
        neighbors{i} = unique(neighbors{i});
        if isempty(neighbors{i})
            neighbors{i}=[];
        end
    end
end