% @file     calc_adjacent_vertices.m
% @author   anna fruehstueck
% @date     18/02/2017
function [S] = calc_adjacent_vertices(V, F)
str = sprintf('generating vertex adjacency matrix...');
fprintf(str);

tic;
num_rows = size(F, 1);
indices = [];

for k=1:num_rows %iterate over all faces
    vertices = F(k, :);
    edges = nchoosek(vertices, 2); %look at all pairs of vertices
    indices = [indices; edges];
end

indices = unique(indices, 'rows'); %remove duplicates
entries = length(indices)

toc;

S = sparse(indices(1:entries,1), indices(1:entries,2), ones(entries, 1), length(V), length(V));
end
