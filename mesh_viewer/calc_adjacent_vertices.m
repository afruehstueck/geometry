% @file     calc_adjacent_vertices.m
% @author   anna fruehstueck
% @date     18/02/2017
%
% calculate the adjacent vertices for each vertex
% calculate the adjacent vertices for each face

function [VxV, FxV] = calc_adjacent_vertices(V, F)
str = sprintf('generating vertex adjacency matrix...');
fprintf(str);

tic;
num_faces = size(F, 1);
vert_ind = [];
face_vert_ind = [];

for k=1:num_faces %iterate over all faces
    vertices = F(k, :);
    edges = nchoosek(vertices, 2); %look at all pairs of vertices
    vert_ind = [vert_ind; edges];
    
    face_vert_ind = [face_vert_ind; [repmat(k, length(vertices), 1), vertices']];
end

vert_ind = [vert_ind; fliplr(vert_ind)]; %symmetry
vert_ind = unique(vert_ind, 'rows'); %remove duplicates
entries = length(vert_ind)
toc;

VxV = sparse(vert_ind(1:entries, 1), vert_ind(1:entries, 2), ones(entries, 1), length(V), length(V));

face_vert_ind = unique(face_vert_ind, 'rows'); %remove duplicates
entries = length(face_vert_ind)
FxV = sparse(face_vert_ind(1:entries,1), face_vert_ind(1:entries,2), ones(entries, 1), length(F), length(V));
end
