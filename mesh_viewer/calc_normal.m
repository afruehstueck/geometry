% @file     calc_normal.m
% @author   anna fruehstueck
% @date     21/02/2017

function [N] = calc_normal(V, F, FxV)
    
% % Face normals
%     num_faces = size(F, 1);
%     N = zeros(num_faces, 3);
%     for k=1:num_faces-1 %iterate over all faces
%          vert_ind = F(k, :);
%          v1 = V(vert_ind(1),:)
%          v2 = V(vert_ind(2),:)
%          v3 = V(vert_ind(3),:)
%          
%          c = cross((v2 - v1), (v3 - v1));
%          N(k, :) = c / norm(c);
%     end
tic;
% Vertex normals
    num_vertices = size(V, 1);
    N = zeros(num_vertices, 3);
    
    for k=1:num_vertices %iterate over all vertices
        [adj_faces,~] = find(FxV(:, k));
        
        num_faces = size(adj_faces, 1);
        Ns = zeros(1, 3);
        %Ns = zeros(num_faces, 3);
        %ws = zeros(num_faces, 1);
        for f=1:num_faces
            verts = F(adj_faces(f), :);
            v1 = V(verts(1),:);
            v2 = V(verts(2),:);
            v3 = V(verts(3),:);
            c = cross((v2 - v1), (v3 - v1));
            Ns = Ns + c;
            %Ns(f, :) = c;% / norm(c);   
        end
        N(k, :) = Ns / norm(Ns);
    end
toc;
end