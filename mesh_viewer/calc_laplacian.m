% @file     calc_laplacian.m
% @author   anna fruehstueck
% @date     21/02/2017

function [L] = calc_laplacian(type, V, F, FxV, VxV) 
    if strcmp(type, 'uniform')
        L = calc_uniform_laplacian(V, F, FxV, VxV);
    elseif strcmp(type, 'cotan')
        L = calc_cotan_laplacian(V, F, FxV, VxV);
    end
end

function [L] = calc_uniform_laplacian(V, F, FxV, VxV) 
tic;
% Vertex normals
    num_vertices = size(V, 1);
    L = zeros(num_vertices, 3);
    
    for i=1:num_vertices %iterate over all vertices      
        v_i = V(i,:);
        [adj_face_indices,~] = find(FxV(:, i));
        adj_faces = F(adj_face_indices, :);
  
        [adj_vertices, ~] = find(VxV(:, i));
        num_verts = size(adj_vertices, 1);
        
        accum = [ 0 0 0 ];
        for f=1:num_verts
            j = adj_vertices(f);
            v_j = V(j,:);
            
            d_ji = v_j - v_i; 
            accum = accum + d_ji;
%             verts = F(adj_faces(f), :);
%             v1 = V(verts(1),:);
%             v2 = V(verts(2),:);
%             v3 = V(verts(3),:);
%             c = cross((v2 - v1), (v3 - v1));
%             Ns = Ns + c;
            %Ns(f, :) = c;% / norm(c);   
        end
        %TODO
        A = 3 / sum(areas);
        L(i,:) = A * accum;
    end
toc;
end

function [L] = calc_cotan_laplacian(V, F, FxV, VxV) 
tic;
% Vertex normals
    num_vertices = size(V, 1);
    L = zeros(num_vertices, 3);
    
    for i=1:num_vertices %iterate over all vertices      
        v_i = V(i,:);
        [adj_face_indices,~] = find(FxV(:, i));
        adj_faces = F(adj_face_indices, :);
  
        [adj_vertices, ~] = find(VxV(:, i));
        num_verts = size(adj_vertices, 1);
        
        areas = 0.0;
        accum = [ 0 0 0 ];
        for f=1:num_verts
            j = adj_vertices(f);
            v_j = V(j,:);
            
            %      v_i o -------o v_k2
            %         /  \     /
            %        /    \   / 
            %       /      \ /
            % v_k1 o--------o v_j
            %
            %find third vertices for faces adjacent to edge [v_i, v_j]
            [rows, ~] = find(adj_faces == j);
            triangles = adj_faces(rows, :);
            lin = find(triangles~=i & triangles~= j); %find adjacent vertices to [v_i, v_j]
            v_inds = triangles(lin); %get indices 
            v_k = V(v_inds,:); %get coordinates
            
            d_ik = v_i - v_k;
            d_jk = v_j - v_k;
            d_ij = v_i - v_j;
            dots = dot(d_ik, d_jk, 2); %dot products from vectors in rows
            lengths_i = (sqrt(sum((d_ik').^2)))'; %||v_i - v_k||
            lengths_j = (sqrt(sum((d_jk').^2)))'; %||v_j - v_k||
            lens = lengths_i .* lengths_j; %||v_i - v_k|| * ||v_j - v_k||
            diffs = dots ./ lens;
            angles = acos(diffs); %angles at vertices v_k
            
            %TODO 
            areas = areas + (0.5 * lens .* sin(angles)); %areas of triangles
 
            accum(:) = accum + sum(cot(angles)) * d_ij;
%             verts = F(adj_faces(f), :);
%             v1 = V(verts(1),:);
%             v2 = V(verts(2),:);
%             v3 = V(verts(3),:);
%             c = cross((v2 - v1), (v3 - v1));
%             Ns = Ns + c;
            %Ns(f, :) = c;% / norm(c);   
        end
        %TODO
        A = 3 / sum(areas);
        L(i,:) = A * accum;
    end
toc;
end