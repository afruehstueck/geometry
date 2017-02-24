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
        [adj_vertices, ~] = find(VxV(:, i));
        num_verts = size(adj_vertices, 1);
        
        accum = [ 0 0 0 ];
        for f=1:num_verts
            j = adj_vertices(f);
            v_j = V(j,:);
            
            d_ji = v_j - v_i; 
            accum = accum + d_ji;  
        end
        L(i,:) = accum / num_verts;
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
        [adj_face_indices, ~] = find(FxV(:, i));
        adj_faces = F(adj_face_indices, :);
        
        %adj_faces = sort(adj_faces, 2)
        %adj_faces = unique(adj_faces, 'rows') %remove duplicates
        
        [adj_vertices, ~] = find(VxV(:, i));
        num_adj_verts = size(adj_vertices, 1);
        
        areas = 0.0;
        accum = [ 0 0 0 ];
        for f=1:num_adj_verts
            j = adj_vertices(f);
            v_j = V(j,:);
            
            %      v_i o -------o v_k2
            %         /  \     /
            %        /    \   / 
            %       /      \ /
            % v_k1 o-------o v_j
            %
            %find third vertices for faces adjacent to edge [v_i, v_j]
            [rows, ~] = find(adj_faces == j);
            triangles = adj_faces(rows, :);
            lin = find(triangles~=i & triangles~= j); %find adjacent vertices to [v_i, v_j]
            v_inds = triangles(lin); %get indices 
            
            if(length(v_inds) ~= 2) 
                disp(length(v_inds))
            end
        
            v_k = V(v_inds,:); %get coordinates
            
            ik = v_i - v_k;
            jk = v_j - v_k;
            
            dots = dot(ik, jk, 2); %dot products from vectors in rows
            lengths_i = (sqrt(sum((ik').^2)))'; %row-wise norm ||v_i - v_k||
            lengths_j = (sqrt(sum((jk').^2)))'; %row-wise norm ||v_j - v_k||
            comb_lenghts = lengths_i .* lengths_j; %||v_i - v_k|| * ||v_j - v_k||
            diffs = dots ./ comb_lenghts;
            angles = acos(diffs); %angles at vertices v_k
            cotan_sum_angles = sum(cot(angles));
            
            %TODO currently only considering triangle1
            %cotangent weight scheme [Meyer '02]
            if true%is_obtuse([v_i; v_j; v_k(1,:)])
                %TODO wrong
                areas = areas + sum ((0.5 * comb_lenghts .* sin(angles)) / 3); %third of area of triangles
            else
                length_ij = norm(v_i - v_j);
                areas = areas + cotan_sum_angles * length_ij * length_ij;
            end
            
            ji = v_j - v_i;
            accum = accum + cotan_sum_angles * ji;
%             verts = F(adj_faces(f), :);
%             v1 = V(verts(1),:);
%             v2 = V(verts(2),:);
%             v3 = V(verts(3),:);
%             c = cross((v2 - v1), (v3 - v1));
%             Ns = Ns + c;
            %Ns(f, :) = c;% / norm(c);   
        end
        %TODO
        A = sum(areas); %took each triangle area twice, cancels out 2* in equation
        %A = 3 / sum(areas);
        L(i,:) = accum / A;
    end
toc;
end

%see http://mathworld.wolfram.com/ObtuseTriangle.html
function result = is_obtuse(vertices)
    ind = [1 2; 1 3; 2 3]; %all combinations of indices
    rem = [3; 2; 1]; %remaining point
    edges = vertices(ind(:,2),:) - vertices(ind(:,1),:);
    edgessq = sum(abs(edges).^2,2); %row-wise norm squared
    compare = (edgessq(ind(:,1),:) + edgessq(ind(:,2),:) < edgessq(rem(:,1),:)); %evaluate all a^2 + b^2 < c^2
    result = any(compare); %is any of the inequalities true
end