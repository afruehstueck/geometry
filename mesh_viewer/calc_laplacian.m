% @file     calc_laplacian.m
% @author   anna fruehstueck
% @date     21/02/2017

function [L, W, K] = calc_laplacian(type, V, F, FxV, VxV) 
    K = 0;
    W = 0;
    if strcmp(type, 'uniform_vector')
        [L] = calc_uniform_laplacian_vector(V, F, FxV, VxV);
    elseif strcmp(type, 'uniform')
        [L, W] = calc_uniform_laplacian(V, F, FxV, VxV);
    elseif strcmp(type, 'cotan')
        [L, W] = calc_cotan_laplacian(V, F, FxV, VxV);
    elseif strcmp(type, 'cotan_vector')
        [L, K] = calc_cotan_laplacian_vector(V, F, FxV, VxV);
    end
end

function [L, W] = calc_uniform_laplacian(V, F, FxV, VxV) 
    num_vertices = size(VxV, 1);
    sp_L = [];
    sp_W = zeros(num_vertices, 1);
    for i=1:num_vertices %iterate over all vertices      
        [adj_verts, ~] = find(VxV(:, i));
        num_adj_verts = numel(adj_verts);
      
        sp_W(i) = 1/num_adj_verts;
        sp_L = [sp_L; [i i -num_adj_verts]];
        sp_L = [sp_L; [repmat(i, num_adj_verts, 1) adj_verts ones(num_adj_verts, 1)]];
    end
    
    L = sparse(sp_L(:, 1), sp_L(:, 2), sp_L(:, 3));
    W = spdiags(sp_W, 0, num_vertices, num_vertices);
end


function [L, W] = calc_cotan_laplacian(V, F, FxV, VxV) 
    num_vertices = size(VxV, 1);

    sp_L = [];
    sp_W = zeros(num_vertices, 1);
    
    for i=1:num_vertices %iterate over all vertices  
        vi = V(i, :);
        [adj_verts, ~] = find(VxV(:, i));
        num_adj_verts = size(adj_verts, 1);
         
        [adj_faces, ~] = find(FxV(:, i));
        adj_triangles = F(adj_faces, :);
        num_adj_triangles = size(adj_triangles, 1);
        
        accum_area = 0;
        accum_voronoi = 0;
        accum_angle = 0;
        
%         for f=1:num_adj_triangles
%             cur_triangle = adj_triangles(f, :);
%             cur_triangle = cur_triangle(cur_triangle~=i); %remove center point
%             o_verts = V(cur_triangle, :);
%             vec1 = o_verts(1, :) - vi;
%             vec2 = o_verts(2, :) - vi;
%             
%             %angle in vi
%             angle_vi = acos(dot(vec1 / norm(vec1), vec2 / norm(vec2)));
%             %triangle area
%             area = norm(cross(vec1, vec2))/2;
%             
%             vec1_o1 = -vec1;
%             vec2_o1 = o_verts(2, :) - o_verts(1, :);
%             angle_o1 = acos(dot(vec1_o1 / norm(vec1_o1), vec2_o1 / norm(vec2_o1), 2));
%             
%             vec1_o2 = -vec2;
%             vec2_o2 = -vec2_o1;
%             angle_o2 = acos(dot(vec1_o2 / norm(vec1_o2), vec2_o2 / norm(vec2_o2), 2));
%             
%             %voronoi area
%             voronoi_area = (1/8)*( norm(vec1)^2/tan(angle_o2) + norm(vec2)^2/tan(angle_o1) ); 
%             
%             accum_voronoi = accum_voronoi + voronoi_area;
%             accum_angle = accum_angle + angle_vi;
%             accum_area = accum_area + area / 3;
%         end
                
        %step through neighborhood of v_i
        sum_w = 0;
        for f=1:num_adj_verts
            j = adj_verts(f);
            vj = V(j, :); 
            
            [rows, ~] = find(adj_triangles == j);
            cur_triangles = adj_triangles(rows, :);
            ind_ks = cur_triangles(cur_triangles~=j & cur_triangles~=i);
            vks = V(ind_ks, :);
            ki = vi - vks;
            kj = vj - vks;
            
            %if is_obtuse([vi, vj, 
%             dot1 = dot(ki(1,:) / norm(ki(1,:)), kj(1,:) / norm(kj(1,:)), 2);
%             dot2 = dot(ki(2,:) / norm(ki(2,:)), kj(2,:) / norm(kj(2,:)), 2);
%             angle1 = acos(dot1);
%             angle2 = acos(dot2);
%             cotans1 = cot(angle1);
%             cotans2 = cot(angle2);
            
            cotans = dot(ki, kj, 2) ./ sqrt(sum(cross(ki, kj, 2).^2,2));
            w = sum(cotans);
            sum_w = sum_w + w;
            
            sp_L = [sp_L; [i j w]];
            
            lenji = norm(vi - vj)^2;
            accum_voronoi = accum_voronoi + w * lenji;
        end
        sp_L = [sp_L; [i i -sum_w]];
        %accum_voronoi = accum_voronoi / sum_w;
        
        %accum_voronoi = accum_voronoi / 8;
        sp_W(i) = 1/(2 * accum_voronoi);
        %sp_W(i) = 1/(2 * accum_voronoi);
        
        L = sparse(sp_L(:, 1), sp_L(:, 2), sp_L(:, 3));
        W = spdiags(sp_W, 0, num_vertices, num_vertices);
        %L(i,:) = accum_vectors / (2 * accum_voronoi);
    end
end

function [L] = calc_uniform_laplacian_vector(V, F, FxV, VxV) 
    num_vertices = size(V, 1);
    L = zeros(num_vertices, 3);
    
    for i=1:num_vertices %iterate over all vertices      
        [adj_vertices, ~] = find(VxV(:, i));
        num_adjacent_vertices = numel(adj_vertices);
        
        accum = sum(V(adj_vertices, :)) - num_adjacent_vertices * V(i,:); %avoid loop at the cost of readability ;)

%         for f=1:num_adjacent_vertices
%             j = adj_vertices(f);
%             %v_ij = V(j,:) - V(i,:); 
%             accum = accum + V(j,:);  
%         end
        L(i,:) = accum;% / num_adjacent_vertices;
    end
end

function [L, K] = calc_cotan_laplacian_vector(V, F, FxV, VxV) 
    num_vertices = size(V, 1);
    
    L = zeros(num_vertices, 3);
    K = zeros(num_vertices, 1);
    
    for i=1:num_vertices %iterate over all vertices  
        vi = V(i, :);
        [adj_vertices, ~] = find(VxV(:, i));
        num_adj_vertices = size(adj_vertices, 1);
         
        %[i, s, j] = find(FxV(:, i))
        
        [adj_faces, ~] = find(FxV(:, i));
        adj_triangles = F(adj_faces, :);
        num_adj_triangles = size(adj_triangles, 1);
        
        accum_vectors = [ 0 0 0 ];
        accum_area = 0;
        accum_voronoi = 0;
        accum_angle = 0;
        
        for f=1:num_adj_triangles
            cur_triangle = adj_triangles(f, :);
            cur_triangle = cur_triangle(cur_triangle~=i); %remove center point
            o_verts = V(cur_triangle, :);
            vec1 = o_verts(1, :) - vi;
            vec2 = o_verts(2, :) - vi;
            
            %angle in vi
            angle_vi = acos(dot(vec1 / norm(vec1), vec2 / norm(vec2)));
            %triangle area
            area = norm(cross(vec1, vec2))/2;
            
            vec1_o1 = -vec1;
            vec2_o1 = o_verts(2, :) - o_verts(1, :);
            angle_o1 = acos(dot(vec1_o1 / norm(vec1_o1), vec2_o1 / norm(vec2_o1), 2));
            
            vec1_o2 = -vec2;
            vec2_o2 = -vec2_o1;
            angle_o2 = acos(dot(vec1_o2 / norm(vec1_o2), vec2_o2 / norm(vec2_o2), 2));
            
            %voronoi area
            voronoi_area = (1/8)*( norm(vec1)^2/tan(angle_o2) + norm(vec2)^2/tan(angle_o1) ); 
            
            accum_voronoi = accum_voronoi + voronoi_area;
            accum_angle = accum_angle + angle_vi;
            accum_area = accum_area + area / 3;
        end
        
        K(i) = (2*pi - accum_angle) / accum_voronoi;
        
        %step through neighborhood of v_i
        for f=1:num_adj_vertices
            j = adj_vertices(f);
            vj = V(j, :);
            ij = vj - vi; 
            
            [rows, ~] = find(adj_triangles == j);
            cur_triangles = adj_triangles(rows, :);
            ind_ks = cur_triangles(cur_triangles~=j & cur_triangles~=i);
            vks = V(ind_ks, :);
            ki = vi - vks;
            kj = vj - vks;
            
%             angles = acos(dot(sqrt(sum(abs(ki).^2,2)), sqrt(sum(abs(kj).^2,2)),2));
%             cotans = cot(angles);

%             cotan1 = dot(ki(1,:), kj(1,:)) / norm(cross(ki(1,:), kj(1,:)))
%             cotan2 = dot(ki(2,:), kj(2,:)) / norm(cross(ki(2,:), kj(2,:)));
            
            cotans = dot(ki, kj, 2) ./ sqrt(sum(cross(ki, kj, 2).^2,2));
            w = sum(cotans);
            
            accum_vectors = accum_vectors + w * ij;  
        end
        
        L(i,:) = accum_vectors / (2 * accum_area);
    end
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