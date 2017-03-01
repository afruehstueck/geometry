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

function [L, W] = calc_cotan_laplacian3(V, F, FxV, VxV) 
    num_vertices = size(VxV, 1);

    sp_L = [];
    sp_W = zeros(num_vertices, 1);
    
    for i=1:num_vertices %iterate over all vertices  
        vi = V(i, :);
        [adj_verts, ~] = find(VxV(:, i));
        num_adj_verts = size(adj_verts, 1);
         
        [adj_faces, ~] = find(FxV(:, i));
        adj_triangles = F(adj_faces, :)';
        num_adj_triangles = size(adj_triangles, 2);
        %remove all occurrences of current index from adjacent triangles
        adj_triangles = reshape(adj_triangles(adj_triangles~=i),[2, num_adj_triangles])';
        
        accum_area = 0;
        %step through neighborhood of v_i
        sum_w = 0;        
        for n = 1:num_adj_triangles
            tri = adj_triangles(n, :);
            j = tri(1);
            k = tri(2);
            vj = V(j, :);
            vk = V(k, :);
            
            ij = vj - vi;
            ik = vk - vi;
            kj = vj - vk;
            l_ij = norm(ij);
            l_ik = norm(ik);
            l_kj = norm(kj);
            
            angle_i = acos(dot(ij / l_ij, ik / l_ik));
            angle_j = acos(dot(-ij / l_ij, -kj / l_kj));
            angle_k = acos(dot(-ik / l_ik, kj / l_kj));
            
            cotan_j = cot(angle_j);
            cotan_k = cot(angle_k);
            %cotans = dot(ki, kj, 2) ./ sqrt(sum(cross(ki, kj, 2).^2,2));
            
            sp_L = [sp_L; [i j 0.5*cotan_j]];
            sp_L = [sp_L; [i k 0.5*cotan_k]];
            sum_w = sum_w + 0.5*(cotan_j + cotan_k);
                    
            tri_area = norm(cross(ij, ik)) / 2;
            if is_obtuse([vi; vj; vk]) %check left triangle
                if angle_i > pi/2
                    area = tri_area / 2;
                else
                    area = tri_area / 4;
                end
            else
                area = (l_ik^2*cotan_j + l_ij^2*cotan_k) / 8;
            end
            
            accum_area = accum_area + area;     
        end
        
        sp_L = [sp_L; [i i -sum_w]];
        sp_W(i) = 1/ (2 * accum_area);
    end
                        
    L = sparse(sp_L(:, 1), sp_L(:, 2), sp_L(:, 3));
    W = spdiags(sp_W, 0, num_vertices, num_vertices);
end

function [L, W] = calc_cotan_laplacian2(V, F, FxV, VxV) 
    num_vertices = size(VxV, 1);

    sp_L = [];
    sp_W = zeros(num_vertices, 1);
    
    for i=1:num_vertices %iterate over all vertices  
        vi = V(i, :);
        [adj_verts, ~] = find(VxV(:, i));
        num_adj_verts = size(adj_verts, 1);
         
        [adj_faces, ~] = find(FxV(:, i));
        adj_triangles = F(adj_faces, :)';
        num_adj_triangles = size(adj_triangles, 2);
        %remove all occurrences of current index from adjacent triangles
        adj_triangles = reshape(adj_triangles(adj_triangles~=i),[2, num_adj_triangles])';
        
        accum_area = 0;
        accum_voronoi = 0;
        accum_angle = 0;
        %step through neighborhood of v_i
        sum_w = 0;
        sum_w_lens = 0;    
        
        ct = num_adj_triangles;
        
        tri_left = adj_triangles(1, :);  
        j = tri_left(1);
        k_left = tri_left(2);   
        while ct > 0
            vj = V(j, :);
            
            [row, ~] = find(adj_triangles == j);
            tri_right = adj_triangles(row, :);
            k_right = tri_right(tri_right ~= j & tri_right ~= k_left);
            
            vks = V([k_left; k_right], :);
            ki = vi - vks;
            kj = vj - vks;

            cotans = dot(ki, kj, 2) ./ sqrt(sum(cross(ki, kj, 2).^2,2));

            w = sum(cotans); % 
            w = max(0, w);
            
            sp_L = [sp_L; [i j w]];
            sum_w = sum_w + w;
                    
            ij = vj - vi;
            ik = vks(1, :) - vi;
            kj = vj - vks(1, :);
            l_ij = norm(ij);
            l_ik = norm(ik);
            l_kj = norm(kj);
            tri_area = norm(cross(ij, ik)) / 2;
            if is_obtuse([vi; vj; vks(1, :)]) %check left triangle
                angle_i = acos(dot(ij / l_ij, ik / l_ik));
                if angle_i > pi/2
                    area = tri_area / 2;
                else
                    area = tri_area / 4;
                end
            else
                angle_j = acos(dot(-ij / l_ij, -kj / l_kj));
                angle_k = acos(dot(-ik / l_ik, kj / l_kj));

                area = (l_ik^2/tan(angle_j) + l_ij^2/tan(angle_k)) / 8;
            end
            
            accum_area = accum_area + area;
            %accum_area = accum_area + area;
            
            %move to right triangle
            k_left = j;         
            j = k_right;
            ct = ct-1;
        end
        %accum_area = accum_area / sum_w;
        
        sp_L = [sp_L; [i i -sum_w]];
        sp_W(i) = 1/ (2 * accum_area);
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
              
        %step through neighborhood of v_i
        sum_w = 0;
        sum_A = 0;
      
        for f=1:num_adj_verts
            j = adj_verts(f);
            vj = V(j, :); 
            [rows, ~] = find(adj_triangles == j);
            cur_triangles = adj_triangles(rows, :);
            ind_ks = cur_triangles(cur_triangles~=j & cur_triangles~=i);
            vks = V(ind_ks, :);
            ij = vj - vi;
            ki = vi - vks;
            kj = vj - vks;
            
            %if is_obtuse([vi, vj, 
            
            cotans = dot(ki, kj, 2) ./ sqrt(sum(cross(ki, kj, 2).^2,2));
                
            tri_areas = sqrt(sum((cross(repmat(ij, size(vks,1), 1), vks, 2)).^2,2)) / 2;
            w = sum(cotans); % 0.5 *   
            w = max(0, w);
%             if w < 0
%                 A = sum(tri_areas) / 6;
%             else    
%                 A = (1/8) * w * norm(vi - vj)^2;
%             end
            %w = max(0, w);
%             if w < 0 %angles obtuse
%                 w = 0;%accum_voronoi + w;% * lenji;
%                 %w = max(0, w); %exclude negative weights
%             end            
            sum_w = sum_w + w;
            sum_A = sum_A + sum(tri_areas) / 6;
                       
            sp_L = [sp_L; [i j w]];
        end
        sp_L = [sp_L; [i i -sum_w]];
        %disp(['1/2sumw = ', num2str(1/(2*sum_w))]);
        %disp(['sumA/sumw = ', num2str(sum_A/sum_w)]);
        %sp_W(i) = 1/(2 * sum_A);
        sp_W(i) = sum_A/sum_w;
    end
    L = sparse(sp_L(:, 1), sp_L(:, 2), sp_L(:, 3));
    W = spdiags(sp_W, 0, num_vertices, num_vertices);
end


       
%         for f=1:num_adj_triangles
%             cur_triangle = adj_triangles(f, :);
%             cur_triangle = cur_triangle(cur_triangle~=i); %remove center point
%             verts = V(cur_triangle, :);
%             A = vi; %current vertex
%             B = verts(1, :);
%             C = verts(2, :);
%             AB = B - A;
%             AC = C - A;
%             BC = C - B;
%             l_AB = norm(AB);
%             l_AC = norm(AC);
%             l_BC = norm(BC);
%             %triangle area
%             tri_area = norm(cross(AB, AC)) / 2;
%             
%             %unit vectors
%             uAB = AB / l_AB;
%             uAC = AC / l_AC;
%             uBC = BC / l_BC; 
%             
%             %angles in A, B, C
%             angle_A = acos(dot(uAB, uAC));
%             angle_B = acos(dot(-uAB, uBC));
%             angle_C = acos(dot(-uAC, -uBC));
%             
%             %voronoi area
%             if is_obtuse([A; B; C])
%                 if angle_A > pi/2
%                     %disp('a obt')
%                     voronoi_area = tri_area / 2;
%                 else
%                     %disp('obt wo a')
%                     voronoi_area = tri_area / 4;
%                 end
%             else
%                 voronoi_area = (1/8)*( l_AB^2 / tan(angle_C) + l_AC^2 / tan(angle_B) ); 
%             end
%             
%             accum_voronoi = accum_voronoi + voronoi_area;
%             accum_angle = accum_angle + angle_A;
%             accum_area = accum_area + (tri_area / 3);
%         end

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