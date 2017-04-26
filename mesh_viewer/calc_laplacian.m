% @file     calc_laplacian.m
% @author   anna fruehstueck
% @date     21/02/2017
%
% calculate the laplacian matrix for data with vertices V and faces F
% 'type' specifies whether the uniform or cotangent laplace-beltrami is
% calculated
% there are several approaches for calculating the cotangent laplacian in
% the different functions I tried out below.
% the most recent working function is calc_cotan_laplacian_per_face

function [M, D, K] = calc_laplacian(type, V, F) 
    K = 0;
    D = 0;
    if strcmp(type, 'uniform')
        [M, D] = calc_uniform_laplacian(V, F);
    elseif strcmp(type, 'cotan')
        [M, D, K] = calc_cotan_laplacian_per_face(V, F);
    end
end

function [M, D] = calc_uniform_laplacian(V, F) 
    num_vertices = size(V, 1);
    sp_M = [];
    sp_W = zeros(num_vertices, 1);
    neighbors = findNeighbors(V, F);
    
    for i=1:num_vertices %iterate over all vertices      
        adj_verts = neighbors{i}';
        num_adj_verts = numel(adj_verts);
      
        sp_W(i) = 1/num_adj_verts;
        sp_M = [sp_M; [i i -num_adj_verts;
                       repmat(i, num_adj_verts, 1) adj_verts ones(num_adj_verts, 1)]];
    end
    
    M = sparse(sp_M(:, 1), sp_M(:, 2), sp_M(:, 3));
    D = spdiags(sp_W, 0, num_vertices, num_vertices);
end

%try to step through incident triangles and calculate voronoi area
%triangle-wise
function [L, W] = calc_cotan_laplacian_per_triangle(V, F, FxV, VxV) 
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
        cur_L = zeros(2*num_adj_triangles, 3);
        for n = 1:num_adj_triangles
            tri = adj_triangles(n, :);
            j = tri(1);
            k = tri(2);
            vj = V(j, :);
            vk = V(k, :);
            
            ij = vj - vi; ji = -ij;
            ik = vk - vi; ki = -ik;
            kj = vj - vk; jk = -kj;
            
            %angle_j = atan2(norm(cross(ji, jk)), dot(ji, jk)); %numerically more stable than calculating the acos
            %angle_k = atan2(norm(cross(ki, kj)), dot(ki, kj));
            %cotan_j = cot(angle_j);
            %cotan_k = cot(angle_k);
            
            cotan_j = dot(ji, jk) ./ norm(cross(ji, jk));
            cotan_k = dot(ki, kj) ./ norm(cross(ki, kj));
            %cotan_j = max(cotan_j, 0);
            %cotan_k = max(cotan_k, 0);

            cur_L((2*n)-1:2*n, :) = [i j cotan_k; i k cotan_j];
            sum_w = sum_w + cotan_j + cotan_k;
               
            area = norm(cross(ij, ik)) / 2;
%             if is_obtuse([vi; vj; vk]) %check left triangle
%                 %angle_i = acos(dot(-(ji / l_ji), -(ki / l_ki)));
%                 nc_ij_ik = norm(cross(ij, ik));
%                 angle_i = atan2(nc_ij_ik, dot(ij, ik));
%                 tri_area = nc_ij_ik / 2;
%                 if angle_i > pi/2
%                     area = tri_area / 2;
%                 else
%                     area = tri_area / 4;
%                 end
%             else
%                 area = (norm(ji)^2*cotan_j + norm(ki)^2*cotan_k) / 8;
%             end
            accum_area = accum_area + (area/3);     
        end
        sp_L = [sp_L; cur_L];
        sp_L = [sp_L; [i i -sum_w]];
%         if(accum_area < 3)
%             sp_W(i) = 0;
%         else
            sp_W(i) = 1 / (2 * accum_area); %gets completely wonky when area gets small
%         end
        
    end
    
%     disp(['cotan L: [', num2str(min(sp_L(:, 3))), ',', num2str(max(sp_L(:, 3))), ']']);   
%     disp(['cotan W: [', num2str(min(sp_W)), ',', num2str(max(sp_W)), ']']);
%     sorted_L = sort(sp_L(:, 3), 'descend');
%     sorted_W = sort(sp_W, 'descend');
%     disp(['cotan uniform L: ', num2str(sorted_L(1:12,:)'), ' ... ', num2str(sorted_L(end-12:end,:)')]);   
%     disp(['cotan uniform W: ', num2str(sorted_W(1:12,:)'), ' ... ', num2str(sorted_W(end-12:end,:)')]);   
%     
    L = sparse(sp_L(:, 1), sp_L(:, 2), sp_L(:, 3));
    W = spdiags(sp_W, 0, num_vertices, num_vertices);
end

%trying to step through neighboring triangles two at a time 
function [M, D, K] = calc_cotan_laplacian_per_triangle_pair(V, F, FxV, VxV) 
    num_vertices = size(V, 1);

    sp_M = [];
    accum_areas = zeros(num_vertices, 1);
    accum_ws = zeros(num_vertices, 1);
    K = zeros(num_vertices, 1);
    
    neighbors = findNeighbors(V, F);
    
    for i=1:num_vertices %iterate over all vertices  
        vi = V(i, :);
        adj_verts = neighbors{i}';
         
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

            cotans = dot(ki, kj, 2) ./ norm_per_row(cross(ki, kj, 2));
            w = sum(cotans); % 
            
            sp_M = [sp_M; [i j w]];
            sum_w = sum_w + w;
                    
            ij = vj - vi;
            ik = vks(1, :) - vi;
            kj = vj - vks(1, :);
            l_ij = norm(ij);
            l_ik = norm(ik);
            l_kj = norm(kj);
            u_ij = ij/l_ij;
            u_ik = ik/l_ik;
            u_kj = kj/l_kj;
            
            angle_i = acos(dot(u_ij, u_ik));
            K(i) = K(i) + angle_i;
            tri_area = norm(cross(ij, ik)) / 2;
            if is_obtuse([vi; vj; vks(1, :)]) %check left triangle
                if angle_i > pi/2
                    area = tri_area / 2;
                else
                    area = tri_area / 4;
                end
            else
                angle_j = acos(dot(-u_ij, -u_kj));
                angle_k = acos(dot(-u_ik, u_kj));

                area = (l_ik^2/tan(angle_j) + l_ij^2/tan(angle_k)) / 8;
            end
            accum_area = accum_area + area;
            
            %move to right triangle
            k_left = j;         
            j = k_right;
            ct = ct-1;
        end        
        sp_M = [sp_M; [i i -sum_w]];
        accum_areas(i) = accum_area;
        accum_ws(i) = sum_w;
        %sp_W(i) = 1 / (2 * accum_area);
    end
                    
    M = sparse(sp_M(:, 1), sp_M(:, 2), sp_M(:, 3));
    K = (2*pi - K) ./ accum_area;
    sp_W = (2 * accum_areas).^(-1);
    %sp_W = (2 * accum_ws).^(-1);
    D = spdiags(sp_W, 0, num_vertices, num_vertices);
end

function [M, D] = calc_cotan_laplacian(V, F, FxV, VxV) 
    num_vertices = size(VxV, 1);
    sp_M = [];%zeros(num_vertices*num_vertices, 3);
    sp_D = zeros(num_vertices, 1);
   
%     min_adj_faces = 10000;
%     max_adj_faces = 0;
    for i=1:num_vertices %iterate over all vertices  
        v_i = V(i, :);
        [adj_verts, ~] = find(VxV(:, i));
        num_adj_verts = size(adj_verts, 1);
         
        [adj_faces, ~] = find(FxV(:, i));
        adj_triangles = F(adj_faces, :);
        num_adj_faces = size(adj_faces, 1);
        
%         min_adj_faces = min(num_adj_faces, min_adj_faces);
%         max_adj_faces = max(num_adj_faces, max_adj_faces);
              
        %step through neighborhood of v_i
        sum_A = 0;
      
        num_new_entries = num_adj_verts + 1;
        cur_M = zeros(num_new_entries, 3);
        
        for f=1:num_adj_verts
            j = adj_verts(f);
            v_j = V(j, :); 
            [rows, ~] = find(adj_triangles == j);
            cur_triangles = adj_triangles(rows, :);
            ind_ks = cur_triangles(cur_triangles~=j & cur_triangles~=i);
            v_ks = V(ind_ks, :);
            vec_ij = v_j - v_i;
            vec_ki = v_i - v_ks;
            vec_kj = v_j - v_ks;
            
            quots = norm_per_row(cross(vec_ki, vec_kj, 2));
            
            cotans = dot(vec_ki, vec_kj, 2) ./ quots;
                
            w = sum(cotans);  
            %w = max(0, w);
            %tri_areas = sqrt(sum((cross(repmat(vec_ij, size(v_ks, 1), 1), v_ks, 2)).^2, 2)) / 2;
            tri_areas = quots / 2;
            
            l_ij_sq = norm(vec_ij)^2;
            voronoi_areas = (cotans .* l_ij_sq);
            %add thirds of triangle areas
            sum_A = sum_A + sum(tri_areas) / 3;
            %sum_A = sum_A + sum(voronoi_areas);
            cur_M(f, :) = [i j w];
        end
        
        sum_w = sum(cur_M(:, 3));
        %cur_L(:, 3) = cur_L(:, 3) ./ sum_w;
        cur_M(num_new_entries, :) = [i i -sum_w];
        
        sp_M = [sp_M; cur_M];              
        
        %sum_A = sum_A / 8;
        
        %sp_W(i) = 1 / sum_A;

        %sp_W(i) = 1 / sum_A; %each triangle area was taken twice, which eliminates the 2*A
        
        sp_D(i) = 1/(2 * sum_w);
    end
    
    M = sparse(sp_M(:, 1), sp_M(:, 2), sp_M(:, 3));
    D = spdiags(sp_D, 0, num_vertices, num_vertices);
end

% step through all faces stored in face list and calculate contributions of
% each face to adjacent vertices
% M = matrix of vertex weights, 
% D = diag matrix of areas, 
% K = curvature
function [M, D, K] = calc_cotan_laplacian_per_face(V, F)  
    num_vertices = size(V, 1);
    num_faces = size(F, 1);
    
    %collect contributions for sparse matrix in sp_L: 
    %3 entries per vertex, 3 vertices per face
    sp_M_col1 = zeros(num_faces * 9, 1);
    sp_M_col2 = zeros(num_faces * 9, 1);
    sp_M_col3 = zeros(num_faces * 9, 1);
    
    K = zeros(num_vertices, 1);
    %ctK = zeros(num_vertices, 1); %counts number of faces contributing to vertex (debug!)
    
    mixed_voronoi_areas = zeros(num_vertices, 1);
    voronoi_areas = zeros(num_vertices, 1);
    simple_triangle_areas = zeros(num_vertices, 1);
   
    permutations = [1 2 3; 
                    2 3 1; 
                    3 1 2]; %all combinations of indices
                 
    Ps = permutations(:, 1); %current points P
    Qs = permutations(:, 2); %adjacent points Q
    Rs = permutations(:, 3); %adjacent points R
        
    max_angle = 0;
    min_angle = 4*pi;
    
    %calculating for all triangle vertices simultaneously:
    %current vertex is called P, adjacent vertices R and Q
    for f=1:num_faces 
        vert_inds = F(f, :)';
        verts = V(vert_inds, :);
        
        %edges pointing from P to Q
        edgesPQ = verts(Qs, :) - verts(Ps, :);
        %edges pointing from P to R
        edgesPR = verts(Rs, :) - verts(Ps, :);
        
        %calculate area for current triangle
        tri_area = norm(cross(edgesPQ(1, :), edgesPR(1, :))) / 2;

        l_edgesPR = norm_per_row(edgesPR);
        l_edgesPQ = norm_per_row(edgesPQ);
        
        %normalize edges
        edgesPR = edgesPR ./ l_edgesPR;
        edgesPQ = edgesPQ ./ l_edgesPQ;
        
        %angles in P
        %angles = acos(dot(edgesPR, edgesPQ, 2));
        %replaced dot product by faster version
        angles = acos(sum(conj(edgesPR) .* edgesPQ, 2));
        
        %cotans = dot(ki, kj, 2) ./ sqrt(sum(cross(edgesPR, edgesPQ, 2).^2,2));
        
        %angles = atan2(sqrt(sum(cross(edgesPR, edgesPQ, 2).^2, 2)), sum(conj(edgesPR) .* edgesPQ, 2));
            
        %evaluate if angles are obtuse
        obtuse_angles = angles > (pi/2);
        obtuse_triangle = any(obtuse_angles);
        
        %DEBUG
        max_angle = max(max(angles), max_angle);
        min_angle = min(min(angles), min_angle);
        
        %add up angles at vertices
        K(vert_inds) = K(vert_inds) + angles;
        %ctK(vert_inds) = ctK(vert_inds) + [1; 1; 1];
        
        %cotangent values for angles at vertices P, Q, R
        cotangents = cot(angles);
        %ARBITRARY CRASH CONDITION FOR DEBUGGING
%         if f>10 
%             K(vert_inds) = K(vert_inds, 1:2) + angles;
%         end
%          sp_L(f*9-8:f*9, :) = [vert_inds(Ps) vert_inds(Ps) -(cotangents(Rs) + cotangents(Qs)); 
%                                vert_inds(Ps) vert_inds(Qs) cotangents(Rs); 
%                                vert_inds(Ps) vert_inds(Rs) cotangents(Qs)];
        %first three: diagonal entries, then three entries for the mixed
        %contributions
        sp_M_col1(f*9-8:f*9, :) = [vert_inds(Ps); vert_inds(Ps); vert_inds(Ps)];
        sp_M_col2(f*9-8:f*9, :) = [vert_inds(Ps); vert_inds(Qs); vert_inds(Rs)];
        sp_M_col3(f*9-8:f*9, :) = [-(cotangents(Rs) + cotangents(Qs)); cotangents(Rs); cotangents(Qs)];
        
        %clamp cotangents to zero for area computation
        %cotangents = max(cotangents, [0; 0; 0]);
        
        voronoi = (l_edgesPQ.^2 .* cotangents(Rs) + l_edgesPR.^2 .* cotangents(Qs)) / 8;

        %put area three times in column vector
        tri_areas = repmat(tri_area, 3, 1);
        
        %for obtuse triangles, use half or quarter area, for non-obtuse triangles, use voronoi area
        areas_mixed = obtuse_triangle * (obtuse_angles .* (tri_areas ./ 2) + ~obtuse_angles .* (tri_areas ./ 4)) + ~obtuse_triangle * voronoi;
        
        mixed_voronoi_areas(vert_inds)   = mixed_voronoi_areas(vert_inds) + areas_mixed;
        voronoi_areas(vert_inds)         = voronoi_areas(vert_inds) + voronoi;
        simple_triangle_areas(vert_inds) = simple_triangle_areas(vert_inds) + (tri_areas ./ 3);
    end
    
    areas = mixed_voronoi_areas;
    
    %clamp small areas?????????
    areas = max(areas, ones(num_vertices, 1));
   
    %disp(['min/max total angle: ', num2str(min(K)), ', ', num2str(max(K))]); 
    K = (2*pi - K) ./ areas;
    
    sp_D = 1 ./ (areas .* 2);
   
    disp(['min/max angle: ', num2str(min_angle), ', ', num2str(max_angle)]); 
    disp(['min/max tri areas: ', num2str(min(simple_triangle_areas)), ', ', num2str(max(simple_triangle_areas))]); 
    disp(['min/max mixed voronoi areas: ', num2str(min(mixed_voronoi_areas)), ', ', num2str(max(mixed_voronoi_areas))]); 
    disp(['min/max voronoi areas: ', num2str(min(voronoi_areas)), ', ', num2str(max(voronoi_areas))]); 
    disp(['min/max curvature: ', num2str(min(K)), ', ', num2str(max(K))]); 
    disp(['min/max weights: ', num2str(min(sp_D)), ', ', num2str(max(sp_D))]); 
    
    M = sparse(sp_M_col1, sp_M_col2, sp_M_col3);
    %sparse diagonal area matrix
    D = spdiags(sp_D, 0, num_vertices, num_vertices);
end

function [L, W] = calc_cotan_laplacian4(V, F, FxV, VxV) 
    num_vertices = size(VxV, 1);

    sp_L = [];
    sp_W = zeros(num_vertices, 1);
    
    for i=1:num_vertices %iterate over all vertices  
        vi = V(i, :);
        [adj_verts, ~] = find(VxV(:, i));
        num_adj_verts = size(adj_verts, 1);
         
        [adj_faces, ~] = find(FxV(:, i));
        adj_triangles = F(adj_faces, :);
              
        %step through neighborhood of v_i
        sum_w = 0;
        sum_A = 0;
        
        cur_L = zeros(num_adj_verts + 1, 3);
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
            
            cotans = dot(ki, kj, 2) ./ norm_per_row(cross(ki, kj, 2));
                
            tri_areas = norm_per_row(cross(repmat(ij, size(vks, 1), 1), vks, 2)) / 2;
            w = sum(cotans);   
            
            sum_A = sum_A + sum(tri_areas) / 3;
                       
            cur_L(f, :) = [i j w];
        end
        
        sum_A = sum_A / 2; %every triangle area was taken twice
        %cur_L(:, 3) = cur_L(:, 3) ./ (2 * sum_A);
        
        sum_w = sum(cur_L(:, 3));
        cur_L(num_adj_verts + 1, :) = [i i -sum_w];
        
        sp_L = [sp_L; cur_L];  
        sp_W(i) = 1 /(2 * sum_A);
    end
    L = sparse(sp_L(:, 1), sp_L(:, 2), sp_L(:, 3));
    W = spdiags(sp_W, 0, num_vertices, num_vertices);
end

% evaluate obtusity of triangle
%see http://mathworld.wolfram.com/ObtuseTriangle.html
function result = is_obtuse(vertices)
    perm = [1 2 3; 2 3 1; 3 1 2]; %all combinations of indices
    edges = vertices(perm(:,3), :) - vertices(perm(:,2), :);
    sq_edge_lengths = sum(abs(edges).^2, 2); %row-wise norm squared
    compare = (sq_edge_lengths(perm(:, 2), :) + sq_edge_lengths(perm(:, 3), :) < sq_edge_lengths(perm(:, 1), :)); %evaluate all a^2 + b^2 < c^2
    result = any(compare); %is any of the inequalities true
end

% calculate row-wise norm in matrix (function not available in MATLAB)
function N = norm_per_row(M) 
    N = sqrt(sum(M.^2, 2));
end


