% @file     rectangle_constrain.m
% @author   afruehstueck
% @date     03/04/2017

function constrained_rectangles = rectangle_constrain(rectangles, rows, cols, align_in_obj_fun)
    objective_variables = 'xywh';
    constraint_variables = 'xywhst';
    %align_in_obj_fun = true;
    %% create hessian
    num_rects = size(rectangles, 1);
    num_vars = numel(objective_variables); %rows+cols;
    % H = matrix with 2s on main diagonal
    %     ( x1                  )
    %     |   x2                |
    %     |      y1             |
    %     |        y2           |
    % H = |           w1        |
    %     |             w2      |
    %     |                h1   |
    %     (                  h2 )
    %diag_squared_terms = ones(num_vars*num_rects, 1);
    
    %squared terms on main diagonal
    diag_squared_terms = [];
    
    %linear terms
    f = [];
    start_y = 1;
    start_w = 1;
    start_h = 1;
    idx_in_matrix = 1;
    if contains(objective_variables, 'x')
        diag_squared_terms = [diag_squared_terms; repmat(rows, num_rects, 1)];
        f = [f; rectangles(:, 1)];
        idx_in_matrix = idx_in_matrix + num_rects;
    end
    
    if contains(objective_variables, 'y')
        start_y = idx_in_matrix;
        diag_squared_terms = [diag_squared_terms; repmat(cols, num_rects, 1)];
        f = [f; rectangles(:, 2)];
        
        idx_in_matrix = idx_in_matrix + num_rects;
    end
    
    if contains(objective_variables, 'w')
        start_w = idx_in_matrix;
        diag_squared_terms = [diag_squared_terms; repmat(rows*cols, num_rects, 1)];%ones(num_rects, 1)];
        f = [f; rectangles(:, 3)];
        idx_in_matrix = idx_in_matrix + num_rects;
    end
 
    if contains(objective_variables, 'h')
        start_h = idx_in_matrix;
        diag_squared_terms = [diag_squared_terms; repmat(rows*cols, num_rects, 1)];%ones(num_rects, 1)];ones(num_rects, 1)];
        f = [f; rectangles(:, 4)];
    end
    
    if ~align_in_obj_fun
        %skip all previous terms and construct diagonal of ones (only
        %squared distances are evaluated in objective function)
        diag_squared_terms = ones(num_vars*num_rects, 1);
    end
    
    H = diag(diag_squared_terms);
    
    if align_in_obj_fun
        %add mixed terms
        if contains(objective_variables, 'y')
            tri = triu((-2) * ones(cols), 1);
            %enforce row-wise minimization of ys
            for r = 1:rows
                y1 = start_y + (r-1) * cols;
                common_row_idxs = y1:y1+cols-1;
                H(common_row_idxs, common_row_idxs) = H(common_row_idxs, common_row_idxs) + tri; 
            end
        end
        
        if contains(objective_variables, 'x')
            %enforce column-wise minimization of xs
            tri = triu((-2) * ones(rows), 1);
            for c = 1:cols
                x1 = c;
                common_col_idxs = x1:cols:cols*(rows-1)+x1;
                H(common_col_idxs, common_col_idxs) = H(common_col_idxs, common_col_idxs) + tri; 
            end
        end

        tri = triu((-2) * ones(rows*cols), 1);
        if contains(objective_variables, 'w')
            %w1 = start_w + (r-1) * cols + 1;
            w_inds = start_w:start_w+rows*cols-1;
            H(w_inds, w_inds) = H(w_inds, w_inds) + tri;
        end
        
        if contains(objective_variables, 'h')
            h_inds = start_h:start_h+rows*cols-1;
            H(h_inds, h_inds) = H(h_inds, h_inds) + tri;
        end
    end
    
    %multiply by 2 to compensate for 1/2
    H = 2 * H;
    
    % f = [-2*x1 -2*x2 ... -2*y1 -2*y2 ... -2*w1 -2*w2 ... -2*h1 -2*h2 ...]'
    f = (-2) * f;
    
    if ~align_in_obj_fun 
        %create equality constraints
        Aeq = zeros(0, num_vars*num_rects);
        index = 1;

        if contains(constraint_variables, 'x')    
            %enforce equality of all xs per column
            dia = diag((-1) * ones(rows-1, 1));
            for c = 1:cols
                x1 = c;
                common_col_idxs = x1+cols:cols:cols*(rows-1)+x1;

                Aeq(index:index+rows-2, x1) = 1;
                Aeq(index:index+rows-2, common_col_idxs) = dia; 
                index = index+rows-1;
                Aeq
            end
        end
                  
        if contains(constraint_variables, 'y')    
            dia = diag((-1) * ones(cols-1, 1));
            %eenforce equality of all ys per row
            for r = 1:rows
                y1 = start_y + (r-1) * cols;
                Aeq(index:index+cols-2, y1) = 1;
                Aeq(index:index+cols-2, y1+1:y1+cols-1) = dia; 
                index = index+cols-1;
                Aeq
            end
        end    
                  
        if contains(constraint_variables, 'w')    
            dia = diag((-1) * ones(rows*cols-1, 1));
            %eenforce equality of all ws
            w1 = start_w;
            Aeq(index:index+rows*cols-2, w1) = 1;
            Aeq(index:index+rows*cols-2, w1+1:w1+rows*cols-1) = dia; 
            index = index+rows*cols-1;
            Aeq
        end
                  
        if contains(constraint_variables, 'h')    
            dia = diag((-1) * ones(rows*cols-1, 1));
            %eenforce equality of all hs
            h1 = start_h;
            Aeq(index:index+rows*cols-2, h1) = 1;
            Aeq(index:index+rows*cols-2, h1+1:h1+rows*cols-1) = dia; 
            index = index+rows*cols-1;
            Aeq
        end
        
        %enforce equal spacing of rows 
        if contains(constraint_variables, 's') && rows > 2
            dia = diag(-ones(rows, 1)) + diag(2 * ones(rows-1, 1), 1) + diag(-ones(rows-2, 1), 2);
            dia(end-1:end, :) = [];
            y1 = start_y;
            common_col_idxs = y1:cols:cols*(rows-1)+y1;

            Aeq(index:index+size(dia, 1) - 1, common_col_idxs) = dia; 
            index = index+size(dia, 1);
            Aeq
        end
        
        %enforce equal spacing of columns
        if contains(constraint_variables, 't') && cols > 2
            dia = diag(-ones(cols, 1)) + diag(2 * ones(cols-1, 1), 1) + diag(-ones(cols-2, 1), 2);
            dia(end-1:end, :) = [];
            x1 = 1;
            
            Aeq(index:index+size(dia, 1)-1, x1:x1+cols-1) = dia; 
            index = index+size(dia, 1);
            Aeq
        end
        
        feq = zeros(size(Aeq, 1), 1);
        
        %minimize using equality constraints
        minim = quadprog(H, f, [], [], Aeq, feq);
    else 
        
        %minimize using only objective function
        minim = quadprog(H, f);
    end

    
    idxs = 1:num_rects;
    if contains(objective_variables, 'x')
        x = reshape(minim(idxs), [num_rects, 1]);
        idxs = idxs(end)+1:idxs(end)+num_rects;
    else
        x = rectangles(:, 1);
    end
    
    if contains(objective_variables, 'y')
        y = reshape(minim(idxs), [num_rects, 1]);
        idxs = idxs(end)+1:idxs(end)+num_rects;
    else
        y = rectangles(:, 2);
    end
    
    if contains(objective_variables, 'w')
        w = reshape(minim(idxs), [num_rects, 1]);
        idxs = idxs(end)+1:idxs(end)+num_rects;
    else
        w = rectangles(:, 3);
    end
 
    if contains(objective_variables, 'h')
        h = reshape(minim(idxs), [num_rects, 1]);
    else
        h = rectangles(:, 4);
    end
    
    constrained_rectangles = [x y w h]
%     x_values = reshape(rectangles(:, 1), [cols, rows])';
%     y_values = reshape(rectangles(:, 2), [cols, rows])';
%     
%     H = diag(2 * [repmat(rows, cols, 1); repmat(cols, rows, 1)]);
%     
%     f = (-2) * [sum(x_values, 1)'; sum(y_values, 2)];
%     
%     minim = quadprog(H, f);
%     
%     x = minim(1:cols);
%     y = minim(cols+1:end);
%     
%     constrained_rectangles = rectangles;
%     constrained_rectangles(:, 1) = repmat(x, rows, 1);
%     constrained_rectangles(:, 2) = repmat(y, cols, 1);
end
