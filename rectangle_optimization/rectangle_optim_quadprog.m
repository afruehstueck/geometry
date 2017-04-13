% @file     rectangle_optim_quadprog.m
% @author   afruehstueck
% @date     03/04/2017

function constrained_rectangles = rectangle_optim_quadprog(rectangles, rows, cols, minimize_variables, alignment, align_in_obj_fun)    
    %% create hessian
    num_rects = size(rectangles, 1);
    % H = matrix with variables on main diagonal minimizes 
    %     distance of pointsto input points
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
    %indices at which entries for this variable start in matrix
    start_y = 1;
    start_w = 1;
    start_h = 1;
    idx_in_matrix = 1;
    if contains(minimize_variables, 'c')
        diag_squared_terms = [diag_squared_terms; repmat(rows, num_rects, 1)];
        f = [f; rectangles(:, 1)];
        idx_in_matrix = idx_in_matrix + num_rects;
    end
    
    if contains(minimize_variables, 'r')
        start_y = idx_in_matrix;
        diag_squared_terms = [diag_squared_terms; repmat(cols, num_rects, 1)];
        f = [f; rectangles(:, 2)];
        idx_in_matrix = idx_in_matrix + num_rects;
    end
    
    if contains(minimize_variables, 'w')
        start_w = idx_in_matrix;
        diag_squared_terms = [diag_squared_terms; repmat(rows*cols, num_rects, 1)];%ones(num_rects, 1)];
        f = [f; rectangles(:, 3)];
        idx_in_matrix = idx_in_matrix + num_rects;
    end
 
    if contains(minimize_variables, 'h')
        start_h = idx_in_matrix;
        diag_squared_terms = [diag_squared_terms; repmat(rows*cols, num_rects, 1)];%ones(num_rects, 1)];ones(num_rects, 1)];
        f = [f; rectangles(:, 4)];
    end
    
    if ~align_in_obj_fun
        %skip all previously calculated terms and construct diagonal of ones (only
        %squared distances are evaluated in objective function)
        diag_squared_terms = ones(numel(f), 1);
    end
    
    H = diag(diag_squared_terms);
    
    if align_in_obj_fun
        %add mixed terms
        if contains(minimize_variables, 'r')
            tri = triu((-2) * ones(cols), 1);
            %enforce row-wise minimization of ys
            for r = 1:rows
                y1 = start_y + (r-1) * cols;
                common_row_idxs = y1:y1+cols-1;
                H(common_row_idxs, common_row_idxs) = H(common_row_idxs, common_row_idxs) + tri; 
            end
        end
        
        if contains(minimize_variables, 'c')
            %enforce column-wise minimization of xs
            tri = triu((-2) * ones(rows), 1);
            for c = 1:cols
                x1 = c;
                common_col_idxs = x1:cols:cols*(rows-1)+x1;
                H(common_col_idxs, common_col_idxs) = H(common_col_idxs, common_col_idxs) + tri; 
            end
        end

        tri = triu((-2) * ones(rows*cols), 1);
        if contains(minimize_variables, 'w')
            w_inds = start_w:start_w+rows*cols-1;
            H(w_inds, w_inds) = H(w_inds, w_inds) + tri;
        end
        
        if contains(minimize_variables, 'h')
            h_inds = start_h:start_h+rows*cols-1;
            H(h_inds, h_inds) = H(h_inds, h_inds) + tri;
        end
    end
    
    % f = [-x1 -x2 ... -y1 -y2 ... -w1 -w2 ... -h1 -h2 ...]'
    f = -f;
    
    if ~align_in_obj_fun 
        %create equality constraints
        Aeq = zeros(0, numel(f));
        beq = [];
        index = 1;

        %enforce equality of all xs per column
        if contains(minimize_variables, 'c')    
            dia = diag(ones(rows-1, 1));
            widths = rectangles(:, 3);
            for c = 1:cols
                x1 = c;
                common_col_idxs = x1+cols:cols:cols*(rows-1)+x1;

                Aeq(index:index+rows-2, x1) = 1;
                Aeq(index:index+rows-2, common_col_idxs) = -dia; 
                
                %determine right side term according to alignment
                %calculate from width difference if widths are not aligned              
                if alignment(2) == 'r' && ~contains(minimize_variables, 'w')
                    beq(index:index+rows-2, 1) = (widths(common_col_idxs) - widths(x1))/2;
                elseif alignment(2) == 'l' && ~contains(minimize_variables, 'w')
                    beq(index:index+rows-2, 1) = (widths(x1) - widths(common_col_idxs))/2;
                else %alignment is centered and/or heights are aligned          
                    beq(index:index+rows-2, 1) = 0;
                end
                
                index = index+rows-1;
            end
        end
                
        %enforce equality of all ys per row  
        if contains(minimize_variables, 'r')    
            dia = diag(ones(cols-1, 1));
            heights = rectangles(:, 4);
            for r = 1:rows
                y1 = start_y + (r-1) * cols;
                h1 = 1 + (r-1) * cols;
                Aeq(index:index+cols-2, y1) = 1;
                Aeq(index:index+cols-2, y1+1:y1+cols-1) = -dia;
                
                %determine right side term according to alignment
                %calculate from height difference if heights are not aligned
                if alignment(1) == 't' && ~contains(minimize_variables, 'h')
                    beq(index:index+cols-2, 1) = (heights(h1+1:h1+cols-1) - heights(h1))/2;
                elseif alignment(1) == 'b' && ~contains(minimize_variables, 'h')
                    beq(index:index+cols-2, 1) = (heights(h1) - heights(h1+1:h1+cols-1))/2;
                else %alignment is centered and/or heights are aligned
                    beq(index:index+cols-2, 1) = 0;
                end
                
                %beq(index:index+cols-2, 1) = 0;
                index = index+cols-1;
            end
        end    
        
        %enforce equality of all ws          
        if contains(minimize_variables, 'w')    
            dia = diag(ones(rows*cols-1, 1));
            w1 = start_w;
            Aeq(index:index+rows*cols-2, w1) = 1;
            Aeq(index:index+rows*cols-2, w1+1:w1+rows*cols-1) = -dia; 
            beq(index:index+rows*cols-2, 1) = 0;
            index = index+rows*cols-1;
        end
                  
        %enforce equality of all hs
        if contains(minimize_variables, 'h')    
            dia = diag(ones(rows*cols-1, 1));
            h1 = start_h;
            Aeq(index:index+rows*cols-2, h1) = 1;
            Aeq(index:index+rows*cols-2, h1+1:h1+rows*cols-1) = -dia; 
            beq(index:index+rows*cols-2, 1) = 0;
            index = index+rows*cols-1;
        end
        
        %enforce equal spacing of rows 
        if contains(minimize_variables, 's') && rows > 2
            dia = diag(-ones(rows, 1)) + diag(2 * ones(rows-1, 1), 1) + diag(-ones(rows-2, 1), 2);
            dia(end-1:end, :) = [];
            y1 = start_y;
            common_col_idxs = y1:cols:cols*(rows-1)+y1;

            Aeq(index:index+size(dia, 1) - 1, common_col_idxs) = dia; 
            beq(index:index+size(dia, 1) - 1, 1) = 0;
            index = index+size(dia, 1);
        end
        
        %enforce equal spacing of columns
        if contains(minimize_variables, 't') && cols > 2
            dia = diag(-ones(cols, 1)) + diag(2 * ones(cols-1, 1), 1) + diag(-ones(cols-2, 1), 2);
            dia(end-1:end, :) = [];
            x1 = 1;
            
            Aeq(index:index+size(dia, 1) - 1, x1:x1+cols-1) = dia; 
            beq(index:index+size(dia, 1) - 1, 1) = 0;
            index = index+size(dia, 1);
        end
        
        % = zeros(size(Aeq, 1), 1);
        
        %minimize using equality constraints
        minim = quadprog(H, f, [], [], Aeq, beq);
    else 
        %minimize using only objective function
        minim = quadprog(H, f);
    end

    
    idxs = 1:num_rects;
    if contains(minimize_variables, 'c')
        x = reshape(minim(idxs), [num_rects, 1]);
        idxs = idxs(end)+1:idxs(end)+num_rects;
    else
        x = rectangles(:, 1);
    end
    
    if contains(minimize_variables, 'r')
        y = reshape(minim(idxs), [num_rects, 1]);
        idxs = idxs(end)+1:idxs(end)+num_rects;
    else
        y = rectangles(:, 2);
    end
    
    if contains(minimize_variables, 'w')
        w = reshape(minim(idxs), [num_rects, 1]);
        idxs = idxs(end)+1:idxs(end)+num_rects;
    else
        w = rectangles(:, 3);
    end
 
    if contains(minimize_variables, 'h')
        h = reshape(minim(idxs), [num_rects, 1]);
    else
        h = rectangles(:, 4);
    end
    
    constrained_rectangles = [x y w h];
end
