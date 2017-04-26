% from http://www.mathworks.com/matlabcentral/fileexchange/34428-voronoilimit
function [V, C] = voronoiBounded(x, y, bounds)
    % add 4 additional edges
    xA = [x; 0;  0; -2; 2];
    yA = [y; -2; 2;  0; 0];

    [vi, ci] = voronoin([xA, yA]);

    % remove the last 4 cells
    C = ci(1:end-4);
    V = vi;
    
    for idx=1:length(C)
            % convert the contour coordinate to clockwise order
            pts = V(C{idx},:);
            K = convhull(pts);
            K = K(end-1:-1:1);
            C{idx} = C{idx}(K);
            X2 = pts(K, 1);
            Y2 = pts(K, 2);
            % polybool restricts polygons to domain
            % if all points are inside the bounding box, then skip it
            if (all((X2 <= max(bounds(:, 1))) & (X2 >= min(bounds(:, 1))) & (Y2 <= max(bounds(:, 2))) & (Y2 >= min(bounds(:, 2))))) continue; end;
        
            [vert_x, vert_y] = polybool('intersection', bounds(:, 1), bounds(:, 2), X2, Y2);
            
            inds = nan(1, length(vert_x));
            for vert_i = 1:length(vert_x)
                if any(V(:, 1) == vert_x(vert_i)) && any(V(:, 2) == vert_y(vert_i))
                    i_x = find(V(:, 1) == vert_x(vert_i));
                    i_y = find(V(:, 2) == vert_y(vert_i));
                    for ib=1:length(i_x)
                        if any(i_x(ib) == i_y)
                            inds(vert_i) = i_x(ib);
                        end
                    end
                    if isnan(inds(vert_i)) == 1
                        lv = length(V);
                        V(lv+1,1) = vert_x(vert_i);
                        V(lv+1,2) = vert_y(vert_i);
                        inds(vert_i) = lv + 1;
                    end
                else
                    lv = length(V);
                    V(lv+1, 1) = vert_x(vert_i);
                    V(lv+1, 2) = vert_y(vert_i);
                    inds(vert_i) = lv + 1;
                end
            end
            C{idx} = inds;
    end
end