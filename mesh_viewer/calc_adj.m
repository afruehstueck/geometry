clc;
clear;
[V, F] = read_obj('../data/shuttle.obj'); 

lPrompt = 10; %length of command window prompt: this may vary depending on MATLAB version
fA = 1;
% str = sprintf('iterating %d of %d faces', fA, size(F, 1));
% fprintf(str);
% %     if mod(fA, 10) == 0
% %         str = sprintf('iterating %d of %d faces', fA, size(F, 1));
% %         lStr = length(str);
% %         [char(8)*ones(1, lStr + lPrompt), str]
% %     end

tic;
nrows = size(F, 1);
pairs = [];
for k=1:nrows
    for j=k+1:nrows
        if sum(ismember(F(k,:), F(j,:))) == 2
            pairs = [pairs; [k, j]; [j, k]];
        end
    end  
end
pairs = [pairs; [(1:nrows)' (1:nrows)']]; %add self-symmetry
entries = size(pairs, 1)
MAT = sparse(pairs(1:entries, 1), pairs(1:entries, 2), ones(entries, 1), entries, entries);
toc;

%https://www.mathworks.com/matlabcentral/answers/1853-search-nx3-array-row-wise-for-2-equal-elements
%Wayyyyyyyy faster - but there is a bug somewhere
tic;
[f_rows, f_cols] = size(F);
ind = [];% [(1:m)' (1:m)']; %put all self-similarities
for k=1:1%f_rows-1
    %u = F(k, :);%unique(F(k, :)); % representation of kth row
    [in, J] = ismember(F(k:end, :), F(k, :));
    in
    J
    l = f_rows - k + 1
    r = repmat((1:l).', f_cols, 1)
    c = accumarray([r(in) J(in)], 1, [l f_cols]) % count
    c = bsxfun(@min, c, c(1,:)) % clip
    rows = sum(c, 2)==2 % check if at least 2 elements are common
    rows(1) = false;
    if any(rows)
        rows(1) = true;
        similar = (k-1) + find(rows);
        
        combos = nchoosek(similar, 2);
        ind = [ind; combos; fliplr(combos)];
    end
end
% sortedres = sortrows(Res);
% sortedind = sortrows(ind);
%compare = [[sortedres; zeros(max(length(ind)-length(Res),0),2)] [sortedind; zeros(max(length(Res)-length(ind), 0),2)]]

% spent = length(ind)
% B = sparse(ind(1:spent,1), ind(1:spent,2), ones(spent,1), spent, spent);
% toc;