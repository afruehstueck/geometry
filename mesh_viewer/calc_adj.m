% @file     calc_adj.m
% @author   anna fruehstueck
% @date     18/02/2017

function [S] = calc_adj(F)
% tic;
% nrows = size(F, 1);
% pairs = [];
% for k=1:nrows
%     for j=k+1:nrows
%         if sum(ismember(F(k,:), F(j,:))) == 2
%             pairs = [pairs; [k, j]; [j, k]];
%         end
%     end  
% end
% %pairs = [pairs; [(1:nrows)' (1:nrows)']]; %add self-symmetry
% entries = size(pairs, 1)
% MAT = sparse(pairs(1:entries, 1), pairs(1:entries, 2), ones(entries, 1), entries, entries);
% toc;

%F = F(1:30, :)

lPrompt = 10; %length of command window prompt: this may vary depending on MATLAB version

fprintf('\n\n\n');

str = sprintf('generating adjacency matrix: %d%%', 0);
fprintf(str);

tic;
num_rows = size(F, 1);
ind = [];
for k=1:num_rows-1 %iterate over all rows
    common_vertices = ismember(F(k:end, :), F(k, :));
    rows = sum(common_vertices, 2) == 2; % check if at least 2 elements are common
    rows(1) = false; %exclude self-similarity
    if any(rows)
        similar = (k-1) + find(rows); %row indices of faces similar to current face
        current = repmat(k, length(similar), 1);
        ind = [ind; [current similar]];
    end
    
    if mod(k, 100) == 0
        lStr = length(str);
        str = sprintf('generating adjacency matrix: %d%%', round(100*k/num_rows));
        [char(8)*ones(1, lStr + lPrompt), str]
    end
end
ind = [ind; fliplr(ind)];
entries = length(ind)
toc;

S = sparse(ind(1:entries,1), ind(1:entries,2), ones(entries, 1), entries, entries);
spy(S)
end
