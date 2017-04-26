% @file     calc_adjacent_faces.m
% @author   anna fruehstueck
% @date     18/02/2017
%
% calculate the indices of all adjacent faces for each face

function [FxF] = calc_adjacent_faces(F)
% UNFORTUNATELY SUPER INEFFICIENT
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

%%%%%%%%%% CALCULATION PROGRESS OUTPUT %%%%%%%%%%%%%%%%%%
lPrompt = 10; %length of command window prompt: this may vary depending on MATLAB version
fprintf('\n\n\n');
str = sprintf('generating face adjacency matrix... %d%%', 0);
fprintf(str);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;
num_faces = size(F, 1);
indices = [];
for k=1:num_faces-1 %iterate over all faces
    common_vertices = ismember(F(k:end, :), F(k, :));
    rows = sum(common_vertices, 2) == 2; % check if at least 2 elements are common
    rows(1) = false; %exclude self-similarity
    if any(rows)
        similar = (k-1) + find(rows); %row indices of faces similar to current face
        current = repmat(k, length(similar), 1);
        indices = [indices; [current similar]];
    end
    
    %calculation progress output
    if mod(k, 100) == 0
        lStr = length(str);
        str = sprintf('generating face adjacency matrix... %d%%', round(100*k/num_faces));
        [char(8)*ones(1, lStr + lPrompt), str]
    end
end

indices = [indices; fliplr(indices)]; %preserve symmetry
entries = length(indices)
toc;

FxF = sparse(indices(1:entries,1), indices(1:entries,2), ones(entries, 1), length(F), length(F));
end
