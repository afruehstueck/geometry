% @file     read_obj.m
% @author   anna fruehstueck
% @origin   alec jacobson
%           modified from source found at
% @source   http://www.alecjacobson.com/weblog/?p=917
%
% Reads a .obj mesh file and outputs the vertex and face list
% assumes a 3D triangle mesh and ignores everything but:
% v x y z and f i j k [l] lines
% Input:
%  filename string of obj file's path
%
% Output:
%  V  number of vertices x 3 array of vertex positions
%  F  number of faces x 3 array of face indices
  
function [V,F] = read_obj(filename)
  V = zeros(0, 3);
  F = zeros(0, 3);
  vertex_index = 1;
  face_index = 1;
 
  fid = fopen(filename, 'rt');
  line = fgets(fid);
  fn = filename;

  %progress output
  vct = 0;
  fct = 0;
  showProgress = true;
  reverseStr = '';
    
  while ischar(line)
    [token, remain] = strtok(line);
    switch token %check first letter
        case 'v' %vertex
            values = sscanf(remain, '%f');
            if ~isempty(values)
              V(vertex_index, :) = values;
              vertex_index = vertex_index + 1;
            end
        case 'f' %face
            values = sscanf(remain, '%d');
            numel = size(values, 1);
            if numel == 3 %triangle
                F(face_index, :) = values;
                face_index = face_index + 1;
            elseif numel == 4 %subdivide quads into two triangles
                F(face_index, :) = values(1:3);
                F(face_index + 1, :) = [values(3:4); values(1)];
                face_index = face_index + 2;
            else
                values6a = sscanf(remain, '%d/%d %d/%d %d/%d');
                values6b = sscanf(remain, '%d//%d %d//%d %d//%d');
                values9 = sscanf(remain, '%d/%d/%d %d/%d/%d %d/%d/%d');
                values8 = sscanf(remain, '%d/%d %d/%d %d/%d %d/%d');
                if size(values6a, 1) == 6 %face with normal indices
                    % remove normal
                    F(face_index,:) = values6a(1:2:end);
                    face_index = face_index + 1;
                elseif size(values6b, 1) == 6 %face with normal and texture indices
                    % remove normal
                    F(face_index,:) = values6b(1:2:end);
                    face_index = face_index + 1;
                elseif size(values9, 1) == 9 %face with normal and texture indices
                    % remove normal and texture indices
                    F(face_index,:) = values9(1:3:end);
                    face_index = face_index + 1;
                elseif size(values8, 1) == 8 %face with normal indices
                    % remove normal and texture indices
                    vals = values8(1:2:end);
                    F(face_index, :) = vals(1:3);
                    F(face_index + 1, :) = [vals(3:4); vals(1)];
                    face_index = face_index + 2;
                end
            end
    end
    
    %output progress of loading in console
    if mod(vertex_index, 100) == 0
        vct = vertex_index;
        showProgress = true;
    end
    if mod(face_index, 100) == 0
        fct = face_index;
        showProgress = true;
    end
    
    if showProgress
        msg = sprintf('loading %s\n%d vertices\n%d faces\n', fn, vct, fct);
        fprintf([reverseStr, msg]);
        %clever trick: \b is the backspace character
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
        showProgress = false;
    end
    
    line = fgets(fid);
  end
  fclose(fid);
end