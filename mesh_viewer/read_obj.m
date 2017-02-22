% @file     read_obj.m
% @author   anna fruehstueck
% @origin   alec jacobson
%           modified from source found at
% @source   http://www.alecjacobson.com/weblog/?p=917

function [V,F] = read_obj(filename)
  % Reads a .obj mesh file and outputs the vertex and face list
  % assumes a 3D triangle mesh and ignores everything but:
  % v x y z and f i j k [l] lines
  % Input:
  %  filename string of obj file's path
  %
  % Output:
  %  V  number of vertices x 3 array of vertex positions
  %  F  number of faces x 3 array of face indices
  %
  V = zeros(0, 3);
  F = zeros(0, 3);
  vertex_index = 1;
  face_index = 1;
 
  fid = fopen(filename, 'rt');
  line = fgets(fid);
  
  %progress output
  %fn = strsplit(filename,'/');
  %fn = fn{end};
  fn = filename;
  lPrompt = 10; %length of command window prompt: this may vary depending on MATLAB version
  vct = 0;
  fct = 0;
  str = sprintf('loading %s\n%d vertices\n%d faces', fn, vct, fct);
  fprintf(str);
  
  while ischar(line)
    [token, remain] = strtok(line);
    switch token %check first letter
        case 'v' %vertex
            values = sscanf(remain, '%f');
            if ~isempty(values)
              V(vertex_index, :) = values;
              vertex_index = vertex_index + 1;
            end
            
            %progress output
            if mod(vertex_index, 100) == 0
                vct = vertex_index;
                str = sprintf('loading %s\n%d vertices\n%d faces',fn, vct, fct);
                lStr = length(str);
                [char(8)*ones(1, lStr + lPrompt), str]
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
                values6 = sscanf(remain, '%d//%d %d//%d %d//%d');
                values9 = sscanf(remain, '%d/%d/%d %d/%d/%d %d/%d/%d');
                if size(values6, 1) == 6 %face with normal indices
                    % remove normal
                    F(face_index,:) = values6(1:2:end);
                    face_index = face_index + 1;
                elseif size(values9, 1) == 9 %face with normal and texture indices
                    % remove normal and texture indices
                    F(face_index,:) = values9(1:3:end);
                    face_index = face_index + 1;
                end
            end
            
            %progress output
            if mod(face_index, 100) == 0
                fct = face_index;
                str = sprintf('loading %s\n%d vertices\n%d faces', fn, vct, fct);
                lStr = length(str);
                [char(8)*ones(1, lStr + lPrompt), str]
            end
        otherwise
%            if ~isempty(line)
%             fprintf('Ignored: %s', line); 
%            end
    end
    line = fgets(fid);
  end
%   elems = all(F);
%   if elems(4) == 0
%       F(:,4)=[]
%   end
  fclose(fid);
end