% @file     matrix.m
% @author   afruehstueck
% @date     30/01/2017

%create special matrices on the basis of various input arguments
function m = matrix(varargin)
    %dimension specifies the size of the matrix
    dimension = varargin{1};
    %type is a string specifying the type of the matrix
    type      = varargin{2};
    
    varargin(1)=[]; %remove type from arguments, pass all other args on
    varargin(1)=[]; %remove dimension from arguments
    
    %check for argument errors
    if (~exist('type', 'var'))
        disp('no matrix type specified.');
        return;
    end
    if (~exist('dimension', 'var'))
        disp('no matrix dimension specified.');
        return;
    end
    if (~isfloat(dimension) || length(dimension) ~= 1 || dimension < 2 || dimension > 4 ) 
        disp('invalid matrix dimension specified.');
        return;
    end
   
    %handle different input types 
    %return matrix from respective function
    switch type
        case 'rotation'
            if dimension == 2
                m = matrix2DRotation(varargin);
            elseif dimension == 3 || dimension == 4
                m = matrix3DRotation(varargin);
            end
        case 'shear'
            if dimension == 2       
                m = matrix2DShear(varargin);
            elseif dimension == 3 || dimension == 4
                m = matrix3DShear(varargin);
            end
        case 'reflection'
            if dimension == 2
                m = matrix2DReflection(varargin);
            elseif dimension == 3 || dimension == 4
                m = matrix3DReflection(varargin);
            end
        case 'scale'
            if dimension == 2
                m = matrix2DScale(varargin);
            elseif dimension == 3 || dimension == 4
                m = matrix3DScale(varargin);
            end
        case 'translation'
            if dimension < 3 ||  dimension > 4
                disp('translation matrix only exists with homogeneous coordinates');
                m = -1;
                return;    
            end
            if dimension == 3
                m = matrix2DTranslation(varargin);
            elseif dimension == 4
                m = matrix3DTranslation(varargin);
            end
        case 'projection'
            if dimension < 3 ||  dimension > 4
                disp('translation matrix only exists with homogeneous coordinates');
                m = -1;
                return;    
            end
            m = matrix3DPerspectiveProjection(varargin);
        case 'householder'
            m = matrix3DHouseholder(varargin);
        otherwise
            disp('unknown matrix type');
    end
    
    %if homogeneous matrix is requested, pad matrix with ones
    if dimension > length(m)
        % pad 3x3 matrix with homogeneous coordinate
        m = [ [ m zeros(length(m), 1) ]; [ zeros(1, length(m)) 1 ] ];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2D ROTATION MATRIX         %
% @args(1)angle              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function m = matrix2DRotation(args) 
    if length(args) ~= 1 || length(args{1}) > 1 || ~isfloat(args{1})
        disp('2D rotation matrix expects 1 input: angle (radians)');
        m = -1;
        return;
    end
    
    a = args{1};
    m = [cos(a) -sin(a); 
         sin(a)  cos(a)];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3D ROTATION MATRIX         %
% @args(1)axis {'x','y','z'} %
%      (2)angle              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function m = matrix3DRotation(args) 
    if length(args) ~= 2 || ~ischar(args{1}) || length(args{2}) > 1 || ~isfloat(args{2})
        disp('3D rotation matrix expects 2 inputs: rotation axis {''x'', ''y'', ''z''} AND angle (radians)');
        m = -1;
        return;
    end
    
    axis = args{1};
    a    = args{2};
    
    switch axis
        case 'x'
            m = [   1     0       0   ;
                    0 cos(a) -sin(a)  ; 
                    0 sin(a)  cos(a) ];
        case 'y'
            m = [   cos(a) 0 sin(a)   ;
                        0  1     0    ; 
                   -sin(a) 0 cos(a)  ];
        case 'z'
            m = [   cos(a) -sin(a) 0  ;
                    sin(a)  cos(a) 0  ; 
                        0       0  1 ];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2D SHEAR MATRIX           %
% @args(1)axis {'x','y'}    %
%      (2)shear factor      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function m = matrix2DShear(args) 
    if (length(args) ~= 2 || ~ischar(args{1}) || length(args{2}) ~= 1 || ~isfloat(args{2}))
        disp('2D shear matrix expects 2 inputs: axis {''x'', ''y''} AND factor');
        m = -1;
        return;
    end
    
    axis = args{1};
    f    = args{2};
    switch axis
        case 'x'
            m = [ 1 f  ; 
                  0 1 ];
        case 'y'
            m = [ 1 0  ; 
                  f 1 ];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3D SHEAR MATRIX           %
% @args(1)axis {'x','y','z'}%
%      (2)shear factor 1    %
%      (3)shear factor 2    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function m = matrix3DShear(args) 
    if (length(args) ~= 3 || ~ischar(args{1}) || length(args{2}) ~= 1 || ~isfloat(args{2}) || length(args{3}) ~= 1 || ~isfloat(args{3}))
        disp('3D shear matrix expects 3 inputs: axis {''x'', ''y'', ''z''} AND factor1 AND factor2');
        m = -1;
        return;
    end
    
    axis = args{1};
    f1   = args{2};
    f2   = args{3};
    switch axis
        case 'x'
            m = [ 1 0 0  ;
                 f1 1 0  ; 
                 f2 0 1 ];
        case 'y'
            m = [ 1 f1 0  ;
                  0  1 0  ; 
                  0 f2 1 ];
        case 'z'
            m = [ 1 0 f1  ;
                  0 1 f2  ; 
                  0 0  1 ];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2D REFLECTION MATRIX          %
% @args(1) axis direction [x,y] %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function m = matrix2DReflection(args)  
    if (length(args) ~= 1 || ~isfloat(args{1}) || length(args{1}) ~= 2)
        disp('2D reflection matrix expects 1 input: vector [x y] of reflection axis');
        m = -1;
        return;
    end
    
    unit = args{1} / norm(args{1});
    x = unit(1);
    y = unit(2);
    m = [ (x*x - y*y)  2*x*y; 
           2*x*y      (y*y - x*x) ];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3D REFLECTION MATRIX          %
% @args(1) axis {'x','y','z'}   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function m = matrix3DReflection(args)  
    if (length(args) ~= 1 || ~ischar(args{1}))
        disp('3D reflection matrix expects 1 input: axis {''x'', ''y'', ''z''}');
        m = -1;
        return;
    end
     
    axis = args{1};
    switch axis
        case 'x'
            m = [ -1 0 0  ;
                   0 1 0  ; 
                   0 0 1 ];
        case 'y'
            m = [ 1  0 0  ;
                  0 -1 0  ; 
                  0  0 1 ];
        case 'z'
            m = [ 1 0  0  ;
                  0 1  0  ; 
                  0 0 -1 ];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2D SCALE MATRIX                    %
% @args(1) scale factor x            %
%      (2) scale factor y (optional) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function m = matrix2DScale(args) 
    if (length(args) < 1 || length(args) > 2)
        disp('2D scale matrix expects 1 OR 2 inputs: scale factor x AND scale factor y (optional for isotropic scale)');
        m = -1;
        return;
    end
    sx = args{1}
    
    if (length(args) == 1)
        sy = sx; %isotropic scale
    else
        sy = args{2}
    end 
    m = [ sx  0  ; 
           0 sy ];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3D SCALE MATRIX                    %
% @args(1) scale factor x            %
%      (2) scale factor y (optional) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function m = matrix3DScale(args) 
    if (length(args) < 1 || length(args) == 2 || length(args) > 3)
        disp('3D scale matrix expects 1 OR 3 inputs: scale factor x AND scale factor y AND scale factor z (y and z optional for isotropic scale)');
        m = -1;
        return;
    end
    sx = args{1};
    
    if length(args) == 1
        sy = sx; %isotropic scale
        sz = sx;
    else 
        sy = args{2};
        sz = args{3};
    end 
    m = [ sx  0  0  ;
           0 sy  0  ; 
           0  0 sz ];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2D TRANSLATION MATRIX              %
% @args(1) translation dx            %
%      (2) translation dy            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function m = matrix2DTranslation(args) 
    if (length(args) ~= 2 )
        disp('2D translation matrix expects 2 inputs: translation factor dx AND dy');
        m = -1;
        return;
    end
    dx = args{1};
    dy = args{2};
    
    m = [ 1 0 dx ;
          0 1 dy ; 
          0 0  1 ];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3D TRANSLATION MATRIX              %
% @args(1) translation dx            %
%      (2) translation dy            %
%      (3) translation dz            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function m = matrix3DTranslation(args) 
    if (length(args) ~= 3 )
        disp('3D translation matrix expects 3 inputs: translation factor dx, dy AND dz');
        m = -1;
        return;
    end
    dx = args{1};
    dy = args{2};
    dz = args{3};
    
    m = [ 1 0 0 dx ;
          0 1 0 dy ; 
          0 0 1 dz ;
          0 0 0  1 ];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3D PROJECTION MATRIX               %
% @args(1) factor       d            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function m = matrix3DPerspectiveProjection(args) 
    if (length(args) ~= 2 )
        disp('3D perspective projection matrix expects 1 inputs: factor d');
        m = -1;
        return;
    end
    d = args{1};
    
    m = [ 1 0  0  0 ;
          0 1  0  0 ; 
          0 0  1  0 ;
          0 0 1/d 1 ];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3D HOUSEHOLDER MATRIX              %
% @args(1) factor       d            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HouseholderMatrix[v_?VectorQ] :=
%     IdentityMatrix[Length[v]]
%       - 2 Transpose[{v}] . {v} / (v.v)
      
function m = matrix3DHouseholder(args) 
    if (length(args) ~= 1 )
        disp('3D householder projection matrix expects 1 inputs: vector [x, y]');
        m = -1;
        return;
    end
    
    unit = args{1} / norm(args{1}); %normalize v
    x = unit(1);
    y = unit(2);
    
    m = [ 1 0  0 ;
          0 1  0 ; 
          0 0  1 ];
      
    m = m - 2 * (unit'*unit); 
end