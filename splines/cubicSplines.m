% @file     cubicSplines.m
% @author   afruehstueck
% @date     07/02/2017

%evaluate arguments and pick corresponding spline function
function [viewer, points] = cubicSplines(x, y, z, closed)
    if ~exist('closed', 'var')
        closed = 0;
    end

    if ~exist('z', 'var') || isempty(z) || ~any(z) %2 dimensional
        cubic2DSpline(x, y, closed);
    else %3 dimensional
        cubic3DSpline(x, y, z, closed);
    end
end

%calculates spline for two-dimensional input
function [] = cubic2DSpline(x, y, closed)    
    abcd_x = evaluate1DSpline(x, closed);
    abcd_y = evaluate1DSpline(y, closed);
    spl = size(abcd_x, 1); %number of splines
    
    u = linspace(0, 1);
    u3 = u.^3;
    u2 = u.^2;
    
    for i = 1:spl %step through segments
        c_x = abcd_x(i,:);
        c_y = abcd_y(i,:);
        
        p_x = c_x(4).*u3 + c_x(3).*u2 + c_x(2).*u + c_x(1);
        p_y = c_y(4).*u3 + c_y(3).*u2 + c_y(2).*u + c_y(1);

        %plot spline segment
        plot(p_x, p_y,'-');
    end
end

%calculates spline for three-dimensional input
function [] = cubic3DSpline(x, y, z, closed)
    if (~exist('closed', 'var'))
        closed = 0;
    end
    
    abcd_x = evaluate1DSpline(x, closed);
    abcd_y = evaluate1DSpline(y, closed);
    abcd_z = evaluate1DSpline(z, closed);
    spl = size(abcd_x, 1); %number of splines
    
    u = linspace(0, 1);
    u3 = u.^3;
    u2 = u.^2;
    
    for i = 1:spl %step through segments
        c_x = abcd_x(i,:);
        c_y = abcd_y(i,:);
        c_z = abcd_z(i,:);
        
        p_x = c_x(4).*u3 + c_x(3).*u2 + c_x(2).*u + c_x(1);
        p_y = c_y(4).*u3 + c_y(3).*u2 + c_y(2).*u + c_y(1);
        p_z = c_z(4).*u3 + c_z(3).*u2 + c_z(2).*u + c_z(1);

        %plot spline segment
        plot3(p_x, p_y, p_z,'-')
    end
end

function abcd = evaluate1DSpline( v, closed ) 
    pts = length(v);
    spl = pts - 1; %number of splines

    M = zeros(pts);
    b = zeros(pts, 1);

%           M       .   D     =         b   
%     [ 2 1        ] [ D(1) ]   [ 3(x(2) - x(1))   ]
%     | 1 4 1      | | D(2) |   | 3(x(3) - x(1))   |
%     |   1 4 1    | |  .   | = |       .          |
%     |     .....  | |  .   |   |       .          |
%     |      1 4 1 | |  .   |   | 3(x(n) - x(n-2)) |
%     [        1 2 ] [ D(n) ]   [ 3(x(n) - x(n-1)) ]
       
    for i = 2:pts-1
        M(i, i-1:i+1) = [1 4 1];
        b(i) = 3*(v(i+1) - v(i-1));
    end

    %different edge conditions for closed or open spline curve
    if closed == 1
        %first row
        M(1, 1:2) = [4 1];
        M(1, pts) = 1;
        b(1) = 3*(v(2) - v(pts));
        %last row
        M(pts, 1) = 1;
        M(pts, pts-1:pts) = [1 4];
        b(pts) = 3*(v(1) - v(pts-1));
    else
        %first row
        M(1, 1:2) = [2 1];
        b(1) = 3*(v(2) - v(1));
        %last row
        M(pts, pts-1:pts) = [1 2];
        b(pts) = 3*(v(pts) - v(pts-1));
    end

    %solve system
    D = M \ b;
    
    %store coefficients
    abcd = zeros(spl, 4);
    for i = 1:spl %step through segments
        a = v(i);
        b = D(i);
        c = 3*(v(i+1) - v(i)) - 2*D(i) - D(i+1);
        d = 2*(v(i) - v(i+1)) + D(i) + D(i+1);
        abcd(i, :) = [a b c d];
    end   
    
    if closed == 1 %add another segment for closing the curve
        a = v(pts);
        b = D(pts);
        c = 3*(v(1) - v(pts)) - 2*D(pts) - D(1);
        d = 2*(v(pts) - v(1)) + D(pts) + D(1);
        
        abcd(spl+1, :) = [a b c d];
    end
end
