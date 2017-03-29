% @file     centroid.m
% @author   afruehstueck
% @date     24/03/2017

% calculate coordinates of centroid of polygon from polygon vertices vX/
% according to https://en.wikipedia.org/wiki/Centroid#Centroid_of_polygon
function [cX, cY] = centroid(vX, vY)
    vXp1 = vX([2:end 1]); %x coordinates with shifted indices (plus 1)
    vYp1 = vY([2:end 1]); %y coordinates with shifted indices (plus 1)

    %term used multiple times in calculation
    mult = vX.*vYp1 - vXp1.*vY;

    A = sum(mult) / 2; %signed area of the polygon

    cX = sum((vX + vXp1) .* mult) / (6 * A);
    cY = sum((vY + vYp1) .* mult) / (6 * A);
end