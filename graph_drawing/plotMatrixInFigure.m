% @file     plotMatrixInFigure.m
% @author   afruehstueck
% @date     20/03/2017
%
% displays a matrix of numbers in a figure window
% inspired by https://www.mathworks.com/matlabcentral/answers/95255-how-do-i-display-a-matrix-in-a-figure-window

function [] = plotMatrixInFigure(figure, A)    
    [m,n] = size(A);
    %get axes of current figure
    xLim = figure.CurrentAxes.XLim;
    yLim = figure.CurrentAxes.YLim;

    %create temporary character to measure its width and height
    tmp = text(.5, .5, '0');
    chsz = get(tmp, 'Extent');
    delete(tmp)
    dy = chsz(4);       % height of single character
    dx = 2.5 * chsz(3); % width of single character
    mx = xLim(1) + dx;  % location of first column

    for mi = 1:n % step through columns
        my = yLim(1) + 0.5 * dy;      
        for mj = m:-1:1 % step through rows (inverse)
          my = my + abs(dy);  % location of row

          color = [0 0 0];
          %print diagonal entries in grey
          if mi == mj
              color = [0.7 0.5 0.5];
          end
          %print single digit output
          text(mx, my, sprintf('%1.0f', A(mj, mi)), 'Color', color);
        end
        mx = mx + dx;
    end
end