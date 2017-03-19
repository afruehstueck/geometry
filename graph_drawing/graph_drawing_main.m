% @file     graph_drawing_main.m
% @author   afruehstueck
% @date     18/03/2017

function graph_drawing_main
    close all;
    clc;
    clear;

    %define labels strings to label vertices in graph
    labels = cell(20, 1);
    for i= 1:20
        labels{i} = strcat('P', num2str(i));
    end
    %% figure 1  
    graphs_fig = figure('Name', 'Spectral Graph Drawing examples', 'NumberTitle', 'off', 'Position', get(0, 'ScreenSize'));

    % define the default values for text
    set(graphs_fig, 'DefaultTextFontName', 'consolas', ...  % So text lines up
                    'DefaultTextFontSize', 7 )

    %4 connected vertices           
    A = [0 1 1 0; 
         1 0 0 1; 
         1 0 0 1; 
         0 1 1 0];
    subplot(2, 3, 1); 
    hold on;
    drawSpectralGraph(A);
    plotMatrixInFigure(graphs_fig, A);

    A = [0 1 1 0 0; 
         1 0 1 1 1; 
         1 1 0 1 0; 
         0 1 1 0 1; 
         0 1 0 1 0];
    subplot(2, 3, 2); 
    hold on;
    drawSpectralGraph(A);
    plotMatrixInFigure(graphs_fig, A);

    %12-vertex ring
    A = diag(ones(1, 11), 1);
    A = A + A';
    A(1, 12) = 1; 
    A(12, 1) = 1;
    subplot(2, 3, 3); 
    drawSpectralGraph(A);
    plotMatrixInFigure(graphs_fig, A);

    %random 5x5
    A = randomAdjacencyMatrix(5);
    subplot(2, 3, 4); 
    drawSpectralGraph(A);
    plotMatrixInFigure(graphs_fig, A);

    %random 8x8
    A = randomAdjacencyMatrix(8);
    subplot(2, 3, 5); 
    drawSpectralGraph(A);
    plotMatrixInFigure(graphs_fig, A);

    %random 12x12
    A = randomAdjacencyMatrix(12);
    subplot(2, 3, 6); 
    drawSpectralGraph(A);
    plotMatrixInFigure(graphs_fig, A);

    %% figure 2
    graphs_fig2 = figure('Name', '20 random 5x5 graphs', 'NumberTitle', 'off', 'Position', get(0, 'ScreenSize'));
    % define the default values for text
    set(graphs_fig2, 'DefaultTextFontName', 'consolas', ...  % So text lines up
                     'DefaultTextHorizontalAlignment', 'left', ...
                     'DefaultTextVerticalAlignment', 'bottom', ...
                     'DefaultTextClipping', 'on', ...
                     'DefaultTextFontSize', 7 )

    %plot 20 randomized 5x5 graphs
    for y = 1:4
        for x = 1:5
            subplot(4, 5, (y - 1)*5 + x);
            A = randomAdjacencyMatrix(5);
            drawSpectralGraph(A);
            plotMatrixInFigure(graphs_fig2, A);
        end
    end

    %% helper functions

    %generate random adjacency matrix
    function [A] = randomAdjacencyMatrix(N)    
        %generate NxN matrix of random 0s and 1s
        G = round(rand(N));
        %take only upper triangular part and replicate on lower triangular
        %while leaving the main diagonal as zeros (no loops)
        randA = triu(G, 1) + triu(G, 1)';

        %check if each vertex has at least one edge
        if any(sum(randA) < 2)
            disp('<2 edges for at least one vertex... generate new matrix');
            A = randomAdjacencyMatrix(N);  
        else
            A = randA;
        end
    end

    % spectral graph drawing of graph with adjacency matrix A
    % following Jean H. Gallier, CIS 515, Spectral Graph Drawing
    function [V, Dvec, L] = drawSpectralGraph(A)    
        %degree matrix (#edges for each vertex)
        D = diag(sum(A));
        %laplacian matrix (degree matrix - adjacency matrix)
        L = D - A;
        %find eigenvectors and eigenvalues of laplacian matrix
        [eVecs, ~] = eig(L);

        hold on
        %fix limits to make all limits equal and make room for matrix in plot
        xlim([-1.05, 1.05])
        ylim([-1.05, 1.05])
        %draw vertices
        gplot(A, eVecs(:, [3 2]), 'o')
        %draw edges
        gplot(A, eVecs(:, [3 2]))
        
        text(eVecs(:, 3) .* 1.1, eVecs(:, 2) .* 1.1, labels(1:size(A, 1)));
    end

    % function inspired by https://www.mathworks.com/matlabcentral/answers/95255-how-do-i-display-a-matrix-in-a-figure-window
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
end