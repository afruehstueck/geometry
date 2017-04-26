% @file     graph_drawing_main.m
% @author   afruehstueck
% @date     18/03/2017
%
% generates several random adjacency matrices and does spectral graph
% drawing by analyzing the eigenvectors of the graph's laplacian matrix
% following tutorial by Jean H. Gallier, CIS 515, Spectral Graph Drawing

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
    drawSpectralGraph(A, labels);
    plotMatrixInFigure(graphs_fig, A);

    A = [0 1 1 0 0; 
         1 0 1 1 1; 
         1 1 0 1 0; 
         0 1 1 0 1; 
         0 1 0 1 0];
    subplot(2, 3, 2); 
    hold on;
    drawSpectralGraph(A, labels);
    plotMatrixInFigure(graphs_fig, A);

    %12-vertex ring
    A = diag(ones(1, 11), 1);
    A = A + A';
    A(1, 12) = 1; 
    A(12, 1) = 1;
    subplot(2, 3, 3); 
    drawSpectralGraph(A, labels);
    plotMatrixInFigure(graphs_fig, A);

    %random 5x5
    A = randomAdjacencyMatrix(5);
    subplot(2, 3, 4); 
    drawSpectralGraph(A, labels);
    plotMatrixInFigure(graphs_fig, A);

    %random 8x8
    A = randomAdjacencyMatrix(8);
    subplot(2, 3, 5); 
    drawSpectralGraph(A, labels);
    plotMatrixInFigure(graphs_fig, A);

    %random 12x12
    A = randomAdjacencyMatrix(12);
    subplot(2, 3, 6); 
    drawSpectralGraph(A, labels);
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
            drawSpectralGraph(A, labels);
            plotMatrixInFigure(graphs_fig2, A);
        end
    end
end