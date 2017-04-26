% @file     regression_inputClick.m
% @author   afruehstueck
% @date     02/02/2017
%
%create a plot and ask user for input via mouse click
%fit a polynomial of degree k-1 through the k points selected by the user
%also, fit polynomials with regularization

function [] = regression_inputClick(allow_input, show_title)    
    if (~exist('show_title', 'var'))
        show_title = 0;
    end
 
    if (~exist('allow_input', 'var')) 
        allow_input = 1;
    end
    
    dimX = [-5, 5];
    dimY = [-5, 5];

    scr = get(0, 'ScreenSize');  
    fig = figure('Name', 'Polynomial Regression Picker', 'NumberTitle', 'off', 'Position', [scr(3)/2 50 scr(3)/3 scr(3)/3]);
    hold on;
    axis([dimX, dimY]);    
  
    x = zeros(1,1);
    y = zeros(1,1);
    r = linspace(-5, 5);
    num = 1;
    
    % Give instructions
    disp('Left mouse button picks points.')
    disp('Right mouse button picks last point.')
    % Stores integer value for which button is pressed
    % 1 for left, 2 for middle, 3 for right
    
    button = 1; 
    try
        while button == 1
            [x(num), y(num), button] = ginput(1);

            A = ones(num);

            for col=1:num-1 
                A(:, col) = x'.^(num - col);
            end

            b = y';
            %X = linsolve(A, b);
            coeff = A \ b;

            cla(fig);

            % DISPLAY POLYNOMIAL AS FIGURE TITLE
            if show_title == 1
            poly_title = '';
            for row=num:-1:1
                c = num2str(abs(coeff(row)), 3);
                if sign(coeff(row)) == 1
                    c = [' + ' c];
                elseif sign(coeff(row)) == -1
                    c = [' - ' c];
                end
                e = '';
                if row > 2 
                    e = [ 'x^' num2str(row-1) ];
                elseif row == 2
                    e = 'x';
                end
                poly_title = [poly_title c e];
            end
            title(poly_title)
            end

            plot(x, y, 'o');
            plot(r, polyval(coeff, r), 'LineWidth', 1.5);

            ln_lambda = -2.3;
            lambda = exp(ln_lambda);
            coeff_reg = ( A' * A + lambda * eye ( size(A, 2) ) ) \ ( A' * b )
            plot(r, polyval(coeff_reg, r));

            ln_lambda = 0;
            lambda = exp(ln_lambda);
            coeff_reg = ( A' * A + lambda * eye ( size(A, 2) ) ) \ ( A' * b )
            plot(r, polyval(coeff_reg, r));

            legend('data points', 'no regularization', 'regularization with ln \lambda = -2.3' ,'regularization with ln \lambda = 0');
            drawnow;

            num = num + 1;
        end

        if allow_input == 1
            while ishandle(fig)    
                prompt = 'Input degree of output polynomial [press ''Enter'' to exit]: ';
                str = input(prompt, 's');
                if isempty(str) 
                    fprintf('Input terminated.');
                    break
                end

                deg = str2num(str);

                num = length(x);
                A = ones(num, deg + 1);

                for col=1:deg 
                    A(:, col) = x'.^(deg + 1 - col);
                end

                b = y';
                %X = linsolve(A, b);
                coeff = A \ b

                cla(fig);

                plot(x, y, 'o');
                plot(r, polyval(coeff, r), 'LineWidth', 1.5);

                ln_lambda = -2.3;
                lambda = exp(ln_lambda);
                coeff_reg = ( A' * A + lambda * eye ( size(A, 2) ) ) \ ( A' * b )
                plot(r, polyval(coeff_reg, r));

                ln_lambda = 0;
                lambda = exp(ln_lambda);
                coeff_reg = ( A' * A + lambda * eye ( size(A, 2) ) ) \ ( A' * b )
                plot(r, polyval(coeff_reg, r));
            end
        end
    catch
        allow_input = 0;
        button = 0;
        fprintf('Figure closed.\n');
    end
end