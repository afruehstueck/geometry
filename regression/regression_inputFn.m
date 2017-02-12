% @file     regression_inputFn.m
% @author   afruehstueck
% @date     06/02/2017

%takes input values x and y and does least squares fitting on the data
function [] = regression_inputFn(x, y, descr)
%     if ~exist('dimX', 'var')
%         dimX = [-5, 5];
%     end
% 
%     if ~exist( 'dimY', 'var')
%         dimY = [-5, 5];
%     end

    scr = get(0, 'ScreenSize');  
    fig = figure('Name', 'Least-Squares Fitting', 'NumberTitle', 'off', 'Position', [scr(3)/2 50 scr(3)/3 scr(3)/3]);
    if exist( 'descr', 'var')
        title(descr);
    end
    hold on;
    axis auto;   

    num = length(x);  
    
    plot(x, y, 'o', 'DisplayName', 'data with noise');
    legend('-DynamicLegend');
    legend('show')
    drawnow;
    
    prompt = 'Input degree of output polynomial [press ''Enter'' to exit]: ';
    while ishandle(fig)    
        str = input(prompt, 's');
        if isempty(str), break, end
        
        deg = str2num(str);
        
        A = ones(num, deg + 1);

        for col=1:deg 
            A(:, col) = x'.^(deg + 1 - col);
        end

        b = y';
        coeff = A \ b

        cla(fig);
        plot(x, y, 'o');
        plot(x, polyval(coeff, x) );
        
        
%         coeff_reg = ( A' * A + lambda * eye ( size(A, 2) ) ) \ ( A' * b )
%         plot(x, polyval(coeff_reg, x), 'b');
%         drawnow;

        ln_lambda = -2.3;
        lambda = exp(ln_lambda);
        An = [ A; (sqrt(lambda) * eye ( size(A, 2) )) ]
        bn = [ b; zeros(size(A, 2), 1) ]
        coeff_reg = An \ bn
        plot(x, polyval(coeff_reg, x));
        
        ln_lambda = 0;
        lambda = exp(ln_lambda);
        An = [ A; (sqrt(lambda) * eye ( size(A, 2) )) ]
        bn = [ b; zeros(size(A, 2), 1) ]
        coeff_reg = An \ bn
        plot(x, polyval(coeff_reg, x));
        
        legend('data with noise', 'no regularization', 'regularization with ln \lambda = -2.3', 'regularization with ln \lambda = 0');
        %drawnow;
    end
end