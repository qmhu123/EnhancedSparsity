function x = f_prox_newton(y, sigma, lambda, iter, x0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function finds the root for     
%   exp(-x^2/sigma^2) sign(x) + lambda (x - y) = 0
% and solves min ERF(x) + lambda/2||x-y||^2
% Inputs: y, lambda, sigma
%         iter: the maximum number of iterations
%         x0: the initial value
% Output: x
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    yabs = abs(y);    % we only consider the absolution values and keep the sign.
    index  = yabs > 1/lambda;  % Find out the zeros 
    x = zeros(size(y)); 

    lambda = sigma*lambda;   % for the nonzero values, we can change x to be x/sigma
    yabs = yabs/lambda;
    x(index) = abs(x0(index)/lambda);  % the initial of abs(x), we only consider the absolute values
    for i = 1:iter        
        fx  = exp(-x(index).^2) + lambda * (x(index) - yabs(index));   % f(x)
%        max(abs(fx))
        if max(abs(fx)) < 1e-10
            break
        end
        fxx = exp(-x(index).^2).*(-2 * x(index)) + lambda;  % f'(x) 
%        min(abs(fxx))
        if i < 20
            x(index)   = x(index) - fx./max(fxx,0.1);  % Newton step     
        else
            x(index)   = x(index) - fx./fxx;  % Newton step
        end
%        plot(fx)
%        [x fx fxx]
    end
    x = sigma * x .* sign(y);
%    i
end