function [x,result] = mReWexp_constrained_ADMM_Proximal(A,b,pm)

result = pm; 

rho = 1000;
sigma = 1;
reltol  = 1e-6; 

if isfield(pm,'rho'); rho = pm.rho; end
if isfield(pm,'sigma'); sigma = pm.sigma; end
if isfield(pm,'reltol'); reltol = pm.reltol; end
% if isfield(pm,'xg'); xg = pm.xg; end
% L1 solution as initial point(if initial is not given)
if isfield(pm,'xr')
    xr = pm.xr;
else
    xr = CS_L1_uncon_ADMM(A,b,pm);
%     xr = mL1_constrained_ADMM_Proximal(A,b,pm);
end

obj = @(t) sum(erf(abs(t)/sigma));

% [M,N] = size(A);
[~,N] = size(A);

AAt = A'*((A*A')\b);
B = eye(N) - A'*((A*A')\A);
iter = 15*N;
eps = 1e-9;
x = xr;
y = zeros(N,1);
u = zeros(N,1);


for i = 1:iter
    xold = x;
    % x update
    x = B*(y - u/rho) + AAt;
    % y update
    %y = shrink(x+u/rho,1/rho);
    y = f_prox_newton(x+u/rho, sigma, rho, 100, x+u/rho);
    
    % u update
    u = u + rho*(x-y);
    if norm(x-xold)/(norm(xold)+eps) < reltol
        break;
    end
    
end


% i is the number of iterations
result.i = i;
% fx is the objective function
fx = obj(x);
result.fx = fx;
% % fxg is the objective function at xg(the true solution)
% fxg = sum(erf(abs(xg)/sigma));
% result.fxg = fxg;

%% Evaluation 
% result.error = norm(x-xg)/norm(xg);

% 
% if result.error < pm.restol % success
%       result.rate = 1;
%     elseif fx+ eps < fxg % model failure
%          result.rate = -1;
%     else    % algorithm failure
%         result.rate = -2;
% end


end


function x = f_prox_newton(y, sigma, lambda, iter, x0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function finds the root for     
%   exp(-x^2/sigma^2) sign(x) + lambda (x - y) = 0
% Inputs: y, lambda, sigma
%         iter: the maximum number of iterations
%         x0: the initial value
% Output: x
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sizey = size(y);
    y = y(:);
    yabs = abs(y);    % we only consider the absolution values and keep the sign.
    index  = yabs > 1/lambda;  % Find out the zeros 
    index2 = ~index;
    x = zeros(size(y)); 

    lambda = sigma*lambda;   % for the nonzero values, we can change x to be x/sigma
    yabs = yabs/sigma;
    x(index) = abs(x0(index)/sigma);  % the initial of abs(x), we only consider the absolute values
    for i = 1:iter        
        fx  = exp(-x.^2) + lambda * (x - yabs);   % f(x)
%        max(abs(fx))
        if max(abs(fx(index))) < 1e-8
            break
        end
        fxx = (-2 * x) .* exp(-x.^2) + lambda;  % f'(x) 
%        min(abs(fxx))
        if i < 20
            x   = x - fx./max(fxx,0.1);  % Newton step     
%            x   = x - fx./fxx;  % Newton step
        else
            x   = x - fx./fxx;  % Newton step
        end
        
%        plot(fx)
%        [x fx fxx]
    end
    x = sigma * x .* sign(y);
    x(index2) = 0;   
    x = reshape(x,sizey);
end