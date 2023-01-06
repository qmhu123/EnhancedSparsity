function [x,output] = CS_uncon_IRL14Lp(A,b,pm)
%min_x .5||Ax-b||^2 + lambda \sum w_i |x_i|

%Input: dictionary A, data b, parameters set pm
%       pm.lambda: regularization paramter
%       pm.delta: penalty parameter for ADMM, default value: 10*lambda
%       pm.maxoit: max outer iterations, default value: 10
%       pm.maxit: max inner iterations: default value: 5000
%       pm.tol: outer tolerace, default value: 1e-3
%       pm.abstol: abs tolerance for ADMM: default value: 1e-7
%       pm.reltol: rel tolerance for ADMM: default value: 1e-5
%Output: computed coefficients x
%       output.relerr: relative error of x_n and x_{n-1}
%       output.obj: objective function of x_n:
%       obj(x) = lambda |x|_1+0.5|Ax-b|^2
%       output.res: residual of x_n: norm(Ax-b)/norm(b)
%       output.err: error to the ground-truth: norm(x-xg)/norm(xg)
%       output.relerr: error between two neighboring iterates
%       output.time: computational time

%% DCA for outer interation
[M,N] = size(A); start_time = tic;


%% parameters
lambda = 1e-5; delta = 100*lambda;
maxit = 5*N; %maxoit = 10;
x0 = zeros(N,1); xg = x0;
reltol = 1e-6; abstol = 1e-6;
eps = 1e-16; %sigma = 1;
epsit = 1;
mu = 0.2;


if isfield(pm,'delta'); delta = pm.delta; end
% if isfield(pm,'sigma'); sigma = pm.sigma; end
if isfield(pm,'lambda'); lambda = pm.lambda; end
if isfield(pm,'maxit'); maxit = pm.maxit; end
if isfield(pm,'x0'); x0 = pm.x0; end
if isfield(pm,'xg'); xg = pm.xg; end
if isfield(pm,'tol'); tol = pm.tol; end
if isfield(pm,'abstol'); abstol = pm.abstol; end
if isfield(pm,'reltol'); reltol = pm.reltol; end


%% pre-computing/initialize
AAt = A*A';
L = chol( speye(M) + 1/delta*AAt, 'lower' );
L = sparse(L);
U = sparse(L');

x = x0;
Atb = A'*b;
y = zeros(N,1); u = y;

obj = @(x,w) .5*norm(A*x-b)^2+ lambda*(w'*abs(x));
w = ones(N,1);

for it = 1:20
    
    wp = w;
    
    for inner = 1:maxit/20
        
        %update x
        xold = x;
        x =shrink(y-u, lambda/delta*w);
        
        
        %update y
        yold = y;
        rhs = Atb  + delta*(x+u);
        y = rhs/delta - (A'*(U\(L\(A*rhs))))/delta^2;
        
        %update u
        u = u + x-y;

        
        % stop conditions & outputs
        relerr = norm(x-y)/max([norm(x),norm(y),eps]);
        residual    = norm(A*x - b);
        
        if  residual < abstol || relerr < reltol
            break;
        end
    end
    
    output.relerr(it) = relerr;
    output.obj(it) = obj(x,w);
    output.time(it) = toc(start_time);
    output.res(it) = residual;
    output.err(it) = norm(x-xg)/norm(xg);
    output.L0(it) = nnz(x);
    
    w =  0.5 ./ sqrt(abs(x)+epsit);
    epsit = epsit*mu;
    if norm(w-wp)/norm(wp) < reltol
        break;
    end
    
end



end

%% Shrinkage
function z = shrink(x, r)
z = sign(x).*max(abs(x)-r,0);
end

