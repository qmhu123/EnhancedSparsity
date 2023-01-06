%--------------------------------------------------------------------------
% Difference of Convex Functions Methods for 1-D Compressed Sensing Problem.
%
% Solves
%           min  P_a(x)
%           s.t. Ax = b
%
% Reference: "Minimization of Transformed L_1 Penalty: Theory,
%             Difference of Convex Function Algorithm,
%             and Robust Application in Compressed Sensing."
%             Shuai Zhang, Jack Xin
%             IEEE transactions on Information Theory, 2014
% Available at:
%             http://arxiv.org/abs/1411.5735
%             ftp://ftp.math.ucla.edu/pub/camreport/cam14-86.pdf
%
% Author: Shuai Zhang
% Date: Feb 07. 2015
%--------------------------------------------------------------------------

% min_x  .5||Ax-b||^2 + lambda(sum  (a+1)|x_i| / (a+|x_i|) )
% here g(x) =  .5||Ax-b||^2 + C||x||^2 + lambda* (a+1)/a * ||x||_1
%      h(x) = lambda*{(a+1)/a *||x||_1 - sum (a+1)|x_i|/(a+|x_i|)} + C||x||^2
%      F(x) = g(x) - h(x)
%Input: dictionary A, data b, parameters option and initial value x0
% recommended value
% x0: initial value, better to chosed as zeros vector, see reference
% pm.C       1.0e-9
% pm.del     1.0e-4
% pm.maxoit  1000
% pm.lambda     1.0e-5
% pm.tol     1.0e-5
% pm.a       1
%Output: computed coefficients x


function [x,output] = CS_TL1_DCA(A,b,pm)
[M,N] = size(A); start_time = tic;
x = zeros(N,1);
if ~exist('pm','var')
    pm = [];
end
C = 1.0e-9;
maxoit = 20;
lambda = 1.0e-6;
tol = 1.0e-5;
del = lambda*10;
a = 1;
if isfield(pm,'C'),        C = pm.C;             end
if isfield(pm,'del'),      del = pm.del;         end
if isfield(pm,'maxoit'),   maxoit = pm.maxoit;   end
if isfield(pm,'lambda'),      lambda = pm.lambda;         end
if isfield(pm,'tol'),      tol = pm.tol;         end
if isfield(pm,'a'),        a = pm.a;             end
if isfield(pm,'x0'),       x = pm.x0;            end

p_init = zeros(N,1);
%precompute inverse to be used in inner problem
L = chol( speye(M) + 1/(2*C+del)*(A*A'), 'lower' ); % cholesky factorization
L = sparse(L);
U = sparse(L');
Atb = A'*b;


obj = @(x) .5*norm(A*x-b)^2+ lambda*sum(((a+1)*abs(x))./(a+abs(x)));

for it = 1:maxoit
    if norm(x) < eps
        f = -Atb;
    else
        V = lambda*(a+1)/a * sign(x) - lambda*sign(x)*(a+1)./( a+abs(x) ) ...
            + lambda*(a+1)*x./( a+abs(x) ).^2 + 2*C*x;
        f = -Atb - V;
    end
    [y,p_init] = CS_TL1_DCA_sub(A,L,U,f,pm,x,p_init);

    
    % stop conditions & outputs
    relerr = norm(x-y)/max([norm(x),norm(y),eps]);
    residual    = norm(A*x - b);
    
    output.relerr(it) = relerr;
    output.obj(it) = obj(x);
    output.time(it) = toc(start_time);
    output.res(it) = residual;
    %     output.err(it) = norm(x-xg)/norm(xg);
    output.L0(it) = nnz(x);
    
    x = y;
    
    if relerr < tol
        disp(['tolerance met after ' num2str(it) ' iterations with a = ' num2str(a)]);
        break;
    end
    
end

if it == maxoit, fprintf('Unconstrained DCATL1 fails ...... \n'); end

end