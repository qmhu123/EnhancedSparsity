% Value of Charlier polynomial at point x, which can be either a scalar or a
% column vector.
%
% Author: Xiu Yang
% Update: May 3rd, 2018

function [val] = charlierF(x, order, lambda)

if order == -1
  val = zeros(length(x), order+1);
end

if order == 0
  val = ones(length(x), order+1);
end

if order > 0
  val = x.*charlierF(x-1, order-1, lambda) - lambda*charlierF(x, order-1, lambda);
end

