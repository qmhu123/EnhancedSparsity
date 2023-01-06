% Value of Laguerre polynomial at point x, which can be either
% a scalar or a column vector.
%
% Author: Xiu Yang
% Update: August 5th, 2012

function [val] = laguerre_norm(x, order)

val = ones(length(x), order+1);
val(:,2) = 1-x;

for i = 2:order
  val(:,i+1) = ((2*i-1-x).*val(:,i) - (i-1)*val(:,i-1))/i;
end
