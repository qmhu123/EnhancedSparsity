% Value of n-th order polynomial at point x, which can be either a scalar or a
% column vector.
%
% Author: Xiu Yang
% Update: July 5th, 2012

function [val] = taylor(x, order)

val = ones(length(x), order+1);
val(:,2) = x;

for i = 1:order
  val(:,i+1) = x.^i;
end
