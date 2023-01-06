% Value of Legendre polynomial at point x, which can be either a scalar or a
% column vector.
%
% Author: Xiu Yang
% Update: July 5th, 2012

function [val] = legendre(x, order)

val = ones(length(x), order+1);
val(:,2) = x;

for i = 2:order
  val(:,i+1) = ((2.0*i-1)*x.*val(:,i) - (i-1.0)*val(:,i-1))/i;
end
