% Value of normalized Chebyshev polynomial of the second kind at point x, which
% can be either a scalar or a column vector.
%
% Author: Xiu Yang
% Update: Feb 28th, 2017

function [val] = cheby2_norm(x, order)

val = ones(length(x), order+1);
val(:,2) = x;

for i = 2:order
  val(:,i+1) = 2.0*x.*val(:,i) - 1.0*val(:,i-1);
end

val = val/sqrt(pi/2.0);
