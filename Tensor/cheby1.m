% Value of Chebyshev polynomial of the first kind at point x, which can be 
% either a scalar or a column vector.
%
% Author: Xiu Yang
% Update: Feb 28th, 2017

function [val] = cheby1(x, order)

val = ones(length(x), order+1);

if order > 0
  val(:,2) = x;

  for i = 2:order
    val(:,i+1) = 2.0*x.*val(:,i) - 1.0*val(:,i-1);
  end
end
