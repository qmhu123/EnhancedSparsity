% Value of Hermite polynomial (prob) at point x, which can be either
% a scalar or a column vector.
%
% Author: Xiu Yang
% Update: July 5th, 2012

function [val] = hermite(x, order)

val = ones(length(x), order+1);

if order > 0
  val(:,2) = x;

  for i = 2:order
    val(:,i+1) = (x.*val(:,i) - (i-1)*val(:,i-1));
  end
end
