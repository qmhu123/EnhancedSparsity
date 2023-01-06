% Value of the derivative of (probabilistic) Hermite polynomial at point x, 
% which can be either a scalar or a column vector.
%
% Author: Xiu Yang
% Update: July 5th, 2012

function [dval] = dhermite(x, order)

val = hermite(x, order-1);
dval = zeros(length(x), order+1);

if order > 0
  dval(:,2) = ones(length(x),1);

  for i = 2:order
    dval(:,i+1) = i*val(:,i);
  end
end
