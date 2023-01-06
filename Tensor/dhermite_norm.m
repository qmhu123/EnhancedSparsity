% Value of the derivative of (probabilistic) Hermite polynomial at point x, 
% which can be either a scalar or a column vector.
%
% Author: Xiu Yang
% Update: July 5th, 2012

function [dval] = dhermite_norm(x, order)

val = hermite(x, order-1);
dval = zeros(length(x), order+1);
dval(:,2) = ones(length(x),1);

for i = 2:order
  dval(:,i+1) = i*val(:,i);
end

for i = 2:order
  dval(:,i+1) = dval(:,i+1)/sqrt(factorial(i));
end
