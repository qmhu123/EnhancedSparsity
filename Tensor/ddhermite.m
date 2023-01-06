% Value of the second order derivative of (probabilistic) Hermite polynomial at
% point x, which can be either a scalar or a column vector.
%
% Author: Xiu Yang
% Update: July 22th, 2015

function [ddval] = ddhermite(x, order)

val = hermite(x, order-1);
ddval = zeros(length(x), order+1);
ddval(:,3) = 2*ones(length(x),1);

for i = 3:order
  ddval(:,i+1) = i*(i-1)*val(:,i-1);
end
