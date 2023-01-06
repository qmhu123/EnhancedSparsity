% Value of the derivative of Chebyshev polynomial of the first kind at point x,
% which can be either a scalar or a column vector.
%
% Author: Xiu Yang
% Update: Feb 28th, 2017

function [dval] = dcheby1_norm(x, order)

val = cheby2(x, order-1);
dval = zeros(length(x), order+1);
dval(:,2) = ones(length(x),1);

if order > 0
  for i = 2:order
    dval(:,i+1) = i*val(:,i);
  end

  dval(:,2:end) = dval(:,2:end)*sqrt(2.0);
end
